/*
 * rtl-sdr, turns your Realtek RTL2832 based DVB dongle into a SDR receiver
 * Copyright (C) 2012 by Steve Markgraf <steve@steve-m.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* Modified for 2-Frequency-Operation by DC9ST for TDOA operation
   samples at the first frequency, switches to the second frequency and then back to the first again
   http://www.panoradio-sdr.de/tdoa-transmitter-localization-with-rtl-sdrs/
   
   main changes:
   samples 3n IQ samples: first n at frequency f, then n at freuqency h and then again n at frequency f
*/

#include <errno.h>
#include <signal.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _WIN32
#include <unistd.h>
#else
#include <windows.h>
#include <io.h>
#include <fcntl.h>
#include "getopt/getopt.h"
#endif

#include "rtl-sdr.h"
#include "convenience/convenience.h"

#define DEFAULT_SAMPLE_RATE		2048000
#define DEFAULT_BUF_LENGTH		(16 * 16384)
#define MINIMAL_BUF_LENGTH		512
#define MAXIMAL_BUF_LENGTH		(256 * 16384)

static int do_exit = 0;
static uint32_t bytes_to_read_per_freq = 2000000 * 2;
static uint32_t bytes_to_read = 6000000 * 2; // sum of i and q samples
static rtlsdr_dev_t *dev = NULL;
static int gain1_required = 0;
static int gain2_required = 0;
uint32_t frequency1 = 100000000;
uint32_t frequency2 = 100000000;


void usage(void)
{
	fprintf(stderr,
		"\n"
		"rtl_sdr, an I/Q recorder for RTL2832 based DVB-T receivers\n"
		"2-Frequency-Mode for TDOA by DC9ST, 2017\n\n"
		"receives 3xn IQ samples:\n"
		"first n at frequency 1, then n at frequency 2, then n at frequency 1 again\n\n"
		"Usage:\t-f frequency_to_tune_to frequency 1/reference [Hz]\n"
		"\t-h frequency_to_tune_to frequency 2/measure [Hz]\n"
		"\t-1 gain1 (gain for frequency 1/reference) [REQUIRED]\n"
		"\t-2 gain2 (gain for frequency 2/target) [REQUIRED]\n"
		"\t[-s samplerate (default: 2048000 Hz)]\n"
		"\t[-d device_index (default: 0)]\n"
		"\t[-p ppm_error (default: 0)]\n"
		"\t[-b output_block_size (default: 16 * 16384)]\n"
		"\t[-n number n of IQ samples to read per frequency (total length = 3x specified)(default: 2e6)]\n"
		"\t[-S force sync output (default: async)]\n"
		"\tfilename (a '-' dumps samples to stdout)\n\n");
	exit(1);
}

#ifdef _WIN32
BOOL WINAPI
sighandler(int signum)
{
	if (CTRL_C_EVENT == signum) {
		fprintf(stderr, "Signal caught, exiting!\n");
		do_exit = 1;
		rtlsdr_cancel_async(dev);
		return TRUE;
	}
	return FALSE;
}
#else
static void sighandler(int signum)
{
	fprintf(stderr, "Signal caught, exiting!\n");
	do_exit = 1;
	rtlsdr_cancel_async(dev);
}
#endif

static void rtlsdr_callback(unsigned char *buf, uint32_t len, void *ctx)
{
	static uint32_t samples_collected = 0;
	
	if (ctx) {
		if (do_exit)
			return;

		if ((bytes_to_read > 0) && (bytes_to_read < len)) {
			len = bytes_to_read;
			do_exit = 1;
			rtlsdr_cancel_async(dev);
		}

		if (bytes_to_read > 0) {
			uint32_t block_start = samples_collected;
			uint32_t block_end = samples_collected + len;
			
			/* Switch to frequency 2 at start of block 2 */
			if (block_start < bytes_to_read_per_freq && block_end > bytes_to_read_per_freq) {
				rtlsdr_set_agc_mode(dev, 0);
				verbose_gain_set(dev, nearest_gain(dev, gain2_required));
				verbose_set_frequency(dev, frequency2);
				verbose_gain_set(dev, nearest_gain(dev, gain2_required));
			}
			
			/* Switch back to frequency 1 at start of block 3 */
			if (block_start < bytes_to_read_per_freq*2 && block_end > bytes_to_read_per_freq*2) {
				rtlsdr_set_agc_mode(dev, 0);
				verbose_gain_set(dev, nearest_gain(dev, gain1_required));
				verbose_set_frequency(dev, frequency1);
				verbose_gain_set(dev, nearest_gain(dev, gain1_required));
			}

			samples_collected += len;
			bytes_to_read -= len;
		}

		if (fwrite(buf, 1, len, (FILE*)ctx) != len) {
			fprintf(stderr, "Short write, samples lost, exiting!\n");
			rtlsdr_cancel_async(dev);
		}
	}
}

int main(int argc, char **argv)
{
#ifndef _WIN32
	struct sigaction sigact;
#endif
	char *filename = NULL;
	int n_read;
	int r, opt;
	int gain1 = 0;  // Gain for frequency 1 (reference)
	int gain2 = 0;  // Gain for frequency 2 (target)
	int ppm_error = 0;
	int sync_mode = 0;
	FILE *file;
	uint8_t *buffer;
	int dev_index = 0;
	int dev_given = 0;
	uint32_t samp_rate = DEFAULT_SAMPLE_RATE;
	uint32_t out_block_size = DEFAULT_BUF_LENGTH;

	while ((opt = getopt(argc, argv, "d:f:h:1:2:s:b:n:p:S")) != -1) {
		switch (opt) {
		case 'd':
			dev_index = verbose_device_search(optarg);
			dev_given = 1;
			break;
		case 'f':
			frequency1 = (uint32_t)atofs(optarg);
			break;
		case 'h':
			frequency2 = (uint32_t)atofs(optarg);
			break;
		case '1':
			gain1 = (int)(atof(optarg) * 10); /* tenths of a dB for frequency 1 */
			break;
		case '2':
			gain2 = (int)(atof(optarg) * 10); /* tenths of a dB for frequency 2 */
			break;
		case 's':
			samp_rate = (uint32_t)atofs(optarg);
			break;
		case 'p':
			ppm_error = atoi(optarg);
			break;
		case 'b':
			out_block_size = (uint32_t)atof(optarg);
			break;
		case 'n':
			bytes_to_read_per_freq = (uint32_t)atof(optarg) * 2;
			bytes_to_read = bytes_to_read_per_freq * 3; // three recordings with different frequencies
			break;
		case 'S':
			sync_mode = 1;
			break;
		default:
			usage();
			break;
		}
	}

	if (argc <= optind) {
		usage();
	} else {
		filename = argv[optind];
	}
	
	/* Require both gain1 and gain2 parameters */
	if (gain1 == 0 || gain2 == 0) {
		fprintf(stderr, "ERROR: Both -1 (gain1) and -2 (gain2) parameters are required.\n");
		usage();
	}

	if(out_block_size < MINIMAL_BUF_LENGTH ||
	   out_block_size > MAXIMAL_BUF_LENGTH ){
		fprintf(stderr,
			"Output block size wrong value, falling back to default\n");
		fprintf(stderr,
			"Minimal length: %u\n", MINIMAL_BUF_LENGTH);
		fprintf(stderr,
			"Maximal length: %u\n", MAXIMAL_BUF_LENGTH);
		out_block_size = DEFAULT_BUF_LENGTH;
	}

	buffer = malloc(out_block_size * sizeof(uint8_t));

	if (!dev_given) {
		dev_index = verbose_device_search("0");
	}

	if (dev_index < 0) {
		exit(1);
	}

	r = rtlsdr_open(&dev, (uint32_t)dev_index);
	if (r < 0) {
		fprintf(stderr, "Failed to open rtlsdr device #%d.\n", dev_index);
		exit(1);
	}
#ifndef _WIN32
	sigact.sa_handler = sighandler;
	sigemptyset(&sigact.sa_mask);
	sigact.sa_flags = 0;
	sigaction(SIGINT, &sigact, NULL);
	sigaction(SIGTERM, &sigact, NULL);
	sigaction(SIGQUIT, &sigact, NULL);
	sigaction(SIGPIPE, &sigact, NULL);
#else
	SetConsoleCtrlHandler( (PHANDLER_ROUTINE) sighandler, TRUE );
#endif
	/* Set the sample rate */
	verbose_set_sample_rate(dev, samp_rate);

	/* Set the frequency */
	verbose_set_frequency(dev, frequency1);

	/* CRITICAL: Disable hardware AGC once at startup for TDOA consistency */
	r = rtlsdr_set_agc_mode(dev, 0);
	if (r != 0) {
		fprintf(stderr, "WARNING: Failed to disable hardware AGC.\n");
	} else {
		fprintf(stderr, "Hardware AGC disabled for TDOA operation.\n");
	}

	/* Set required gain variables for callback function */
	gain1_required = gain1;
	gain2_required = gain2;
	
	/* Set initial gain for frequency 1 */
	verbose_gain_set(dev, nearest_gain(dev, gain1));

	verbose_ppm_set(dev, ppm_error);

	if(strcmp(filename, "-") == 0) { /* Write samples to stdout */
		file = stdout;
#ifdef _WIN32
		_setmode(_fileno(stdin), _O_BINARY);
#endif
	} else {
		file = fopen(filename, "wb");
		if (!file) {
			fprintf(stderr, "Failed to open %s\n", filename);
			goto out;
		}
	}

	/* Reset endpoint before we start reading from it (mandatory) */
	verbose_reset_buffer(dev);

	if (sync_mode) {
		fprintf(stderr, "Reading samples in sync mode...\n");
		while (!do_exit) {
			r = rtlsdr_read_sync(dev, buffer, out_block_size, &n_read);
			if (r < 0) {
				fprintf(stderr, "WARNING: sync read failed.\n");
				break;
			}

			if ((bytes_to_read > 0) && (bytes_to_read < (uint32_t)n_read)) {
				n_read = bytes_to_read;
				do_exit = 1;
			}

			if (fwrite(buffer, 1, n_read, file) != (size_t)n_read) {
				fprintf(stderr, "Short write, samples lost, exiting!\n");
				break;
			}

			if ((uint32_t)n_read < out_block_size) {
				fprintf(stderr, "Short read, samples lost, exiting!\n");
				break;
			}

			if (bytes_to_read > 0)
				bytes_to_read -= n_read;
		}
	} else {
		fprintf(stderr, "Reading samples in async mode...\n");
		r = rtlsdr_read_async(dev, rtlsdr_callback, (void *)file,
				      0, out_block_size);
	}

	if (do_exit)
		fprintf(stderr, "\nUser cancel, exiting...\n");
	else
		fprintf(stderr, "\nLibrary error %d, exiting...\n", r);

	if (file != stdout)
		fclose(file);

	rtlsdr_close(dev);
	free (buffer);
out:
	return r >= 0 ? r : -r;
}
