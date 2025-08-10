package main

import (
	"flag"
	"fmt"
	"os"
	"os/exec"
	"strconv"
	"strings"
	"time"
)

func main() {
	// Parse flags
	duration := flag.Int("duration", 30, "Collection duration in seconds (max 100)")
	gain := flag.Float64("gain", 0, "RTL-SDR gain in dB (0 = auto gain)")
	gain1 := flag.Float64("gain1", 0, "RTL-SDR gain for reference frequency in dB (0 = use --gain)")
	gain2 := flag.Float64("gain2", 0, "RTL-SDR gain for target frequency in dB (0 = use --gain)")
	flag.Parse()

	args := flag.Args()
	if len(args) < 4 {
		fmt.Printf("Usage: %s [--duration=seconds] [--gain=dB] [--gain1=dB] [--gain2=dB] <reference_freq_hz> <target_freq_hz> <start_epoch_seconds> <station_id>\n", os.Args[0])
		fmt.Printf("Example: %s --duration=60 --gain1=20.0 --gain2=45.0 96900000 162550000 1723234567 kx0u\n", os.Args[0])
		fmt.Printf("Gain: 0 = auto gain, or specify dB value (e.g., 49.6)\n")
		fmt.Printf("gain1/gain2: Separate gains for reference/target frequencies (overrides --gain)\n")
		fmt.Printf("Output: [station_id]-[start_epoch].dat\n")
		os.Exit(1)
	}

	if *duration > 100 {
		fmt.Printf("Error: Duration exceeds maximum of 100 seconds\n")
		os.Exit(1)
	}

	// Parse reference frequency
	refFreq, err := strconv.ParseUint(args[0], 10, 32)
	if err != nil {
		fmt.Printf("Error: Invalid reference frequency: %s\n", args[0])
		os.Exit(1)
	}

	// Parse target frequency  
	targetFreq, err := strconv.ParseUint(args[1], 10, 32)
	if err != nil {
		fmt.Printf("Error: Invalid target frequency: %s\n", args[1])
		os.Exit(1)
	}

	// Parse start time
	startTime, err := strconv.ParseInt(args[2], 10, 64)
	if err != nil {
		fmt.Printf("Error: Invalid start time: %s\n", args[2])
		os.Exit(1)
	}

	// Station ID
	stationID := args[3]
	
	// Output filename
	filename := fmt.Sprintf("%s-%d.dat", stationID, startTime)

	fmt.Printf("Collector starting:\n")
	fmt.Printf("Reference: %d Hz\n", refFreq)
	fmt.Printf("Target: %d Hz\n", targetFreq)
	fmt.Printf("Start: %s\n", time.Unix(startTime, 0))
	fmt.Printf("Duration: %d seconds\n", *duration)
	if *gain1 != 0 || *gain2 != 0 {
		fmt.Printf("Gain1 (ref): %.1f dB", *gain1)
		if *gain1 == 0 { fmt.Printf(" (auto)") }
		fmt.Printf(", Gain2 (target): %.1f dB", *gain2)  
		if *gain2 == 0 { fmt.Printf(" (auto)") }
		fmt.Printf("\n")
	} else {
		fmt.Printf("Gain: %.1f dB", *gain)
		if *gain == 0 { fmt.Printf(" (auto)") }
		fmt.Printf("\n")
	}
	fmt.Printf("Station: %s\n", stationID)
	fmt.Printf("Output: %s\n", filename)

	// Constants for sample collection
	const sampleRate = 2000000 // 2 Msps
	const bytesPerSample = 2   // I + Q bytes
	const switchInterval = 2   // seconds per frequency
	
	// Calculate total samples needed
	totalSamples := sampleRate * (*duration)
	totalBytes := totalSamples * bytesPerSample
	
	fmt.Printf("Allocating memory: %d samples (%d MB)\n", totalSamples, totalBytes/1024/1024)
	
	// Pre-allocate sample buffer
	sampleBuffer := make([]byte, totalBytes)
	if len(sampleBuffer) != totalBytes {
		fmt.Printf("Error: Failed to allocate %d bytes of memory\n", totalBytes)
		os.Exit(1)
	}
	
	fmt.Printf("Memory allocated successfully\n")
	
	// Calculate switching pattern
	samplesPerSwitch := sampleRate * switchInterval
	numSwitches := (*duration) / switchInterval
	if (*duration) % switchInterval != 0 {
		numSwitches++ // Handle remainder
	}
	
	fmt.Printf("Switch pattern: %d samples per frequency (%d seconds)\n", samplesPerSwitch, switchInterval)
	fmt.Printf("Total switches: %d\n", numSwitches)

	// Wait for start time
	fmt.Printf("Waiting for start time...\n")
	for time.Now().Unix() < startTime {
		time.Sleep(100 * time.Millisecond)
	}
	
	fmt.Printf("Starting collection at %s\n", time.Now().Format(time.RFC3339))
	
	// Calculate samples needed for rtl_sdr (it collects 3×n samples total)
	samplesPerFreq := totalSamples / 3
	
	// Build rtl_sdr command
	rtlSdrPath := "librtlsdr-2freq/build/src/rtl_sdr" // Path to our built binary
	cmd := []string{
		rtlSdrPath,
		"-f", fmt.Sprintf("%d", refFreq),
		"-h", fmt.Sprintf("%d", targetFreq), 
		"-s", fmt.Sprintf("%d", sampleRate),
	}
	
	// Add gain parameters based on what was specified
	if *gain1 != 0 {
		cmd = append(cmd, "-1", fmt.Sprintf("%.1f", *gain1))
	}
	if *gain2 != 0 {
		cmd = append(cmd, "-2", fmt.Sprintf("%.1f", *gain2))
	}
	if *gain1 == 0 && *gain2 == 0 && *gain != 0 {
		cmd = append(cmd, "-g", fmt.Sprintf("%.1f", *gain))
	}
	
	cmd = append(cmd, "-n", fmt.Sprintf("%d", samplesPerFreq))
	cmd = append(cmd, filename)
	
	fmt.Printf("Executing: %s\n", strings.Join(cmd, " "))
	
	// Execute rtl_sdr command
	rtlCmd := exec.Command(cmd[0], cmd[1:]...)
	err = rtlCmd.Run()
	if err != nil {
		fmt.Printf("Error executing rtl_sdr: %v\n", err)
		os.Exit(1)
	}
	
	// Print completion timestamp with millisecond precision
	completionTime := time.Now()
	epochMillis := completionTime.UnixMilli()
	fmt.Printf("Collection completed at: %d (epoch milliseconds)\n", epochMillis)
	fmt.Printf("Data saved to: %s\n", filename)
	
	_ = sampleBuffer // Remove this once we integrate better
}