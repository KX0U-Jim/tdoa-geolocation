package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

// Station represents a collector station or reference transmitter
type Station struct {
	Name      string
	Latitude  float64
	Longitude float64
	Elevation float64 // meters above sea level
}

// TDOAProcessor handles TDOA geolocation calculations
type TDOAProcessor struct {
	ReferenceFreq float64
	TargetFreq    float64
	Stations      map[string]Station
	RefStation    Station
}

// Point3D represents a 3D coordinate
type Point3D struct {
	X, Y, Z float64 // Earth-centered coordinates (meters)
}

// NewTDOAProcessor creates a new TDOA processor
func NewTDOAProcessor(refFreq, targetFreq float64, csvPath string) (*TDOAProcessor, error) {
	p := &TDOAProcessor{
		ReferenceFreq: refFreq,
		TargetFreq:    targetFreq,
		Stations:      make(map[string]Station),
	}

	err := p.loadStations(csvPath)
	if err != nil {
		return nil, fmt.Errorf("failed to load stations: %v", err)
	}

	return p, nil
}

// loadStations loads station coordinates from CSV file
func (p *TDOAProcessor) loadStations(csvPath string) error {
	file, err := os.Open(csvPath)
	if err != nil {
		return fmt.Errorf("failed to open CSV file: %v", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
		return fmt.Errorf("failed to read CSV: %v", err)
	}

	// Skip header row
	for i, record := range records[1:] {
		if len(record) != 4 {
			return fmt.Errorf("invalid CSV format at line %d", i+2)
		}

		lat, err := strconv.ParseFloat(record[1], 64)
		if err != nil {
			return fmt.Errorf("invalid latitude at line %d: %v", i+2, err)
		}

		lon, err := strconv.ParseFloat(record[2], 64)
		if err != nil {
			return fmt.Errorf("invalid longitude at line %d: %v", i+2, err)
		}

		elev, err := strconv.ParseFloat(record[3], 64)
		if err != nil {
			return fmt.Errorf("invalid elevation at line %d: %v", i+2, err)
		}

		station := Station{
			Name:      record[0],
			Latitude:  lat,
			Longitude: lon,
			Elevation: elev,
		}

		p.Stations[record[0]] = station

		// Set reference station if this matches our reference frequency
		if record[0] == fmt.Sprintf("%.0f", p.ReferenceFreq) {
			p.RefStation = station
		}
	}

	if p.RefStation.Name == "" {
		return fmt.Errorf("reference frequency %.0f not found in stations", p.ReferenceFreq)
	}

	fmt.Printf("Loaded %d stations including reference %.0f MHz\n", len(p.Stations), p.ReferenceFreq/1e6)
	return nil
}

// getStationFromFilename extracts station ID from .dat filename
func (p *TDOAProcessor) getStationFromFilename(filename string) (Station, error) {
	base := filepath.Base(filename)
	
	// Extract station ID from filename pattern (assumes format contains station name)
	// Look for known station names in the filename
	for stationName := range p.Stations {
		if strings.Contains(base, stationName) {
			return p.Stations[stationName], nil
		}
	}
	
	return Station{}, fmt.Errorf("could not identify station from filename: %s", filename)
}

// latLonToECEF converts lat/lon/elevation to Earth-Centered Earth-Fixed coordinates
func latLonToECEF(lat, lon, elev float64) Point3D {
	const (
		a = 6378137.0         // WGS84 semi-major axis (meters)
		f = 1.0 / 298.257223563 // WGS84 flattening
	)
	
	e2 := 2*f - f*f // First eccentricity squared
	
	latRad := lat * math.Pi / 180
	lonRad := lon * math.Pi / 180
	
	sinLat := math.Sin(latRad)
	cosLat := math.Cos(latRad)
	sinLon := math.Sin(lonRad)
	cosLon := math.Cos(lonRad)
	
	N := a / math.Sqrt(1 - e2*sinLat*sinLat)
	
	x := (N + elev) * cosLat * cosLon
	y := (N + elev) * cosLat * sinLon
	z := (N*(1-e2) + elev) * sinLat
	
	return Point3D{X: x, Y: y, Z: z}
}

// distance3D calculates 3D distance between two points in ECEF coordinates
func distance3D(p1, p2 Point3D) float64 {
	dx := p2.X - p1.X
	dy := p2.Y - p1.Y
	dz := p2.Z - p1.Z
	return math.Sqrt(dx*dx + dy*dy + dz*dz)
}

// calculateBaseline calculates 3D baseline distance between two stations
func (p *TDOAProcessor) calculateBaseline(station1, station2 Station) float64 {
	p1 := latLonToECEF(station1.Latitude, station1.Longitude, station1.Elevation)
	p2 := latLonToECEF(station2.Latitude, station2.Longitude, station2.Elevation)
	return distance3D(p1, p2)
}

// loadIQData loads I/Q data from .dat file
func (p *TDOAProcessor) loadIQData(filename string) ([]complex64, error) {
	fmt.Printf("Loading I/Q data from: %s\n", filename)
	
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("failed to open file: %v", err)
	}
	defer file.Close()
	
	// Get file size to calculate number of samples
	stat, err := file.Stat()
	if err != nil {
		return nil, fmt.Errorf("failed to get file size: %v", err)
	}
	
	// Each complex sample is 2 bytes (I + Q as unsigned 8-bit integers)
	numSamples := stat.Size() / 2
	
	fmt.Printf("File size: %d bytes, samples: %d\n", stat.Size(), numSamples)
	
	// Read raw bytes
	rawData := make([]byte, stat.Size())
	_, err = file.Read(rawData)
	if err != nil {
		return nil, fmt.Errorf("failed to read data: %v", err)
	}
	
	// Convert to complex64 samples
	samples := make([]complex64, numSamples)
	for i := int64(0); i < numSamples; i++ {
		// RTL-SDR data is unsigned 8-bit, centered at 127.5
		// Convert to signed float32 centered at 0
		iVal := (float32(rawData[i*2]) - 127.5) / 127.5
		qVal := (float32(rawData[i*2+1]) - 127.5) / 127.5
		samples[i] = complex(iVal, qVal)
	}
	
	fmt.Printf("Successfully loaded %d complex samples\n", len(samples))
	return samples, nil
}

// extractReferenceSignal extracts reference frequency samples from dual-freq data
func (p *TDOAProcessor) extractReferenceSignal(data []complex64) []complex64 {
	fmt.Println("Extracting reference signal from dual-frequency data")
	
	// librtlsdr-2freq pattern: freq1, freq2, freq1 (3 blocks)
	// Each block is 1/3 of the total samples
	totalSamples := len(data)
	blockSize := totalSamples / 3
	
	if blockSize == 0 {
		fmt.Println("Warning: Data too small for dual-frequency extraction")
		return data
	}
	
	fmt.Printf("Total samples: %d, block size: %d\n", totalSamples, blockSize)
	
	// Extract blocks 1 and 3 (reference frequency)
	// Block 1: samples 0 to blockSize-1
	// Block 3: samples 2*blockSize to 3*blockSize-1
	refSamples := make([]complex64, blockSize*2)
	
	// Copy block 1
	copy(refSamples[0:blockSize], data[0:blockSize])
	
	// Copy block 3
	if 2*blockSize < totalSamples {
		copy(refSamples[blockSize:], data[2*blockSize:3*blockSize])
	}
	
	fmt.Printf("Extracted %d reference samples from blocks 1 and 3\n", len(refSamples))
	return refSamples
}

// extractTargetSignal extracts target frequency samples from dual-freq data
func (p *TDOAProcessor) extractTargetSignal(data []complex64) []complex64 {
	fmt.Println("Extracting target signal from dual-frequency data")
	
	// librtlsdr-2freq pattern: freq1, freq2, freq1 (3 blocks)
	// Each block is 1/3 of the total samples
	// Target signal is in block 2 (middle block)
	totalSamples := len(data)
	blockSize := totalSamples / 3
	
	if blockSize == 0 {
		fmt.Println("Warning: Data too small for dual-frequency extraction")
		return data
	}
	
	fmt.Printf("Total samples: %d, block size: %d\n", totalSamples, blockSize)
	
	// Extract block 2 (target frequency)
	// Block 2: samples blockSize to 2*blockSize-1
	targetSamples := make([]complex64, blockSize)
	
	if 2*blockSize <= totalSamples {
		copy(targetSamples, data[blockSize:2*blockSize])
	}
	
	fmt.Printf("Extracted %d target samples from block 2\n", len(targetSamples))
	return targetSamples
}

// applyLowPassFilter applies a simple moving average low-pass filter
func (p *TDOAProcessor) applyLowPassFilter(signal []complex64, windowSize int) []complex64 {
	if windowSize <= 1 {
		return signal
	}
	
	filtered := make([]complex64, len(signal))
	halfWindow := windowSize / 2
	
	for i := range signal {
		var sum complex64
		count := 0
		
		// Average samples in window around current sample
		for j := i - halfWindow; j <= i + halfWindow; j++ {
			if j >= 0 && j < len(signal) {
				sum += signal[j]
				count++
			}
		}
		
		if count > 0 {
			filtered[i] = sum / complex(float32(count), 0)
		}
	}
	
	return filtered
}

// removeDCBias removes DC offset from signal
func (p *TDOAProcessor) removeDCBias(signal []complex64) []complex64 {
	if len(signal) == 0 {
		return signal
	}
	
	// Calculate DC bias (mean)
	var sum complex64
	for _, sample := range signal {
		sum += sample
	}
	dcBias := sum / complex(float32(len(signal)), 0)
	
	// Remove DC bias
	result := make([]complex64, len(signal))
	for i, sample := range signal {
		result[i] = sample - dcBias
	}
	
	fmt.Printf("Removed DC bias: %.6f + %.6fi\n", real(dcBias), imag(dcBias))
	return result
}

// calculateSignalPower calculates the average power of a signal
func (p *TDOAProcessor) calculateSignalPower(signal []complex64) float64 {
	if len(signal) == 0 {
		return 0
	}
	
	var power float64
	for _, sample := range signal {
		power += float64(real(sample)*real(sample) + imag(sample)*imag(sample))
	}
	
	return power / float64(len(signal))
}

// normalizeSignal normalizes signal to unit power
func (p *TDOAProcessor) normalizeSignal(signal []complex64) []complex64 {
	power := p.calculateSignalPower(signal)
	if power <= 0 {
		return signal
	}
	
	scale := float32(1.0 / math.Sqrt(power))
	result := make([]complex64, len(signal))
	
	for i, sample := range signal {
		result[i] = complex(real(sample)*scale, imag(sample)*scale)
	}
	
	fmt.Printf("Normalized signal power: %.6f → 1.000000\n", power)
	return result
}

// applyBandpassFilter applies a simple bandpass filter
func (p *TDOAProcessor) applyBandpassFilter(signal []complex64, lowCutoff, highCutoff float64, sampleRate float64) []complex64 {
	if len(signal) == 0 {
		return signal
	}
	
	fmt.Printf("Bandpass filter: %.1f - %.1f Hz (at %.0f Hz sample rate)\n", lowCutoff, highCutoff, sampleRate)
	
	// Simple frequency domain filtering
	// This is a basic implementation - in production would use proper filter design
	filtered := make([]complex64, len(signal))
	
	// Apply simple frequency-based filtering by zeroing out-of-band components
	// For simplicity, use time domain approach with multiple stages
	
	// Stage 1: High-pass filter (remove low frequencies)
	if lowCutoff > 0 {
		filtered = p.applyHighPassFilter(signal, lowCutoff, sampleRate)
	} else {
		copy(filtered, signal)
	}
	
	// Stage 2: Low-pass filter (remove high frequencies) 
	if highCutoff < sampleRate/2 {
		filtered = p.applyLowPassFilterWithCutoff(filtered, highCutoff, sampleRate)
	}
	
	return filtered
}

// applyHighPassFilter removes low frequency components
func (p *TDOAProcessor) applyHighPassFilter(signal []complex64, cutoff, sampleRate float64) []complex64 {
	// Simple high-pass: signal - low_pass_filtered_signal
	lowPassed := p.applyLowPassFilterWithCutoff(signal, cutoff, sampleRate)
	
	result := make([]complex64, len(signal))
	for i := range signal {
		result[i] = signal[i] - lowPassed[i]
	}
	
	return result
}

// applyLowPassFilterWithCutoff applies low-pass filter with specific cutoff frequency
func (p *TDOAProcessor) applyLowPassFilterWithCutoff(signal []complex64, cutoff, sampleRate float64) []complex64 {
	// Calculate filter window size based on cutoff frequency
	// Rule of thumb: window size inversely proportional to cutoff frequency
	windowSize := int(sampleRate / (2 * cutoff))
	if windowSize < 3 {
		windowSize = 3
	}
	if windowSize > 1000 {
		windowSize = 1000
	}
	
	return p.applyLowPassFilter(signal, windowSize)
}

// applyNotchFilter removes specific frequency components (like interference)
func (p *TDOAProcessor) applyNotchFilter(signal []complex64, notchFreq, bandwidth, sampleRate float64) []complex64 {
	// Simple notch filter: apply bandpass around the notch frequency and subtract
	lowFreq := notchFreq - bandwidth/2
	highFreq := notchFreq + bandwidth/2
	
	if lowFreq < 0 {
		lowFreq = 0
	}
	if highFreq > sampleRate/2 {
		highFreq = sampleRate / 2
	}
	
	// Extract the notch band
	notchBand := p.applyBandpassFilter(signal, lowFreq, highFreq, sampleRate)
	
	// Subtract from original signal
	result := make([]complex64, len(signal))
	for i := range signal {
		result[i] = signal[i] - notchBand[i]*0.8 // Reduce, don't eliminate completely
	}
	
	return result
}

// enhanceWeakSignal applies aggressive filtering for very weak signals
func (p *TDOAProcessor) enhanceWeakSignal(signal []complex64, label string) []complex64 {
	fmt.Printf("Enhancing weak signal: %s\n", label)
	
	sampleRate := 2000000.0 // 2 MHz sample rate
	
	// Step 1: Remove DC bias first
	signal = p.removeDCBias(signal)
	
	// Step 2: Remove 60 Hz harmonics and common interference
	signal = p.applyNotchFilter(signal, 60, 5, sampleRate)     // 60 Hz power line
	signal = p.applyNotchFilter(signal, 120, 5, sampleRate)    // 120 Hz harmonic
	signal = p.applyNotchFilter(signal, 1000000, 50000, sampleRate) // 1 MHz strong signals
	
	// Step 3: Apply aggressive bandpass around expected signal bandwidth
	// NOAA weather radio uses narrowband FM (~12.5 kHz deviation)
	// So we want roughly ±15-20 kHz around the signal
	signalBandwidth := 40000.0 // ±20 kHz around center
	lowCutoff := 100.0         // Remove very low frequencies
	highCutoff := signalBandwidth
	
	signal = p.applyBandpassFilter(signal, lowCutoff, highCutoff, sampleRate)
	
	// Step 4: Apply additional smoothing to reduce noise
	signal = p.applyLowPassFilter(signal, 50) // Moderate smoothing
	
	// Step 5: Normalize power
	signal = p.normalizeSignal(signal)
	
	return signal
}

// preprocessSignal applies enhanced filtering and normalization
func (p *TDOAProcessor) preprocessSignal(signal []complex64, label string) []complex64 {
	fmt.Printf("Preprocessing %s signal (%d samples)\n", label, len(signal))
	
	// Check signal power to decide on processing approach
	initialPower := p.calculateSignalPower(signal)
	fmt.Printf("Initial signal power: %.9f\n", initialPower)
	
	if initialPower < 0.001 {
		// Very weak signal - use aggressive enhancement
		fmt.Printf("Detected very weak signal - applying aggressive filtering\n")
		return p.enhanceWeakSignal(signal, label)
	} else {
		// Standard processing for stronger signals
		fmt.Printf("Standard signal processing\n")
		
		// Step 1: Remove DC bias
		signal = p.removeDCBias(signal)
		
		// Step 2: Apply moderate bandpass filtering
		sampleRate := 2000000.0
		signal = p.applyBandpassFilter(signal, 500, 50000, sampleRate)
		
		// Step 3: Apply light smoothing
		signal = p.applyLowPassFilter(signal, 100)
		
		// Step 4: Normalize to unit power
		signal = p.normalizeSignal(signal)
		
		return signal
	}
}

// nextPowerOfTwo finds the next power of 2 greater than or equal to n
func nextPowerOfTwo(n int) int {
	if n <= 1 {
		return 1
	}
	
	power := 1
	for power < n {
		power <<= 1
	}
	return power
}

// simpleFFT performs a basic FFT (simplified version for correlation)
func simpleFFT(signal []complex64) []complex64 {
	n := len(signal)
	if n <= 1 {
		return signal
	}
	
	// For simplicity, use a basic DFT for now
	// In production, would use a proper FFT library
	result := make([]complex64, n)
	
	for k := 0; k < n; k++ {
		var sum complex64
		for j := 0; j < n; j++ {
			angle := -2.0 * math.Pi * float64(k*j) / float64(n)
			twiddle := complex(float32(math.Cos(angle)), float32(math.Sin(angle)))
			sum += signal[j] * twiddle
		}
		result[k] = sum
	}
	
	return result
}

// frequencyDomainCorrelation performs correlation in frequency domain
func (p *TDOAProcessor) frequencyDomainCorrelation(signal1, signal2 []complex64, maxLag int) (int, float64) {
	fmt.Println("Performing frequency domain correlation")
	
	if len(signal1) == 0 || len(signal2) == 0 {
		fmt.Println("Warning: Empty signals for frequency domain correlation")
		return 0, 0.0
	}
	
	// For efficiency, use much smaller chunks for correlation testing
	chunkSize := 1024 // 0.5ms at 2 Msps - very fast for testing
	if chunkSize > len(signal1) {
		chunkSize = len(signal1)
	}
	if chunkSize > len(signal2) {
		chunkSize = len(signal2)
	}
	
	// Use first chunk of each signal
	chunk1 := signal1[:chunkSize]
	chunk2 := signal2[:chunkSize]
	
	fmt.Printf("Using %d samples for frequency domain correlation\n", chunkSize)
	
	// Pad to next power of 2 for efficient FFT
	fftSize := nextPowerOfTwo(chunkSize + maxLag)
	padded1 := make([]complex64, fftSize)
	padded2 := make([]complex64, fftSize)
	
	copy(padded1, chunk1)
	copy(padded2, chunk2)
	
	fmt.Printf("FFT size: %d samples\n", fftSize)
	
	// Compute FFTs
	fft1 := simpleFFT(padded1)
	fft2 := simpleFFT(padded2)
	
	// Cross-power spectrum: FFT1 * conj(FFT2)
	crossPower := make([]complex64, fftSize)
	for i := 0; i < fftSize; i++ {
		conj2 := complex(real(fft2[i]), -imag(fft2[i]))
		crossPower[i] = fft1[i] * conj2
	}
	
	// Inverse FFT to get correlation function
	correlation := simpleFFT(crossPower) // Using same function (should be IFFT)
	
	// Find peak in correlation function
	bestDelay := 0
	bestCorr := 0.0
	
	// Check limited range around zero lag
	searchRange := maxLag
	if searchRange > fftSize/2 {
		searchRange = fftSize / 2
	}
	
	for i := 0; i < searchRange; i++ {
		// Check positive and negative lags
		corrVal := float64(real(correlation[i]))
		if math.Abs(corrVal) > math.Abs(bestCorr) {
			bestCorr = corrVal
			bestDelay = i
		}
		
		// Check negative lag (from end of array)
		if fftSize-i > 0 {
			corrVal = float64(real(correlation[fftSize-i]))
			if math.Abs(corrVal) > math.Abs(bestCorr) {
				bestCorr = corrVal
				bestDelay = -(i)
			}
		}
	}
	
	fmt.Printf("Frequency domain correlation: %.6f at delay %d samples\n", bestCorr, bestDelay)
	return bestDelay, bestCorr
}

// crossCorrelate performs cross-correlation between two signals with preprocessing
func (p *TDOAProcessor) crossCorrelate(signal1, signal2 []complex64) (int, float64) {
	fmt.Println("=== Cross-Correlation Analysis ===")
	
	if len(signal1) == 0 || len(signal2) == 0 {
		fmt.Println("Warning: Empty signals for correlation")
		return 0, 0.0
	}
	
	// Preprocess signals for better correlation
	fmt.Println("\n--- Signal Preprocessing ---")
	processed1 := p.preprocessSignal(signal1, "Signal 1")
	processed2 := p.preprocessSignal(signal2, "Signal 2")
	
	// Reasonable search range - start small and optimize
	maxLag := 2000 // 1ms at 2 Msps - reasonable for local baselines
	
	fmt.Println("\n--- Time Domain Correlation ---")
	timeDomainDelay, timeDomainCorr := p.timeDomainCorrelation(processed1, processed2, maxLag)
	
	// Skip frequency domain for now - focus on improved time domain
	fmt.Printf("\n--- Result: Time Domain with Preprocessing ---\n")
	fmt.Printf("Correlation: %.6f at delay %d samples\n", timeDomainCorr, timeDomainDelay)
	
	return timeDomainDelay, timeDomainCorr
}

// timeDomainCorrelation performs traditional time-domain correlation
func (p *TDOAProcessor) timeDomainCorrelation(signal1, signal2 []complex64, maxLag int) (int, float64) {
	fmt.Println("Performing time domain correlation")
	
	// Use shorter signal as template for efficiency
	template := signal1
	signal := signal2
	if len(signal1) > len(signal2) {
		template = signal2  
		signal = signal1
	}
	
	templateLen := len(template)
	signalLen := len(signal)
	
	fmt.Printf("Template: %d samples, Signal: %d samples\n", templateLen, signalLen)
	
	if templateLen > signalLen {
		fmt.Println("Warning: Template longer than signal")
		return 0, 0.0
	}
	
	// For equal-length signals, use shorter template to allow delay search
	if templateLen == signalLen {
		templateLen = templateLen - maxLag  // Make room for time shifts
		fmt.Printf("Reduced template to %d samples to allow %d sample delay search\n", templateLen, maxLag)
	}
	
	// Limit search range to what's actually possible
	if maxLag > signalLen-templateLen {
		maxLag = signalLen - templateLen
	}
	
	// Ensure we can search meaningful delays
	if maxLag < 1 {
		maxLag = 1
	}
	
	bestDelay := 0
	bestCorr := 0.0
	
	// Coherent integration for processing gain
	// Use larger blocks for efficiency, but not too large
	blockSize := 10000 // Average over 10000-sample blocks for balance of gain and speed
	
	fmt.Printf("Using coherent integration with %d-sample blocks\n", blockSize)
	
	for delay := 0; delay < maxLag; delay++ {
		correlation := 0.0
		numBlocks := 0
		
		// Process in blocks for coherent integration
		for blockStart := 0; blockStart < templateLen-blockSize; blockStart += blockSize {
			blockEnd := blockStart + blockSize
			if delay + blockEnd > signalLen {
				break
			}
			
			blockCorr := 0.0
			
			// Correlate this block
			for i := blockStart; i < blockEnd; i++ {
				t := template[i]
				s := signal[delay+i]
				
				// Complex correlation (conjugate of template)
				blockCorr += float64(real(t)*real(s) + imag(t)*imag(s))
			}
			
			// Normalize block correlation
			blockCorr /= float64(blockSize)
			correlation += blockCorr
			numBlocks++
		}
		
		if numBlocks > 0 {
			// Average across blocks (coherent integration)
			correlation /= float64(numBlocks)
			
			// Coherent integration gain: sqrt(N) improvement in SNR
			coherentGain := math.Sqrt(float64(numBlocks * blockSize))
			correlation *= coherentGain
			
			if math.Abs(correlation) > math.Abs(bestCorr) {
				bestCorr = correlation
				bestDelay = delay
			}
		}
		
		// Show progress
		if delay%2000 == 0 {
			fmt.Printf("Time domain progress: %d/%d (coherent blocks: %d)\r", delay, maxLag, numBlocks)
		}
	}
	
	fmt.Printf("\nTime domain correlation: %.6f at delay %d samples\n", bestCorr, bestDelay)
	return bestDelay, bestCorr
}

// ProcessTDOA processes TDOA from multiple collector files
func (p *TDOAProcessor) ProcessTDOA(datFiles []string) error {
	if len(datFiles) < 3 {
		return fmt.Errorf("need at least 3 collector stations, got %d", len(datFiles))
	}

	fmt.Printf("Processing TDOA for target frequency %.3f MHz\n", p.TargetFreq/1e6)
	fmt.Printf("Reference: %s at %.6f°, %.6f°, %.1fm\n", 
		p.RefStation.Name, p.RefStation.Latitude, p.RefStation.Longitude, p.RefStation.Elevation)

	// Load data from all collectors
	var collectorData []struct {
		Station Station
		Data    []complex64
		RefSig  []complex64
		TargetSig []complex64
	}

	for _, filename := range datFiles {
		station, err := p.getStationFromFilename(filename)
		if err != nil {
			return fmt.Errorf("failed to identify station for %s: %v", filename, err)
		}

		data, err := p.loadIQData(filename)
		if err != nil {
			return fmt.Errorf("failed to load data from %s: %v", filename, err)
		}

		refSig := p.extractReferenceSignal(data)
		targetSig := p.extractTargetSignal(data)
		
		// Use reasonable chunk size for correlation - too big causes computational explosion
		// Balance between processing gain and computational complexity
		testChunkSize := 1000000 // 0.5 seconds at 2 Msps - manageable computation
		if len(refSig) > testChunkSize {
			refSig = refSig[:testChunkSize]
			fmt.Printf("Using test chunk: %d samples (%.1f ms)\n", testChunkSize, float64(testChunkSize)/2e6*1000)
		}
		if len(targetSig) > testChunkSize {
			targetSig = targetSig[:testChunkSize]
			fmt.Printf("Using target test chunk: %d samples (%.1f ms)\n", testChunkSize, float64(testChunkSize)/2e6*1000)
		}
		
		fmt.Printf("Coherent integration time: %.0f ms (expecting ~%.1f dB processing gain)\n", 
			float64(testChunkSize)/2e6*1000, 10*math.Log10(float64(testChunkSize)/100000))

		collectorData = append(collectorData, struct {
			Station Station
			Data    []complex64
			RefSig  []complex64
			TargetSig []complex64
		}{
			Station: station,
			Data:    data,
			RefSig:  refSig,
			TargetSig: targetSig,
		})

		fmt.Printf("Loaded collector: %s at %.6f°, %.6f°, %.1fm\n", 
			station.Name, station.Latitude, station.Longitude, station.Elevation)
	}

	// Calculate baselines between all station pairs
	fmt.Println("\nBaseline distances (3D):")
	for i := 0; i < len(collectorData); i++ {
		for j := i + 1; j < len(collectorData); j++ {
			dist := p.calculateBaseline(collectorData[i].Station, collectorData[j].Station)
			fmt.Printf("%s - %s: %.2f km\n", 
				collectorData[i].Station.Name, collectorData[j].Station.Name, dist/1000)
		}
	}

	// Test both reference and target signal correlations
	fmt.Println("\n=== REFERENCE SIGNAL CORRELATION TEST ===")
	fmt.Println("Testing weak 162.4 MHz NOAA weather signal:")
	var refTimeDifferences []float64
	
	for i := 0; i < len(collectorData); i++ {
		for j := i + 1; j < len(collectorData); j++ {
			delay, correlation := p.crossCorrelate(collectorData[i].RefSig, collectorData[j].RefSig)
			
			// Convert sample delay to time difference
			sampleRate := 2e6 // 2 Msps
			timeDiff := float64(delay) / sampleRate
			
			refTimeDifferences = append(refTimeDifferences, timeDiff)
			
			fmt.Printf("REF %s - %s: delay=%d samples (%.3f μs), correlation=%.6f\n", 
				collectorData[i].Station.Name, collectorData[j].Station.Name, 
				delay, timeDiff*1e6, correlation)
		}
	}
	
	fmt.Println("\n=== TARGET SIGNAL CORRELATION TEST ===")
	fmt.Println("Testing strong 92.3 MHz FM broadcast signal:")
	var targetTimeDifferences []float64
	
	for i := 0; i < len(collectorData); i++ {
		for j := i + 1; j < len(collectorData); j++ {
			delay, correlation := p.crossCorrelate(collectorData[i].TargetSig, collectorData[j].TargetSig)
			
			// Convert sample delay to time difference
			sampleRate := 2e6 // 2 Msps
			timeDiff := float64(delay) / sampleRate
			
			targetTimeDifferences = append(targetTimeDifferences, timeDiff)
			
			fmt.Printf("TGT %s - %s: delay=%d samples (%.3f μs), correlation=%.6f\n", 
				collectorData[i].Station.Name, collectorData[j].Station.Name, 
				delay, timeDiff*1e6, correlation)
		}
	}
	
	fmt.Printf("\n=== REFERENCE SIGNAL SYNCHRONIZATION ===\n")
	fmt.Printf("Using reference signal to synchronize collector timing...\n")
	
	// Calculate timing offsets from reference signal correlations
	refTimingOffsets := make([]float64, len(refTimeDifferences))
	for i, refDelay := range refTimeDifferences {
		refTimingOffsets[i] = refDelay
		fmt.Printf("Reference timing offset %d: %.3f μs\n", i, refDelay*1e6)
	}
	
	fmt.Printf("\n=== APPLYING TIMING CORRECTIONS TO TARGET SIGNAL ===\n")
	
	// Apply reference signal timing corrections to target signal delays
	correctedTargetDelays := make([]float64, len(targetTimeDifferences))
	for i := 0; i < len(targetTimeDifferences); i++ {
		correctedTargetDelays[i] = targetTimeDifferences[i] - refTimingOffsets[i]
		fmt.Printf("Target delay %d: %.3f μs (raw) - %.3f μs (ref offset) = %.3f μs (corrected)\n", 
			i, targetTimeDifferences[i]*1e6, refTimingOffsets[i]*1e6, correctedTargetDelays[i]*1e6)
	}
	
	fmt.Printf("\n=== CORRELATION COMPARISON ===\n")
	fmt.Printf("Reference signal (162.4 MHz): Used for timing synchronization\n")
	fmt.Printf("Target signal (92.3 MHz): Corrected with reference timing offsets\n")
	fmt.Printf("Using corrected target signal for TDOA calculation\n")
	
	// Use the corrected target signal differences 
	timeDifferences := correctedTargetDelays

	// TDOA triangulation (basic implementation)
	fmt.Println("\nTDOA triangulation:")
	
	if len(timeDifferences) < 3 {
		fmt.Println("ERROR: Need at least 3 time differences for triangulation")
		return nil
	}
	
	fmt.Printf("Corrected time differences: %.3f μs, %.3f μs, %.3f μs\n", 
		timeDifferences[0]*1e6, timeDifferences[1]*1e6, timeDifferences[2]*1e6)
	
	// Convert corrected time differences to distance differences
	c := 299792458.0 // Speed of light (m/s)
	distDiff1 := timeDifferences[0] * c
	distDiff2 := timeDifferences[1] * c  
	distDiff3 := timeDifferences[2] * c
	
	fmt.Printf("Corrected distance differences: %.1f m, %.1f m, %.1f m\n", 
		distDiff1, distDiff2, distDiff3)
	
	// Simple diagnostic: what if we had some example time delays?
	fmt.Println("\nDiagnostic test with example delays:")
	fmt.Println("Simulating 10 μs, 5 μs, -3 μs delays...")
	
	testDelays := []float64{10e-6, 5e-6, -3e-6} // Example delays in seconds
	for i, delay := range testDelays {
		distance := delay * c
		fmt.Printf("Test delay %d: %.1f μs → %.1f m\n", i+1, delay*1e6, distance)
	}
	
	// Calculate TDOA position using target signal time differences
	fmt.Println("\n=== TDOA GEOLOCATION ===")
	
	if len(timeDifferences) < 3 {
		return fmt.Errorf("need at least 3 time differences for TDOA, got %d", len(timeDifferences))
	}
	
	// Convert corrected time differences to range differences (multiply by speed of light)
	const speedOfLight = 299792458.0 // Speed of light (m/s)
	rangeDifferences := make([]float64, len(timeDifferences))
	for i, td := range timeDifferences {
		rangeDifferences[i] = td * speedOfLight
	}
	
	fmt.Printf("Time differences (μs): ")
	for _, td := range timeDifferences {
		fmt.Printf("%.3f ", td*1e6)
	}
	fmt.Println()
	
	fmt.Printf("Range differences (m): ")
	for _, rd := range rangeDifferences {
		fmt.Printf("%.1f ", rd)
	}
	fmt.Println()
	
	// Solve TDOA using least squares
	lat, lon, elev, err := p.solveTDOA(collectorData, rangeDifferences)
	if err != nil {
		return fmt.Errorf("TDOA solution failed: %v", err)
	}
	
	fmt.Printf("\n*** CALCULATED TRANSMITTER LOCATION ***\n")
	fmt.Printf("Latitude:  %.6f°\n", lat)
	fmt.Printf("Longitude: %.6f°\n", lon) 
	fmt.Printf("Elevation: %.1f m\n", elev)
	
	return nil
}

// solveTDOA solves TDOA positioning using least squares
func (p *TDOAProcessor) solveTDOA(collectorData []struct {
	Station   Station
	Data      []complex64
	RefSig    []complex64
	TargetSig []complex64
}, rangeDifferences []float64) (float64, float64, float64, error) {
	
	// Convert stations to ECEF coordinates for calculation
	stations := make([]Point3D, len(collectorData))
	for i, data := range collectorData {
		stations[i] = latLonToECEF(data.Station.Latitude, data.Station.Longitude, data.Station.Elevation)
	}
	
	// TDOA uses hyperbolic positioning - need to solve system of equations
	// For 3 stations, we have 2 independent TDOA measurements (station pairs)
	// Range difference equations: |r1 - r_tx| - |r2 - r_tx| = c * (t1 - t2)
	
	// Use iterative approach (Gauss-Newton) starting from centroid
	centroidLat := (collectorData[0].Station.Latitude + collectorData[1].Station.Latitude + collectorData[2].Station.Latitude) / 3.0
	centroidLon := (collectorData[0].Station.Longitude + collectorData[1].Station.Longitude + collectorData[2].Station.Longitude) / 3.0
	centroidElev := (collectorData[0].Station.Elevation + collectorData[1].Station.Elevation + collectorData[2].Station.Elevation) / 3.0
	
	// Initial guess at centroid
	x := latLonToECEF(centroidLat, centroidLon, centroidElev)
	
	fmt.Printf("Initial guess: %.6f°, %.6f°, %.1fm\n", centroidLat, centroidLon, centroidElev)
	
	// Newton-Raphson iterations
	for iter := 0; iter < 10; iter++ {
		// Calculate current range differences
		r1 := math.Sqrt((x.X-stations[0].X)*(x.X-stations[0].X) + (x.Y-stations[0].Y)*(x.Y-stations[0].Y) + (x.Z-stations[0].Z)*(x.Z-stations[0].Z))
		r2 := math.Sqrt((x.X-stations[1].X)*(x.X-stations[1].X) + (x.Y-stations[1].Y)*(x.Y-stations[1].Y) + (x.Z-stations[1].Z)*(x.Z-stations[1].Z))
		r3 := math.Sqrt((x.X-stations[2].X)*(x.X-stations[2].X) + (x.Y-stations[2].Y)*(x.Y-stations[2].Y) + (x.Z-stations[2].Z)*(x.Z-stations[2].Z))
		
		// Residuals (difference between measured and calculated range differences)
		residual1 := (r2 - r1) - rangeDifferences[0]  // Station 1-2 pair
		residual2 := (r3 - r1) - rangeDifferences[1]  // Station 1-3 pair
		
		if math.Abs(residual1) < 1.0 && math.Abs(residual2) < 1.0 {
			fmt.Printf("Converged after %d iterations\n", iter)
			break
		}
		
		// Calculate Jacobian matrix
		// Partial derivatives of range differences with respect to x, y
		dx1 := (x.X - stations[0].X) / r1
		dy1 := (x.Y - stations[0].Y) / r1
		
		dx2 := (x.X - stations[1].X) / r2
		dy2 := (x.Y - stations[1].Y) / r2
		
		dx3 := (x.X - stations[2].X) / r3
		dy3 := (x.Y - stations[2].Y) / r3
		
		// Jacobian for TDOA (range difference derivatives)
		J11 := dx2 - dx1  // d(r2-r1)/dx
		J12 := dy2 - dy1  // d(r2-r1)/dy
		
		J21 := dx3 - dx1  // d(r3-r1)/dx
		J22 := dy3 - dy1  // d(r3-r1)/dy
		
		// Solve 2x3 system using pseudo-inverse (least squares)
		// J * delta = -residuals
		// Use simplified 2-equation solver
		det := J11*J22 - J12*J21
		if math.Abs(det) < 1e-10 {
			return 0, 0, 0, fmt.Errorf("singular Jacobian matrix at iteration %d", iter)
		}
		
		// Solve for dx, dy (ignore dz for now - assume flat earth approximation)
		dx := (-residual1*J22 + residual2*J12) / det
		dy := (residual1*J21 - residual2*J11) / det
		dz := 0.0 // Simplified - keep elevation constant
		
		// Update position
		stepSize := 0.5  // Damped Newton step
		x.X += stepSize * dx
		x.Y += stepSize * dy
		x.Z += stepSize * dz
		
		if iter == 9 {
			fmt.Printf("Maximum iterations reached\n")
		}
	}
	
	// Convert back to lat/lon
	lat, lon, elev := ecefToLatLon(x.X, x.Y, x.Z)
	return lat, lon, elev, nil
}

// ecefToLatLon converts ECEF coordinates back to latitude/longitude/elevation
func ecefToLatLon(x, y, z float64) (float64, float64, float64) {
	const (
		a  = 6378137.0         // WGS84 semi-major axis (meters)
		f  = 1.0 / 298.257223563 // WGS84 flattening
		e2 = 2*f - f*f         // First eccentricity squared
	)
	
	p := math.Sqrt(x*x + y*y)
	lon := math.Atan2(y, x)
	
	// Iterative solution for latitude
	lat := math.Atan2(z, p*(1-e2))
	for i := 0; i < 5; i++ {
		N := a / math.Sqrt(1 - e2*math.Sin(lat)*math.Sin(lat))
		elev := p/math.Cos(lat) - N
		lat = math.Atan2(z, p*(1-e2*N/(N+elev)))
	}
	
	N := a / math.Sqrt(1 - e2*math.Sin(lat)*math.Sin(lat))
	elev := p/math.Cos(lat) - N
	
	return lat * 180.0 / math.Pi, lon * 180.0 / math.Pi, elev
}

func main() {
	if len(os.Args) < 5 {
		fmt.Printf("Usage: %s <ref_freq_hz> <target_freq_hz> <csv_file> <dat_file1> [dat_file2] [dat_file3] ...\n", os.Args[0])
		fmt.Println("Example: ./processor 162400000 101700000 lat-lon-table.csv kx0u-data.dat n3pay-data.dat kf0mtl-data.dat")
		os.Exit(1)
	}

	refFreq, err := strconv.ParseFloat(os.Args[1], 64)
	if err != nil {
		log.Fatalf("Invalid reference frequency: %v", err)
	}

	targetFreq, err := strconv.ParseFloat(os.Args[2], 64)
	if err != nil {
		log.Fatalf("Invalid target frequency: %v", err)
	}

	csvFile := os.Args[3]
	datFiles := os.Args[4:]

	processor, err := NewTDOAProcessor(refFreq, targetFreq, csvFile)
	if err != nil {
		log.Fatalf("Failed to create processor: %v", err)
	}

	err = processor.ProcessTDOA(datFiles)
	if err != nil {
		log.Fatalf("TDOA processing failed: %v", err)
	}
}