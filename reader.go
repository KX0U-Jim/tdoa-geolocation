package main

import (
	"fmt"
	"os"
	"strconv"
)

func main() {
	if len(os.Args) < 2 {
		fmt.Printf("Usage: %s <data_file.dat> [expected_duration_seconds]\n", os.Args[0])
		fmt.Printf("Example: %s kx0u-1723234567.dat 30\n", os.Args[0])
		os.Exit(1)
	}

	filename := os.Args[1]
	
	// Try to parse duration from command line or filename
	expectedDuration := 30 // default
	if len(os.Args) > 2 {
		if dur, err := strconv.Atoi(os.Args[2]); err == nil {
			expectedDuration = dur
		}
	}
	
	fmt.Printf("Analyzing data file: %s\n", filename)
	fmt.Printf("Expected duration: %d seconds\n", expectedDuration)
	
	// Read and analyze the file
	err := analyzeDataFile(filename, expectedDuration)
	if err != nil {
		fmt.Printf("Error analyzing file: %v\n", err)
		os.Exit(1)
	}
}

func analyzeDataFile(filename string, expectedDuration int) error {
	// Open file
	file, err := os.Open(filename)
	if err != nil {
		return fmt.Errorf("failed to open file: %v", err)
	}
	defer file.Close()

	// Get file size
	fileInfo, err := file.Stat()
	if err != nil {
		return fmt.Errorf("failed to get file info: %v", err)
	}
	
	fileSize := fileInfo.Size()
	fmt.Printf("File size: %d bytes (%d MB)\n", fileSize, fileSize/1024/1024)
	
	// Constants from collector
	const sampleRate = 2000000 // 2 Msps
	const bytesPerSample = 2   // I + Q bytes
	
	// Calculate expected values
	expectedTotalSamples := sampleRate * expectedDuration
	expectedTotalBytes := int64(expectedTotalSamples * bytesPerSample)
	expectedSamplesPerFreq := expectedTotalSamples / 3 // librtlsdr-2freq pattern
	
	fmt.Printf("Expected total samples: %d\n", expectedTotalSamples)
	fmt.Printf("Expected total bytes: %d\n", expectedTotalBytes)
	fmt.Printf("Expected samples per frequency: %d\n", expectedSamplesPerFreq)
	
	// Validate file size
	actualSamples := fileSize / bytesPerSample
	fmt.Printf("Actual samples in file: %d\n", actualSamples)
	
	if fileSize == expectedTotalBytes {
		fmt.Printf("✓ File size matches expected exactly\n")
	} else {
		fmt.Printf("✗ File size mismatch. Expected: %d, Got: %d (diff: %d bytes)\n", 
			expectedTotalBytes, fileSize, fileSize-expectedTotalBytes)
	}
	
	// Check if file size is multiple of sample size
	if fileSize%bytesPerSample == 0 {
		fmt.Printf("✓ File size is valid multiple of sample size\n")
	} else {
		fmt.Printf("✗ File size is not a multiple of sample size - data may be corrupted\n")
	}
	
	// Validate librtlsdr-2freq pattern (3×n samples)
	if actualSamples%3 == 0 {
		samplesPerBlock := actualSamples / 3
		fmt.Printf("✓ Sample count follows 3×n pattern: %d samples per frequency block\n", samplesPerBlock)
	} else {
		fmt.Printf("✗ Sample count doesn't follow 3×n pattern - switching may have failed\n")
	}
	
	// Analyze sample data
	err = analyzeSampleData(file, int(actualSamples))
	if err != nil {
		return fmt.Errorf("failed to analyze sample data: %v", err)
	}
	
	return nil
}

func analyzeSampleData(file *os.File, totalSamples int) error {
	fmt.Printf("\nAnalyzing sample data quality...\n")
	
	// Read a subset of samples for analysis
	sampleSize := 10000 // Analyze first 10k samples
	if totalSamples < sampleSize {
		sampleSize = totalSamples
	}
	
	samples := make([]byte, sampleSize*2) // 2 bytes per sample
	_, err := file.ReadAt(samples, 0)
	if err != nil {
		return fmt.Errorf("failed to read sample data: %v", err)
	}
	
	// Analyze I/Q values
	var iSum, qSum int64
	var iMin, iMax, qMin, qMax byte = 255, 0, 255, 0
	
	for i := 0; i < sampleSize; i++ {
		iVal := samples[i*2]     // I component
		qVal := samples[i*2+1]   // Q component
		
		iSum += int64(iVal)
		qSum += int64(qVal)
		
		if iVal < iMin { iMin = iVal }
		if iVal > iMax { iMax = iVal }
		if qVal < qMin { qMin = qVal }
		if qVal > qMax { qMax = qVal }
	}
	
	// Calculate averages
	iAvg := float64(iSum) / float64(sampleSize)
	qAvg := float64(qSum) / float64(sampleSize)
	
	fmt.Printf("I component: min=%d, max=%d, avg=%.1f\n", iMin, iMax, iAvg)
	fmt.Printf("Q component: min=%d, max=%d, avg=%.1f\n", qMin, qMax, qAvg)
	
	// Check for reasonable signal characteristics
	iRange := iMax - iMin
	qRange := qMax - qMin
	
	if iRange > 10 && qRange > 10 {
		fmt.Printf("✓ I/Q components show good dynamic range\n")
	} else {
		fmt.Printf("✗ I/Q components have poor dynamic range - possible hardware issue\n")
	}
	
	// Check for DC bias (should be around 127.5 for 8-bit unsigned)
	expectedDC := 127.5
	if iAvg > expectedDC-20 && iAvg < expectedDC+20 && 
	   qAvg > expectedDC-20 && qAvg < expectedDC+20 {
		fmt.Printf("✓ DC bias appears normal (I: %.1f, Q: %.1f)\n", iAvg, qAvg)
	} else {
		fmt.Printf("⚠ DC bias may be off (I: %.1f, Q: %.1f, expected ~%.1f)\n", iAvg, qAvg, expectedDC)
	}
	
	// Check for all zeros (dead receiver)
	allZeros := true
	for _, b := range samples[:1000] { // Check first 1000 bytes
		if b != 0 {
			allZeros = false
			break
		}
	}
	
	if allZeros {
		fmt.Printf("✗ Data appears to be all zeros - receiver may not be working\n")
	} else {
		fmt.Printf("✓ Data contains non-zero values\n")
	}
	
	return nil
}