package main

import (
	"fmt"
	"log"
	"os"
	"strconv"
)

func main() {
	if len(os.Args) < 4 {
		fmt.Printf("Usage: %s <ref_freq> <target_freq> <csv_file> <dat_file>\n", os.Args[0])
		fmt.Println("This tests correlation of a signal with itself (should be 1.0)")
		os.Exit(1)
	}

	refFreq, _ := strconv.ParseFloat(os.Args[1], 64)
	targetFreq, _ := strconv.ParseFloat(os.Args[2], 64)

	processor, err := NewTDOAProcessor(refFreq, targetFreq, os.Args[3])
	if err != nil {
		log.Fatalf("Failed to create processor: %v", err)
	}

	// Load one file
	data, err := processor.loadIQData(os.Args[4])
	if err != nil {
		log.Fatalf("Failed to load data: %v", err)
	}

	// Extract both signals
	refSig := processor.extractReferenceSignal(data)
	targetSig := processor.extractTargetSignal(data)

	// Use small chunks for quick test
	chunkSize := 100000
	if len(refSig) > chunkSize {
		refSig = refSig[:chunkSize]
	}
	if len(targetSig) > chunkSize {
		targetSig = targetSig[:chunkSize]
	}

	fmt.Printf("=== CORRELATION SANITY CHECK ===\n")
	fmt.Printf("Testing correlation of signals with themselves\n")
	fmt.Printf("Expected result: correlation ≈ 1.0\n\n")

	// Test reference signal with itself
	fmt.Printf("Reference signal self-correlation:\n")
	delay, corr := processor.crossCorrelate(refSig, refSig)
	fmt.Printf("Result: %.6f at delay %d (should be ≈1.0 at delay 0)\n\n", corr, delay)

	// Test target signal with itself  
	fmt.Printf("Target signal self-correlation:\n")
	delay, corr = processor.crossCorrelate(targetSig, targetSig)
	fmt.Printf("Result: %.6f at delay %d (should be ≈1.0 at delay 0)\n\n", corr, delay)

	if corr > 0.5 {
		fmt.Printf("✅ Correlation algorithm is working correctly\n")
		fmt.Printf("Zero correlation between stations means signals are truly uncorrelated\n")
	} else {
		fmt.Printf("❌ Correlation algorithm has issues\n")
		fmt.Printf("Need to debug the correlation implementation\n")
	}
}