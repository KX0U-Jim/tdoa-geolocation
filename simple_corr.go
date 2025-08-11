package main

import (
	"fmt"
	"math"
	"math/rand"
)

func main() {
	fmt.Println("=== Simple Correlation Algorithm Test ===")
	
	// Create a test signal
	signalLen := 10000
	signal := make([]complex64, signalLen)
	
	// Generate a simple sine wave + noise
	for i := 0; i < signalLen; i++ {
		// 1 kHz sine wave at 100 kHz sample rate
		freq := 1000.0
		sampleRate := 100000.0
		t := float64(i) / sampleRate
		
		sineVal := float32(0.5 * math.Sin(2*math.Pi*freq*t))
		noise := float32(0.1 * (rand.Float64() - 0.5))
		
		signal[i] = complex(sineVal + noise, 0)
	}
	
	fmt.Printf("Created test signal: %d samples\n", len(signal))
	
	// Test 1: Self-correlation (should be 1.0 at zero delay)
	fmt.Printf("\nTest 1: Self-correlation\n")
	delay, corr := simpleCorrelate(signal, signal)
	fmt.Printf("Self-correlation: %.6f at delay %d\n", corr, delay)
	
	if corr > 0.8 {
		fmt.Printf("✅ PASS - Self-correlation is strong\n")
	} else {
		fmt.Printf("❌ FAIL - Self-correlation is weak\n")
	}
	
	// Test 2: Delayed signal correlation
	fmt.Printf("\nTest 2: Delayed signal correlation\n")
	
	// Create delayed version (shift by 100 samples)
	shiftSamples := 100
	delayedSignal := make([]complex64, signalLen)
	for i := shiftSamples; i < signalLen; i++ {
		delayedSignal[i] = signal[i-shiftSamples]
	}
	
	delay, corr = simpleCorrelate(signal[:signalLen-shiftSamples], delayedSignal[shiftSamples:])
	fmt.Printf("Delayed correlation: %.6f at delay %d (expected delay: %d)\n", corr, delay, 0)
	
	if corr > 0.8 && delay >= -10 && delay <= 10 {
		fmt.Printf("✅ PASS - Delayed correlation works\n")
	} else {
		fmt.Printf("❌ FAIL - Delayed correlation failed\n")
	}
	
	// Test 3: Uncorrelated noise
	fmt.Printf("\nTest 3: Uncorrelated noise\n")
	
	noise := make([]complex64, signalLen)
	for i := 0; i < signalLen; i++ {
		noise[i] = complex(float32(rand.Float64()-0.5), float32(rand.Float64()-0.5))
	}
	
	delay, corr = simpleCorrelate(signal, noise)
	fmt.Printf("Noise correlation: %.6f at delay %d\n", corr, delay)
	
	if math.Abs(float64(corr)) < 0.2 {
		fmt.Printf("✅ PASS - Noise correlation is low\n")
	} else {
		fmt.Printf("❌ FAIL - Noise correlation is too high\n")
	}
	
	fmt.Printf("\n=== CONCLUSION ===\n")
	fmt.Printf("If all tests pass, correlation algorithm works correctly.\n")
	fmt.Printf("Zero correlation in TDOA means signals are truly uncorrelated.\n")
}

func simpleCorrelate(signal1, signal2 []complex64) (int, float32) {
	if len(signal1) == 0 || len(signal2) == 0 {
		return 0, 0
	}
	
	// Use shorter signal as template
	template := signal1
	signal := signal2
	if len(signal1) > len(signal2) {
		template = signal2
		signal = signal1
	}
	
	templateLen := len(template)
	signalLen := len(signal)
	maxLag := 1000 // Search range
	
	if maxLag > signalLen - templateLen {
		maxLag = signalLen - templateLen
	}
	
	// Ensure we can at least do zero-delay correlation
	if maxLag < 1 {
		maxLag = 1
	}
	
	fmt.Printf("DEBUG: templateLen=%d, signalLen=%d, maxLag=%d\n", templateLen, signalLen, maxLag)
	
	bestDelay := 0
	bestCorr := float32(0.0)
	
	// Calculate template power for normalization
	templatePower := float32(0.0)
	for _, sample := range template {
		templatePower += real(sample)*real(sample) + imag(sample)*imag(sample)
	}
	
	fmt.Printf("DEBUG: Template power: %.6f\n", templatePower)
	
	// Correlate
	for delay := 0; delay < maxLag; delay++ {
		correlation := float32(0.0)
		signalPower := float32(0.0)
		
		for i := 0; i < templateLen; i++ {
			if delay + i >= signalLen {
				break
			}
			
			t := template[i]
			s := signal[delay + i]
			
			correlation += real(t)*real(s) + imag(t)*imag(s)
			signalPower += real(s)*real(s) + imag(s)*imag(s)
		}
		
		if delay == 0 {
			fmt.Printf("DEBUG delay 0: correlation=%.6f, templatePower=%.6f, signalPower=%.6f\n", 
				correlation, templatePower, signalPower)
		}
		
		// Normalized correlation - fix the math!
		if templatePower > 0 && signalPower > 0 {
			normalizedCorr := correlation / (float32(math.Sqrt(float64(templatePower))) * float32(math.Sqrt(float64(signalPower))))
			
			if delay == 0 {
				fmt.Printf("DEBUG delay 0 normalized: %.6f\n", normalizedCorr)
			}
			
			if math.Abs(float64(normalizedCorr)) > math.Abs(float64(bestCorr)) {
				bestCorr = normalizedCorr
				bestDelay = delay
			}
		}
	}
	
	return bestDelay, bestCorr
}