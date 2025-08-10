package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"sort"
	"strconv"
)

const (
	sampleRate      = 2000000.0 // 2 Msps
	bytesPerSample  = 2         // I + Q bytes
	expectedDCBias  = 127.5     // Expected DC for 8-bit unsigned
)

type SignalAnalysis struct {
	// Basic stats
	TotalSamples   int
	IAvg, QAvg     float64
	IMin, IMax     byte
	QMin, QMax     byte
	IStdDev, QStdDev float64
	
	// Signal quality metrics
	SNREstimate    float64
	PowerLevel     float64
	DCOffset       float64
	IQImbalance    float64
	
	// Frequency analysis
	PeakFrequency  float64
	BandwidthUsed  float64
	SpectralPurity float64
	
	// Quality flags
	HasClipping    bool
	HasOverload    bool
	HasDeadZones   bool
	HasNoise       bool
}

func main() {
	if len(os.Args) < 2 {
		fmt.Printf("Usage: %s <data_file.dat> [expected_duration_seconds]\n", os.Args[0])
		fmt.Printf("Advanced signal quality analyzer with gain and parameter recommendations\n")
		fmt.Printf("Example: %s kx0u-1723234567.dat 30\n", os.Args[0])
		os.Exit(1)
	}

	filename := os.Args[1]
	expectedDuration := 30 // default
	if len(os.Args) > 2 {
		if dur, err := strconv.Atoi(os.Args[2]); err == nil {
			expectedDuration = dur
		}
	}

	fmt.Printf("=== Advanced Signal Quality Analysis ===\n")
	fmt.Printf("File: %s\n", filename)
	fmt.Printf("Expected Duration: %d seconds\n\n", expectedDuration)

	// Perform separate analysis for reference and target signals
	refAnalysis, targetAnalysis, err := analyzeDualFrequencyFile(filename, expectedDuration)
	if err != nil {
		fmt.Printf("Error analyzing file: %v\n", err)
		os.Exit(1)
	}

	// Print analysis results
	fmt.Printf("=== REFERENCE SIGNAL ANALYSIS ===\n")
	printAnalysisResults(refAnalysis)
	generateRecommendations(refAnalysis, filename, "Reference")
	
	fmt.Printf("\n=== TARGET SIGNAL ANALYSIS ===\n")
	printAnalysisResults(targetAnalysis)
	generateRecommendations(targetAnalysis, filename, "Target")
	
	// Compare signals
	compareSignals(refAnalysis, targetAnalysis)
}

func analyzeDualFrequencyFile(filename string, expectedDuration int) (*SignalAnalysis, *SignalAnalysis, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, fmt.Errorf("failed to open file: %v", err)
	}
	defer file.Close()

	fileInfo, err := file.Stat()
	if err != nil {
		return nil, nil, fmt.Errorf("failed to get file info: %v", err)
	}

	totalBytes := fileInfo.Size()
	totalSamples := int(totalBytes / bytesPerSample)

	// Read all samples
	samples := make([]byte, totalBytes)
	_, err = file.Read(samples)
	if err != nil {
		return nil, nil, fmt.Errorf("failed to read samples: %v", err)
	}

	// librtlsdr-2freq pattern: 3×n samples (freq1, freq2, freq1)
	samplesPerBlock := totalSamples / 3
	
	fmt.Printf("=== File Structure Analysis ===\n")
	fmt.Printf("Total samples: %d\n", totalSamples)
	fmt.Printf("Samples per frequency block: %d\n", samplesPerBlock)
	fmt.Printf("Reference samples: %d (blocks 1+3)\n", samplesPerBlock*2)
	fmt.Printf("Target samples: %d (block 2)\n\n", samplesPerBlock)

	// Extract reference signal samples (blocks 1 and 3)
	refSamples := make([]byte, samplesPerBlock*2*bytesPerSample)
	copy(refSamples[0:samplesPerBlock*bytesPerSample], samples[0:samplesPerBlock*bytesPerSample])
	copy(refSamples[samplesPerBlock*bytesPerSample:], samples[samplesPerBlock*2*bytesPerSample:])

	// Extract target signal samples (block 2)
	targetSamples := samples[samplesPerBlock*bytesPerSample : samplesPerBlock*2*bytesPerSample]

	// Analyze both signals
	refAnalysis := analyzeSamples(refSamples, samplesPerBlock*2, "Reference")
	targetAnalysis := analyzeSamples(targetSamples, samplesPerBlock, "Target")

	return refAnalysis, targetAnalysis, nil
}

func analyzeSamples(samples []byte, totalSamples int, signalType string) *SignalAnalysis {
	analysis := &SignalAnalysis{
		TotalSamples: totalSamples,
		IMin: 255, IMax: 0,
		QMin: 255, QMax: 0,
	}

	// Basic statistics
	var iSum, qSum float64
	var iSumSq, qSumSq float64

	for i := 0; i < totalSamples; i++ {
		iVal := float64(samples[i*2])
		qVal := float64(samples[i*2+1])

		iSum += iVal
		qSum += qVal
		iSumSq += iVal * iVal
		qSumSq += qVal * qVal

		if samples[i*2] < analysis.IMin { analysis.IMin = samples[i*2] }
		if samples[i*2] > analysis.IMax { analysis.IMax = samples[i*2] }
		if samples[i*2+1] < analysis.QMin { analysis.QMin = samples[i*2+1] }
		if samples[i*2+1] > analysis.QMax { analysis.QMax = samples[i*2+1] }
	}

	// Calculate averages and standard deviations
	n := float64(totalSamples)
	analysis.IAvg = iSum / n
	analysis.QAvg = qSum / n
	analysis.IStdDev = math.Sqrt((iSumSq/n) - (analysis.IAvg * analysis.IAvg))
	analysis.QStdDev = math.Sqrt((qSumSq/n) - (analysis.QAvg * analysis.QAvg))

	// DC offset calculation
	analysis.DCOffset = math.Sqrt(math.Pow(analysis.IAvg-expectedDCBias, 2) + 
								  math.Pow(analysis.QAvg-expectedDCBias, 2))

	// IQ imbalance
	analysis.IQImbalance = math.Abs(analysis.IStdDev - analysis.QStdDev) / 
						   math.Max(analysis.IStdDev, analysis.QStdDev)

	// Power level calculation (RMS)
	analysis.PowerLevel = 20 * math.Log10(math.Sqrt(analysis.IStdDev*analysis.IStdDev + 
													 analysis.QStdDev*analysis.QStdDev))

	// Quality flags
	analysis.HasClipping = (analysis.IMin == 0 || analysis.IMax == 255 || 
							analysis.QMin == 0 || analysis.QMax == 255)
	
	analysis.HasOverload = (analysis.IStdDev < 2 || analysis.QStdDev < 2) // Too little variation
	
	analysis.HasDeadZones = checkForDeadZones(samples, totalSamples)
	
	analysis.HasNoise = analysis.IStdDev > 60 || analysis.QStdDev > 60 // Excessive variation

	// Calculate proper SNR using FFT-based frequency domain analysis
	analysis.SNREstimate = calculateProperSNR(samples, totalSamples)

	// Spectral analysis (simplified)
	analysis.PeakFrequency, analysis.BandwidthUsed, analysis.SpectralPurity = 
		analyzeSpectrum(samples, totalSamples)

	return analysis
}

func checkForDeadZones(samples []byte, totalSamples int) bool {
	consecutiveZeros := 0
	maxZeros := 0
	
	for i := 0; i < totalSamples*2; i++ {
		if samples[i] == 0 {
			consecutiveZeros++
		} else {
			if consecutiveZeros > maxZeros {
				maxZeros = consecutiveZeros
			}
			consecutiveZeros = 0
		}
	}
	
	return maxZeros > 1000 // More than 1000 consecutive zero samples
}

func calculateProperSNR(samples []byte, totalSamples int) float64 {
	// Use a reasonable sample size for analysis (max 16384 samples = 8ms at 2Msps)
	analysisSize := 16384
	if totalSamples < analysisSize {
		analysisSize = totalSamples
	}
	
	// Take samples from middle of capture to avoid startup transients
	startIdx := (totalSamples - analysisSize) / 2
	analysisBytes := samples[startIdx*2 : (startIdx+analysisSize)*2]
	
	// Convert samples to complex numbers with DC offset and IQ correction
	complexSamples := preprocessSamples(analysisBytes, analysisSize)
	
	// Apply Blackman-Harris window to reduce spectral leakage
	windowed := applyBlackmanHarrisWindow(complexSamples)
	
	// Compute FFT using DFT 
	fft := computeDFT(windowed)
	
	// Compute power spectral density
	psd := make([]float64, len(fft))
	for i, c := range fft {
		psd[i] = cmplx.Abs(c) * cmplx.Abs(c)
	}
	
	// Use Welch's method approach: separate signal from noise
	// Sort PSD values to find percentile thresholds
	sortedPSD := make([]float64, len(psd))
	copy(sortedPSD, psd)
	sort.Float64s(sortedPSD)
	
	// Signal power: top 10% of spectral bins (where signal energy concentrates)
	signalThresholdIdx := int(0.9 * float64(len(sortedPSD)))
	signalPower := 0.0
	signalBins := 0
	for _, power := range psd {
		if power >= sortedPSD[signalThresholdIdx] {
			signalPower += power
			signalBins++
		}
	}
	if signalBins > 0 {
		signalPower /= float64(signalBins)
	}
	
	// Noise power: bottom 50% of spectral bins (noise floor)
	noiseEndIdx := int(0.5 * float64(len(sortedPSD)))
	noisePower := 0.0
	for i := 0; i < noiseEndIdx; i++ {
		noisePower += sortedPSD[i]
	}
	noisePower /= float64(noiseEndIdx)
	
	// Calculate SNR in dB
	if noisePower > 0 && signalPower > noisePower {
		return 10 * math.Log10(signalPower/noisePower)
	}
	
	// Fallback for edge cases
	return -20.0
}

func preprocessSamples(samples []byte, totalSamples int) []complex128 {
	// Convert uint8 samples to complex128 with DC offset correction
	complexSamples := make([]complex128, totalSamples)
	
	// First pass: calculate DC offset
	var iSum, qSum float64
	for i := 0; i < totalSamples; i++ {
		iSum += float64(samples[i*2])
		qSum += float64(samples[i*2+1])
	}
	iDC := iSum / float64(totalSamples)
	qDC := qSum / float64(totalSamples)
	
	// Second pass: apply DC correction and convert to complex
	for i := 0; i < totalSamples; i++ {
		// Convert from uint8 [0,255] to float64 [-1,1] with DC correction
		iVal := (float64(samples[i*2]) - iDC) / 127.5
		qVal := (float64(samples[i*2+1]) - qDC) / 127.5
		complexSamples[i] = complex(iVal, qVal)
	}
	
	return complexSamples
}

func applyBlackmanHarrisWindow(samples []complex128) []complex128 {
	// Apply Blackman-Harris window for better spectral analysis
	n := len(samples)
	windowed := make([]complex128, n)
	
	for i := 0; i < n; i++ {
		// Blackman-Harris coefficients
		a0 := 0.35875
		a1 := 0.48829
		a2 := 0.14128
		a3 := 0.01168
		
		w := a0 - a1*math.Cos(2*math.Pi*float64(i)/float64(n-1)) +
			 a2*math.Cos(4*math.Pi*float64(i)/float64(n-1)) -
			 a3*math.Cos(6*math.Pi*float64(i)/float64(n-1))
		
		windowed[i] = complex(w, 0) * samples[i]
	}
	
	return windowed
}

func computeDFT(samples []complex128) []complex128 {
	// Simple DFT implementation (not optimized FFT, but correct)
	n := len(samples)
	dft := make([]complex128, n)
	
	for k := 0; k < n; k++ {
		sum := complex(0, 0)
		for i := 0; i < n; i++ {
			angle := -2 * math.Pi * float64(k) * float64(i) / float64(n)
			sum += samples[i] * complex(math.Cos(angle), math.Sin(angle))
		}
		dft[k] = sum
	}
	
	return dft
}

func analyzeSpectrum(samples []byte, totalSamples int) (peakFreq, bandwidth, purity float64) {
	// Simplified spectral analysis - in real implementation would use FFT
	// For now, estimate based on signal characteristics
	
	peakFreq = 0.0 // Would need FFT to determine actual peak frequency
	bandwidth = 2000000.0 // Full RTL-SDR bandwidth
	purity = 0.8 // Placeholder - would calculate from actual spectrum
	
	return peakFreq, bandwidth, purity
}

func printAnalysisResults(analysis *SignalAnalysis) {
	fmt.Printf("=== Signal Statistics ===\n")
	fmt.Printf("Total Samples: %d\n", analysis.TotalSamples)
	fmt.Printf("I Channel: min=%d, max=%d, avg=%.1f, σ=%.1f\n", 
		analysis.IMin, analysis.IMax, analysis.IAvg, analysis.IStdDev)
	fmt.Printf("Q Channel: min=%d, max=%d, avg=%.1f, σ=%.1f\n", 
		analysis.QMin, analysis.QMax, analysis.QAvg, analysis.QStdDev)
	
	fmt.Printf("\n=== Signal Quality Metrics ===\n")
	fmt.Printf("DC Offset: %.1f (should be ~0)\n", analysis.DCOffset)
	fmt.Printf("IQ Imbalance: %.3f (should be <0.1)\n", analysis.IQImbalance)
	fmt.Printf("Estimated SNR: %.1f dB\n", analysis.SNREstimate)
	fmt.Printf("Power Level: %.1f dB\n", analysis.PowerLevel)
	
	fmt.Printf("\n=== Quality Flags ===\n")
	printFlag("Clipping/Saturation", analysis.HasClipping, "⚠️", "✅")
	printFlag("Overload (too low variation)", analysis.HasOverload, "⚠️", "✅")
	printFlag("Dead zones detected", analysis.HasDeadZones, "⚠️", "✅")
	printFlag("Excessive noise", analysis.HasNoise, "⚠️", "✅")
}

func printFlag(name string, condition bool, failIcon, passIcon string) {
	if condition {
		fmt.Printf("%s %s: DETECTED\n", failIcon, name)
	} else {
		fmt.Printf("%s %s: OK\n", passIcon, name)
	}
}

func generateRecommendations(analysis *SignalAnalysis, filename string, signalType string) {
	fmt.Printf("\n=== %s SIGNAL RECOMMENDATIONS ===\n", signalType)
	
	// Gain recommendations
	generateGainRecommendations(analysis)
	
	// Hardware recommendations  
	generateHardwareRecommendations(analysis)
	
	// Collection parameter recommendations
	generateCollectionRecommendations(analysis)
	
	// System recommendations
	generateSystemRecommendations(analysis)
	
	// Future enhancement recommendations
	generateEnhancementRecommendations(analysis)
}

func compareSignals(refAnalysis, targetAnalysis *SignalAnalysis) {
	fmt.Printf("\n=== SIGNAL COMPARISON ===\n")
	
	// SNR Comparison
	fmt.Printf("SNR Comparison:\n")
	fmt.Printf("  Reference: %.1f dB\n", refAnalysis.SNREstimate)
	fmt.Printf("  Target:    %.1f dB\n", targetAnalysis.SNREstimate)
	if refAnalysis.SNREstimate > targetAnalysis.SNREstimate+10 {
		fmt.Printf("  ⚠️  Reference significantly stronger - consider reducing reference gain\n")
	} else if targetAnalysis.SNREstimate > refAnalysis.SNREstimate+10 {
		fmt.Printf("  ⚠️  Target significantly stronger - consider reducing target gain\n")
	} else {
		fmt.Printf("  ✅ Signal levels reasonably balanced\n")
	}
	
	// Power Level Comparison
	fmt.Printf("\nPower Level Comparison:\n")
	fmt.Printf("  Reference: %.1f dB\n", refAnalysis.PowerLevel)
	fmt.Printf("  Target:    %.1f dB\n", targetAnalysis.PowerLevel)
	
	// Quality Issues Comparison
	fmt.Printf("\nQuality Issues:\n")
	refIssues := countQualityIssues(refAnalysis)
	targetIssues := countQualityIssues(targetAnalysis)
	fmt.Printf("  Reference: %d issues detected\n", refIssues)
	fmt.Printf("  Target:    %d issues detected\n", targetIssues)
	
	if refIssues == 0 && targetIssues == 0 {
		fmt.Printf("  ✅ Both signals appear suitable for TDOA processing\n")
	} else if refIssues > targetIssues {
		fmt.Printf("  ⚠️  Reference signal needs more attention\n")
	} else if targetIssues > refIssues {
		fmt.Printf("  ⚠️  Target signal needs more attention\n")
	}
	
	// TDOA Suitability Assessment
	fmt.Printf("\n=== TDOA SUITABILITY ASSESSMENT ===\n")
	refSuitable := assessTDOASuitability(refAnalysis)
	targetSuitable := assessTDOASuitability(targetAnalysis)
	
	if refSuitable && targetSuitable {
		fmt.Printf("✅ EXCELLENT: Both signals suitable for TDOA correlation\n")
	} else if !refSuitable && !targetSuitable {
		fmt.Printf("❌ POOR: Both signals need improvement before TDOA processing\n")
	} else if !refSuitable {
		fmt.Printf("⚠️  MARGINAL: Reference signal needs improvement\n")
	} else {
		fmt.Printf("⚠️  MARGINAL: Target signal needs improvement\n")
	}
}

func countQualityIssues(analysis *SignalAnalysis) int {
	issues := 0
	if analysis.HasClipping { issues++ }
	if analysis.HasOverload { issues++ }
	if analysis.HasDeadZones { issues++ }
	if analysis.HasNoise { issues++ }
	if analysis.DCOffset > 10 { issues++ }
	if analysis.IQImbalance > 0.1 { issues++ }
	return issues
}

func assessTDOASuitability(analysis *SignalAnalysis) bool {
	if analysis.HasClipping || analysis.HasOverload || analysis.HasDeadZones {
		return false
	}
	if analysis.SNREstimate < 15 {
		return false
	}
	if analysis.DCOffset > 15 || analysis.IQImbalance > 0.15 {
		return false
	}
	return true
}

func generateGainRecommendations(analysis *SignalAnalysis) {
	fmt.Printf("\n--- Gain Recommendations ---\n")
	
	if analysis.HasClipping {
		fmt.Printf("🔻 REDUCE GAIN: Signal clipping detected\n")
		fmt.Printf("   Try --gain=10 to --gain=30 (reduce by 10-20 dB)\n")
		fmt.Printf("   Clipping causes distortion and poor correlation\n")
	} else if analysis.HasOverload {
		fmt.Printf("🔻 REDUCE GAIN: Signal appears overloaded\n")
		fmt.Printf("   Try --gain=20 to --gain=40 (reduce by 5-15 dB)\n")
		fmt.Printf("   Very strong signals can saturate the ADC\n")
	} else if analysis.PowerLevel < -60 {
		fmt.Printf("🔺 INCREASE GAIN: Signal level very low\n")
		fmt.Printf("   Try --gain=40 to --gain=49.6 (increase by 10-20 dB)\n")
		fmt.Printf("   Low signals will have poor SNR for correlation\n")
	} else if analysis.PowerLevel < -40 {
		fmt.Printf("🔺 INCREASE GAIN: Signal level low\n")
		fmt.Printf("   Try --gain=35 to --gain=45 (increase by 5-10 dB)\n")
		fmt.Printf("   Optimize for best SNR without overload\n")
	} else if analysis.IStdDev > 50 && analysis.QStdDev > 50 && !analysis.HasClipping {
		fmt.Printf("✅ GAIN OK: Good signal level, no clipping\n")
		fmt.Printf("   Current gain appears optimal\n")
	} else {
		fmt.Printf("🔧 FINE-TUNE GAIN: Signal usable but could be optimized\n")
		fmt.Printf("   Try adjusting gain by ±5 dB and compare results\n")
	}
	
	// SNR-based recommendations
	if analysis.SNREstimate < 10 {
		fmt.Printf("📡 SNR TOO LOW (%.1f dB): Increase gain or improve antenna\n", analysis.SNREstimate)
	} else if analysis.SNREstimate > 40 {
		fmt.Printf("📡 SNR HIGH (%.1f dB): Consider reducing gain to prevent overload\n", analysis.SNREstimate)
	}
}

func generateHardwareRecommendations(analysis *SignalAnalysis) {
	fmt.Printf("\n--- Hardware Recommendations ---\n")
	
	if analysis.DCOffset > 10 {
		fmt.Printf("⚡ DC OFFSET ISSUE: DC bias = %.1f (should be ~0)\n", analysis.DCOffset)
		fmt.Printf("   • Check RTL-SDR connection and USB power\n")
		fmt.Printf("   • Try different USB port or powered USB hub\n")
		fmt.Printf("   • RTL-SDR may need DC offset correction in software\n")
		fmt.Printf("   • Future: Add --dc-offset parameter to collector\n")
	}
	
	if analysis.IQImbalance > 0.1 {
		fmt.Printf("⚖️ IQ IMBALANCE: %.3f (should be <0.1)\n", analysis.IQImbalance)
		fmt.Printf("   • RTL-SDR hardware calibration issue\n")
		fmt.Printf("   • Try different RTL-SDR dongle if available\n")
		fmt.Printf("   • Future: Add IQ balance correction in software\n")
	}
	
	if analysis.HasDeadZones {
		fmt.Printf("☠️ DEAD ZONES: Extended periods of zero samples\n")
		fmt.Printf("   • USB bandwidth issue - try different USB port\n")
		fmt.Printf("   • System overload - close other applications\n")
		fmt.Printf("   • Sample rate too high for system capabilities\n")
		fmt.Printf("   • Future: Add --sample-rate parameter\n")
	}
	
	if analysis.HasNoise {
		fmt.Printf("🌩️ EXCESSIVE NOISE: High signal variation\n")
		fmt.Printf("   • RF interference - check antenna placement\n")
		fmt.Printf("   • Poor antenna system - improve antenna/feedline\n")
		fmt.Printf("   • EMI from computers/switching supplies nearby\n")
	}
}

func generateCollectionRecommendations(analysis *SignalAnalysis) {
	fmt.Printf("\n--- Collection Parameter Recommendations ---\n")
	
	// Sample rate recommendations
	fmt.Printf("📊 SAMPLE RATE: Currently 2.0 Msps\n")
	if analysis.BandwidthUsed < 500000 {
		fmt.Printf("   • Consider --sample-rate=1200000 for narrowband signals\n")
		fmt.Printf("   • Reduces data size and processing load\n")
	} else if analysis.BandwidthUsed > 1800000 {
		fmt.Printf("   • May need --sample-rate=2400000 for wideband signals\n")
		fmt.Printf("   • Check if signal is being aliased\n")
	}
	
	// Duration recommendations
	fmt.Printf("⏱️ DURATION: Based on signal stability\n")
	if analysis.SNREstimate < 15 {
		fmt.Printf("   • Use longer duration (60-100s) for weak signals\n")
		fmt.Printf("   • More data improves correlation reliability\n")
	} else if analysis.SNREstimate > 30 {
		fmt.Printf("   • Can use shorter duration (10-30s) for strong signals\n")
		fmt.Printf("   • Reduces chance of timing drift\n")
	}
	
	// Frequency recommendations
	fmt.Printf("📻 FREQUENCY: Check signal characteristics\n")
	if analysis.PowerLevel < -50 {
		fmt.Printf("   • Verify target frequency is correct\n")
		fmt.Printf("   • Signal may be off-frequency or not present\n")
		fmt.Printf("   • Future: Add frequency scanning capability\n")
	}
}

func generateSystemRecommendations(analysis *SignalAnalysis) {
	fmt.Printf("\n--- System Optimization ---\n")
	
	fmt.Printf("💻 PERFORMANCE:\n")
	fmt.Printf("   • Close unnecessary applications during collection\n")
	fmt.Printf("   • Use wired Ethernet instead of WiFi if possible\n")
	fmt.Printf("   • Ensure adequate free disk space (>1GB)\n")
	
	if analysis.HasDeadZones {
		fmt.Printf("🚀 USB OPTIMIZATION:\n")
		fmt.Printf("   • Use USB 3.0 port if available\n")
		fmt.Printf("   • Try powered USB hub to improve signal quality\n")
		fmt.Printf("   • Avoid USB hubs shared with other devices\n")
	}
	
	fmt.Printf("🔧 RTLSDR OPTIMIZATION:\n")
	fmt.Printf("   • Future: Add --ppm-error parameter for frequency correction\n")
	fmt.Printf("   • Future: Add --bandwidth parameter for tuner filtering\n")
	fmt.Printf("   • Future: Add --buffer-size parameter for USB transfers\n")
}

func generateEnhancementRecommendations(analysis *SignalAnalysis) {
	fmt.Printf("\n--- Future Enhancements Needed ---\n")
	
	fmt.Printf("🔮 COLLECTOR PARAMETERS TO ADD:\n")
	fmt.Printf("   • --sample-rate=Hz (currently fixed at 2 Msps)\n")
	fmt.Printf("   • --ppm-error=float (frequency correction)\n")
	fmt.Printf("   • --dc-offset (DC bias correction)\n")
	fmt.Printf("   • --iq-balance (IQ imbalance correction)\n")
	fmt.Printf("   • --buffer-size=bytes (USB transfer optimization)\n")
	fmt.Printf("   • --bandwidth=Hz (anti-aliasing filter)\n")
	
	fmt.Printf("📈 ANALYSIS IMPROVEMENTS:\n")
	fmt.Printf("   • Real FFT spectral analysis\n")
	fmt.Printf("   • Automatic gain control recommendations\n")
	fmt.Printf("   • Signal type detection (FM/AM/Digital)\n")
	fmt.Printf("   • Multi-frequency block analysis\n")
	fmt.Printf("   • Correlation preview between stations\n")
	
	fmt.Printf("🎯 QUALITY METRICS:\n")
	fmt.Printf("   • Phase noise analysis\n")
	fmt.Printf("   • Spurious signal detection\n")
	fmt.Printf("   • Image rejection measurement\n")
	fmt.Printf("   • Dynamic range assessment\n")
	
	fmt.Printf("\n=== SUMMARY ===\n")
	if analysis.HasClipping || analysis.HasOverload {
		fmt.Printf("❌ CRITICAL: Adjust gain immediately - signal distortion present\n")
	} else if analysis.PowerLevel < -50 {
		fmt.Printf("⚠️  WARNING: Signal very weak - increase gain or check antenna\n")
	} else if analysis.DCOffset > 10 || analysis.IQImbalance > 0.1 {
		fmt.Printf("🔧 HARDWARE: RTL-SDR calibration issues detected\n")
	} else {
		fmt.Printf("✅ ACCEPTABLE: Signal quality adequate for TDOA processing\n")
	}
}