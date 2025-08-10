package main

import (
	"fmt"
	"os"
	"os/exec"
	"strconv"
	"strings"
	"time"
)

const (
	minGain           = 5.0   // Minimum safe gain
	maxGain           = 45.0  // Maximum RTL-SDR gain
	targetSNR         = 25.0  // Target SNR in dB
	minAcceptableSNR  = 18.0  // Minimum acceptable SNR
	maxAcceptableSNR  = 40.0  // Maximum acceptable SNR (avoid unnecessary gain)
	convergenceTolerance = 2.0 // Stop when gain range < 2 dB
	testDuration      = 2     // Test collection duration in seconds
	maxIterations     = 8     // Safety limit on iterations
)

type CalibrationResult struct {
	Frequency    int64
	OptimalGain  float64
	AchievedSNR  float64
	HasClipping  bool
	HasOverload  bool
	PowerLevel   float64
	Iterations   int
	Success      bool
}

func main() {
	if len(os.Args) < 3 {
		fmt.Printf("Usage: %s <reference_freq_hz> <target_freq_hz>\n", os.Args[0])
		fmt.Printf("Automatic gain calibration for TDOA collectors\n")
		fmt.Printf("Example: %s 96900000 162550000\n", os.Args[0])
		fmt.Printf("\nOutputs optimal --gain1=X --gain2=Y for collector\n")
		os.Exit(1)
	}

	refFreq, err := strconv.ParseInt(os.Args[1], 10, 64)
	if err != nil {
		fmt.Printf("Error: Invalid reference frequency: %s\n", os.Args[1])
		os.Exit(1)
	}

	targetFreq, err := strconv.ParseInt(os.Args[2], 10, 64)
	if err != nil {
		fmt.Printf("Error: Invalid target frequency: %s\n", os.Args[2])
		os.Exit(1)
	}

	fmt.Printf("=== TDOA Automatic Gain Calibration ===\n")
	fmt.Printf("Reference Frequency: %.3f MHz\n", float64(refFreq)/1e6)
	fmt.Printf("Target Frequency:    %.3f MHz\n", float64(targetFreq)/1e6)
	fmt.Printf("Target SNR Range:    %.1f - %.1f dB\n", minAcceptableSNR, maxAcceptableSNR)
	fmt.Printf("Test Duration:       %d seconds per test\n", testDuration)
	fmt.Printf("\n")

	// Check required binaries
	if err := checkRequiredBinaries(); err != nil {
		fmt.Printf("Error: %v\n", err)
		os.Exit(1)
	}

	// Calibrate reference frequency
	fmt.Printf("=== Calibrating Reference Frequency (%.1f MHz) ===\n", float64(refFreq)/1e6)
	refResult := calibrateFrequency(refFreq, "ref")
	
	// Calibrate target frequency  
	fmt.Printf("\n=== Calibrating Target Frequency (%.1f MHz) ===\n", float64(targetFreq)/1e6)
	targetResult := calibrateFrequency(targetFreq, "target")

	// Print results
	printResults(refResult, targetResult)
}

func checkRequiredBinaries() error {
	binaries := []string{"./collector", "./fast_analyzer"}
	for _, binary := range binaries {
		if _, err := os.Stat(binary); os.IsNotExist(err) {
			return fmt.Errorf("%s not found", binary)
		}
	}
	return nil
}

func calibrateFrequency(freq int64, label string) CalibrationResult {
	result := CalibrationResult{
		Frequency: freq,
		Success:   false,
	}

	minG := minGain
	maxG := maxGain
	iteration := 0

	fmt.Printf("Starting binary search calibration...\n")

	for iteration < maxIterations && (maxG-minG) > convergenceTolerance {
		iteration++
		testGain := (minG + maxG) / 2.0

		fmt.Printf("Iteration %d: Testing gain %.1f dB (range: %.1f - %.1f)\n", 
			iteration, testGain, minG, maxG)

		// Run test collection
		analysis, err := runTestCollection(freq, testGain, label, iteration)
		if err != nil {
			fmt.Printf("  ERROR: %v\n", err)
			// On error, try slightly higher gain
			minG = testGain + 1
			continue
		}

		fmt.Printf("  Results: SNR=%.1f dB, Power=%.1f dB, Clipping=%t, Overload=%t\n",
			analysis.SNR, analysis.PowerLevel, analysis.HasClipping, analysis.HasOverload)

		// Decision logic
		if analysis.HasClipping {
			fmt.Printf("  → Clipping detected, reducing gain\n")
			maxG = testGain - 1.0
		} else if analysis.HasOverload {
			fmt.Printf("  → Overload detected, increasing gain\n")
			minG = testGain + 1.0
		} else if analysis.SNR < minAcceptableSNR {
			fmt.Printf("  → SNR too low (%.1f < %.1f), increasing gain\n", analysis.SNR, minAcceptableSNR)
			minG = testGain + 1.0
		} else if analysis.SNR > maxAcceptableSNR {
			fmt.Printf("  → SNR too high (%.1f > %.1f), reducing gain\n", analysis.SNR, maxAcceptableSNR)
			maxG = testGain - 1.0
		} else {
			fmt.Printf("  ✅ Optimal gain found!\n")
			result.OptimalGain = testGain
			result.AchievedSNR = analysis.SNR
			result.HasClipping = analysis.HasClipping
			result.HasOverload = analysis.HasOverload
			result.PowerLevel = analysis.PowerLevel
			result.Iterations = iteration
			result.Success = true
			return result
		}
	}

	// Convergence reached - pick best compromise
	finalGain := (minG + maxG) / 2.0
	fmt.Printf("Convergence reached after %d iterations\n", iteration)
	fmt.Printf("Final validation at gain %.1f dB...\n", finalGain)

	analysis, err := runTestCollection(freq, finalGain, label, iteration+1)
	if err != nil {
		fmt.Printf("Final validation failed: %v\n", err)
		// Return best guess
		result.OptimalGain = finalGain
		result.Iterations = iteration + 1
		return result
	}

	result.OptimalGain = finalGain
	result.AchievedSNR = analysis.SNR
	result.HasClipping = analysis.HasClipping
	result.HasOverload = analysis.HasOverload
	result.PowerLevel = analysis.PowerLevel
	result.Iterations = iteration + 1
	result.Success = !analysis.HasClipping && analysis.SNR >= minAcceptableSNR

	if result.Success {
		fmt.Printf("✅ Calibration successful\n")
	} else {
		fmt.Printf("⚠️  Calibration completed with compromise\n")
	}

	return result
}

type TestAnalysis struct {
	SNR         float64
	PowerLevel  float64
	HasClipping bool
	HasOverload bool
}

func runTestCollection(freq int64, gain float64, label string, iteration int) (*TestAnalysis, error) {
	// Generate unique filename
	timestamp := time.Now().Unix()
	filename := fmt.Sprintf("cal_%s_g%.0f_%d_%d.dat", label, gain, timestamp, iteration)

	// Calculate start time (5 seconds from now)
	startTime := time.Now().Unix() + 5

	// Run collector
	fmt.Printf("  Collecting %d seconds at gain %.1f dB...\n", testDuration, gain)
	
	cmd := exec.Command("./collector",
		fmt.Sprintf("--duration=%d", testDuration),
		fmt.Sprintf("--gain=%.1f", gain),
		fmt.Sprintf("%d", freq),
		fmt.Sprintf("%d", freq), // Same frequency for both (single freq test)
		fmt.Sprintf("%d", startTime),
		fmt.Sprintf("cal_%s_g%.0f", label, gain))

	output, err := cmd.CombinedOutput()
	if err != nil {
		return nil, fmt.Errorf("collector failed: %v\nOutput: %s", err, string(output))
	}

	// Find the actual output file
	actualFile, err := findCollectorOutput(filename, label, gain, startTime)
	if err != nil {
		return nil, fmt.Errorf("output file not found: %v", err)
	}

	// Run fast analyzer
	fmt.Printf("  Analyzing %s...\n", actualFile)
	
	cmd = exec.Command("./fast_analyzer", actualFile)
	output, err = cmd.Output()
	if err != nil {
		return nil, fmt.Errorf("analyzer failed: %v", err)
	}

	// Parse analyzer output
	analysis, err := parseAnalyzerOutput(string(output))
	if err != nil {
		return nil, fmt.Errorf("failed to parse analysis: %v", err)
	}

	// Clean up test file
	os.Remove(actualFile)

	return analysis, nil
}

func findCollectorOutput(expectedFile, label string, gain float64, startTime int64) (string, error) {
	// Try various possible filenames the collector might have used
	patterns := []string{
		expectedFile,
		fmt.Sprintf("cal_%s_g%.0f-*.dat", label, gain),
		fmt.Sprintf("cal_%s_g%.0f_%d.dat", label, gain, startTime),
	}

	for _, pattern := range patterns {
		if strings.Contains(pattern, "*") {
			// Use shell expansion for wildcard patterns
			cmd := exec.Command("sh", "-c", fmt.Sprintf("ls %s 2>/dev/null | head -1", pattern))
			output, err := cmd.Output()
			if err == nil && len(strings.TrimSpace(string(output))) > 0 {
				return strings.TrimSpace(string(output)), nil
			}
		} else {
			// Check exact filename
			if _, err := os.Stat(pattern); err == nil {
				return pattern, nil
			}
		}
	}

	return "", fmt.Errorf("no output file found matching patterns: %v", patterns)
}

func parseAnalyzerOutput(output string) (*TestAnalysis, error) {
	lines := strings.Split(output, "\n")
	
	// Find the first signal line (either REF or TGT)
	for _, line := range lines {
		line = strings.TrimSpace(line)
		if strings.HasPrefix(line, "REF,") || strings.HasPrefix(line, "TGT,") {
			parts := strings.Split(line, ",")
			if len(parts) < 5 {
				continue
			}

			snr, err1 := strconv.ParseFloat(parts[1], 64)
			power, err2 := strconv.ParseFloat(parts[2], 64)
			clipping, err3 := strconv.ParseBool(parts[3])
			overload, err4 := strconv.ParseBool(parts[4])

			if err1 != nil || err2 != nil || err3 != nil || err4 != nil {
				continue
			}

			return &TestAnalysis{
				SNR:         snr,
				PowerLevel:  power,
				HasClipping: clipping,
				HasOverload: overload,
			}, nil
		}
	}

	return nil, fmt.Errorf("no valid analysis data found in output: %s", output)
}

func printResults(refResult, targetResult CalibrationResult) {
	fmt.Printf("\n" + strings.Repeat("=", 60) + "\n")
	fmt.Printf("CALIBRATION RESULTS\n")
	fmt.Printf(strings.Repeat("=", 60) + "\n")

	// Reference results
	fmt.Printf("\nReference Frequency (%.1f MHz):\n", float64(refResult.Frequency)/1e6)
	fmt.Printf("  Optimal Gain:  %.1f dB\n", refResult.OptimalGain)
	fmt.Printf("  Achieved SNR:  %.1f dB\n", refResult.AchievedSNR)
	fmt.Printf("  Power Level:   %.1f dB\n", refResult.PowerLevel)
	fmt.Printf("  Iterations:    %d\n", refResult.Iterations)
	
	status := "✅ SUCCESS"
	if !refResult.Success {
		status = "⚠️ COMPROMISE"
		if refResult.HasClipping {
			status += " (clipping)"
		}
		if refResult.HasOverload {
			status += " (overload)"
		}
	}
	fmt.Printf("  Status:        %s\n", status)

	// Target results
	fmt.Printf("\nTarget Frequency (%.1f MHz):\n", float64(targetResult.Frequency)/1e6)
	fmt.Printf("  Optimal Gain:  %.1f dB\n", targetResult.OptimalGain)
	fmt.Printf("  Achieved SNR:  %.1f dB\n", targetResult.AchievedSNR)
	fmt.Printf("  Power Level:   %.1f dB\n", targetResult.PowerLevel)
	fmt.Printf("  Iterations:    %d\n", targetResult.Iterations)
	
	status = "✅ SUCCESS"
	if !targetResult.Success {
		status = "⚠️ COMPROMISE"
		if targetResult.HasClipping {
			status += " (clipping)"
		}
		if targetResult.HasOverload {
			status += " (overload)"
		}
	}
	fmt.Printf("  Status:        %s\n", status)

	// Collector command recommendation
	fmt.Printf("\n" + strings.Repeat("=", 60) + "\n")
	fmt.Printf("RECOMMENDED COLLECTOR COMMAND:\n")
	fmt.Printf(strings.Repeat("=", 60) + "\n")
	fmt.Printf("./collector --gain1=%.1f --gain2=%.1f %d %d <start_time> <station>\n", 
		refResult.OptimalGain, targetResult.OptimalGain,
		refResult.Frequency, targetResult.Frequency)

	// Quality assessment
	fmt.Printf("\nQUALITY ASSESSMENT:\n")
	overallSuccess := refResult.Success && targetResult.Success
	if overallSuccess {
		fmt.Printf("✅ Both frequencies calibrated successfully\n")
		fmt.Printf("   System ready for TDOA collection\n")
	} else {
		fmt.Printf("⚠️ One or more frequencies required compromise\n")
		if !refResult.Success {
			fmt.Printf("   Reference frequency: check antenna/signal strength\n")
		}
		if !targetResult.Success {
			fmt.Printf("   Target frequency: check antenna/signal strength\n")
		}
		fmt.Printf("   Consider improving antenna system for better results\n")
	}

	fmt.Printf("\nTotal calibration time: ~%d seconds\n", 
		(refResult.Iterations+targetResult.Iterations)*testDuration)
}