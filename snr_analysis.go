package main

import (
	"fmt"
	"math"
)

func main() {
	fmt.Println("=== SNR Analysis for TDOA Reference Signal ===")
	
	// Measured reference signal power levels from processor output
	fmt.Println("\nMeasured reference signal power levels:")
	kx0u_ref := 0.002718990
	n3pay_ref := 0.000075721
	kf0mtl_ref := 0.005146538
	
	fmt.Printf("kx0u:   %.9f\n", kx0u_ref)
	fmt.Printf("n3pay:  %.9f\n", n3pay_ref)
	fmt.Printf("kf0mtl: %.9f\n", kf0mtl_ref)
	
	// Convert to dB (relative to 1.0)
	fmt.Println("\nSignal power in dB (relative to full scale):")
	kx0u_db := 10 * math.Log10(kx0u_ref)
	n3pay_db := 10 * math.Log10(n3pay_ref)
	kf0mtl_db := 10 * math.Log10(kf0mtl_ref)
	
	fmt.Printf("kx0u:   %.1f dB\n", kx0u_db)
	fmt.Printf("n3pay:  %.1f dB\n", n3pay_db)
	fmt.Printf("kf0mtl: %.1f dB\n", kf0mtl_db)
	
	// For RTL-SDR, noise floor is typically around -50 to -60 dB
	noiseFloor := -55.0 // dB
	
	fmt.Printf("\nEstimated noise floor: %.1f dB\n", noiseFloor)
	
	// Calculate SNR
	fmt.Println("\nEstimated SNR:")
	fmt.Printf("kx0u:   %.1f - (%.1f) = %.1f dB\n", kx0u_db, noiseFloor, kx0u_db - noiseFloor)
	fmt.Printf("n3pay:  %.1f - (%.1f) = %.1f dB\n", n3pay_db, noiseFloor, n3pay_db - noiseFloor)
	fmt.Printf("kf0mtl: %.1f - (%.1f) = %.1f dB\n", kf0mtl_db, noiseFloor, kf0mtl_db - noiseFloor)
	
	fmt.Println("\n=== SNR Requirements for Different Applications ===")
	fmt.Println("Minimum SNR typically needed:")
	fmt.Println("- Voice communications (intelligible): 6-10 dB")
	fmt.Println("- Data communications (low error): 10-15 dB") 
	fmt.Println("- Precise timing/correlation: 15-20 dB")
	fmt.Println("- High-precision TDOA: 20-25 dB")
	fmt.Println("- Sub-sample timing accuracy: 25-30 dB")
	
	fmt.Println("\n=== Analysis ===")
	minCorrelationSNR := 15.0
	minPreciseTDOA := 20.0
	
	fmt.Printf("For basic correlation, need: %.1f dB SNR\n", minCorrelationSNR)
	fmt.Printf("For precise TDOA, need: %.1f dB SNR\n", minPreciseTDOA)
	
	stations := []string{"kx0u", "n3pay", "kf0mtl"}
	snrs := []float64{kx0u_db - noiseFloor, n3pay_db - noiseFloor, kf0mtl_db - noiseFloor}
	
	fmt.Println("\nStation analysis:")
	for i, station := range stations {
		status := "❌ TOO WEAK"
		if snrs[i] >= minPreciseTDOA {
			status = "✅ EXCELLENT"
		} else if snrs[i] >= minCorrelationSNR {
			status = "⚠️  MARGINAL"
		}
		
		fmt.Printf("%s: %.1f dB SNR - %s\n", station, snrs[i], status)
	}
	
	fmt.Println("\n=== Recommendations ===")
	if snrs[1] < minCorrelationSNR { // n3pay is weakest
		deficit := minCorrelationSNR - snrs[1]
		fmt.Printf("n3pay needs %.1f dB more signal strength\n", deficit)
		fmt.Printf("This could be achieved by:\n")
		fmt.Printf("- Better antenna (%.1f dB gain)\n", deficit)
		fmt.Printf("- Lower noise figure receiver\n")
		fmt.Printf("- Different reference frequency\n")
		fmt.Printf("- Signal processing gain (coherent integration)\n")
	}
	
	fmt.Println("\nCoherent integration gain:")
	integrationTimes := []float64{1, 10, 100, 1000} // ms
	for _, t := range integrationTimes {
		gain := 10 * math.Log10(t)
		fmt.Printf("%.0f ms integration: %.1f dB processing gain\n", t, gain)
	}
}