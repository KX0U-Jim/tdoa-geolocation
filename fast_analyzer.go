package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"sort"
)

const (
	sampleRate     = 2000000.0
	bytesPerSample = 2
	expectedDCBias = 127.5
)

type FastAnalysis struct {
	TotalSamples int
	IAvg, QAvg   float64
	IStdDev, QStdDev float64
	SNREstimate  float64
	PowerLevel   float64
	HasClipping  bool
	HasOverload  bool
}

func main() {
	if len(os.Args) < 2 {
		fmt.Printf("Usage: %s <data_file.dat>\n", os.Args[0])
		fmt.Printf("Fast signal quality analyzer for gain sweeps\n")
		os.Exit(1)
	}

	filename := os.Args[1]

	// Fast analysis
	refAnalysis, targetAnalysis, err := fastAnalyzeDualFrequencyFile(filename)
	if err != nil {
		fmt.Printf("Error: %v\n", err)
		os.Exit(1)
	}

	// Compact output for CSV processing
	fmt.Printf("REF,%.1f,%.1f,%t,%t\n", 
		refAnalysis.SNREstimate, refAnalysis.PowerLevel,
		refAnalysis.HasClipping, refAnalysis.HasOverload)
	
	fmt.Printf("TGT,%.1f,%.1f,%t,%t\n",
		targetAnalysis.SNREstimate, targetAnalysis.PowerLevel,
		targetAnalysis.HasClipping, targetAnalysis.HasOverload)
}

func fastAnalyzeDualFrequencyFile(filename string) (*FastAnalysis, *FastAnalysis, error) {
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

	// For speed: only analyze first 65536 samples of each signal
	maxSamples := 65536
	samplesPerBlock := totalSamples / 3
	
	if samplesPerBlock*2 < maxSamples {
		maxSamples = samplesPerBlock * 2
	}
	if samplesPerBlock < maxSamples/2 {
		maxSamples = samplesPerBlock
	}

	// Read just what we need
	samples := make([]byte, maxSamples*3*bytesPerSample)
	_, err = file.Read(samples)
	if err != nil {
		return nil, nil, fmt.Errorf("failed to read samples: %v", err)
	}

	// Extract signals for fast analysis
	refSamples := make([]byte, maxSamples*bytesPerSample*2)
	copy(refSamples[0:maxSamples*bytesPerSample], samples[0:maxSamples*bytesPerSample])
	copy(refSamples[maxSamples*bytesPerSample:], samples[maxSamples*2*bytesPerSample:maxSamples*3*bytesPerSample])

	targetSamples := samples[maxSamples*bytesPerSample : maxSamples*2*bytesPerSample]

	// Fast analysis
	refAnalysis := fastAnalyzeSamples(refSamples, maxSamples*2)
	targetAnalysis := fastAnalyzeSamples(targetSamples, maxSamples)

	return refAnalysis, targetAnalysis, nil
}

func fastAnalyzeSamples(samples []byte, totalSamples int) *FastAnalysis {
	analysis := &FastAnalysis{TotalSamples: totalSamples}

	// Basic statistics - optimized loop
	var iSum, qSum, iSumSq, qSumSq float64
	var iMin, iMax, qMin, qMax byte = 255, 0, 255, 0

	for i := 0; i < totalSamples; i++ {
		iVal := samples[i*2]
		qVal := samples[i*2+1]
		
		iFloat := float64(iVal)
		qFloat := float64(qVal)
		
		iSum += iFloat
		qSum += qFloat
		iSumSq += iFloat * iFloat
		qSumSq += qFloat * qFloat

		if iVal < iMin { iMin = iVal }
		if iVal > iMax { iMax = iVal }
		if qVal < qMin { qMin = qVal }
		if qVal > qMax { qMax = qVal }
	}

	// Calculate stats
	n := float64(totalSamples)
	analysis.IAvg = iSum / n
	analysis.QAvg = qSum / n
	analysis.IStdDev = math.Sqrt((iSumSq/n) - (analysis.IAvg * analysis.IAvg))
	analysis.QStdDev = math.Sqrt((qSumSq/n) - (analysis.QAvg * analysis.QAvg))

	// Power level
	analysis.PowerLevel = 20 * math.Log10(math.Sqrt(analysis.IStdDev*analysis.IStdDev + analysis.QStdDev*analysis.QStdDev))

	// Quality flags
	analysis.HasClipping = (iMin == 0 || iMax == 255 || qMin == 0 || qMax == 255)
	analysis.HasOverload = (analysis.IStdDev < 2 || analysis.QStdDev < 2)

	// Fast SNR calculation using smaller FFT
	analysis.SNREstimate = fastSNRCalculation(samples, totalSamples)

	return analysis
}

func fastSNRCalculation(samples []byte, totalSamples int) float64 {
	// Use only 8192 samples for very fast FFT
	analysisSize := 8192
	if totalSamples < analysisSize {
		analysisSize = totalSamples
	}

	// Take samples from middle
	startIdx := (totalSamples - analysisSize) / 2
	analysisBytes := samples[startIdx*2 : (startIdx+analysisSize)*2]

	// Fast conversion to complex
	complexSamples := make([]complex128, analysisSize)
	for i := 0; i < analysisSize; i++ {
		iVal := (float64(analysisBytes[i*2]) - 127.5) / 127.5
		qVal := (float64(analysisBytes[i*2+1]) - 127.5) / 127.5
		complexSamples[i] = complex(iVal, qVal)
	}

	// Simple windowing (Hanning)
	for i, sample := range complexSamples {
		window := 0.5 - 0.5*math.Cos(2*math.Pi*float64(i)/float64(analysisSize-1))
		complexSamples[i] = complex(window, 0) * sample
	}

	// Fast DFT (smaller size)
	fft := fastDFT(complexSamples)

	// Power spectral density
	psd := make([]float64, len(fft))
	for i, c := range fft {
		psd[i] = cmplx.Abs(c) * cmplx.Abs(c)
	}

	// Quick signal/noise separation
	sortedPSD := make([]float64, len(psd))
	copy(sortedPSD, psd)
	sort.Float64s(sortedPSD)

	// Signal: top 10%, Noise: bottom 40%
	signalThreshold := sortedPSD[int(0.9*float64(len(sortedPSD)))]
	noiseThreshold := sortedPSD[int(0.4*float64(len(sortedPSD)))]

	var signalPower, noisePower float64
	signalCount, noiseCount := 0, 0

	for _, power := range psd {
		if power >= signalThreshold {
			signalPower += power
			signalCount++
		} else if power <= noiseThreshold {
			noisePower += power
			noiseCount++
		}
	}

	if signalCount > 0 { signalPower /= float64(signalCount) }
	if noiseCount > 0 { noisePower /= float64(noiseCount) }

	if noisePower > 0 && signalPower > noisePower {
		return 10 * math.Log10(signalPower/noisePower)
	}

	return -20.0
}

func fastDFT(samples []complex128) []complex128 {
	// Optimized DFT for small sizes
	n := len(samples)
	if n > 8192 { n = 8192 } // Limit for speed
	
	dft := make([]complex128, n)
	
	// Pre-calculate twiddle factors
	twiddle := make([]complex128, n)
	for k := 0; k < n; k++ {
		angle := -2 * math.Pi * float64(k) / float64(n)
		twiddle[k] = complex(math.Cos(angle), math.Sin(angle))
	}
	
	for k := 0; k < n; k++ {
		sum := complex(0, 0)
		for i := 0; i < n; i++ {
			idx := (k * i) % n
			sum += samples[i] * twiddle[idx]
		}
		dft[k] = sum
	}
	
	return dft
}