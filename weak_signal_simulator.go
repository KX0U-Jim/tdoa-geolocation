package main

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"strconv"
	"time"
)

// WeakSignalSimulator creates realistic weak/strong signal combinations for TDOA testing
type WeakSignalSimulator struct {
	SampleRate      float64
	TotalSamples    int
	BlockSamples    int
	RefFrequency    float64
	TargetFrequency float64
	Stations        map[string]Station
	Transmitter     SimulatedTransmitter
}

// SimulatedTransmitter represents a transmitter to locate
type SimulatedTransmitter struct {
	Latitude  float64 // degrees
	Longitude float64 // degrees
	Elevation float64 // meters
	RefPower  float64 // reference signal power (arbitrary units)
	TgtPower  float64 // target signal power (arbitrary units)
}

// Station represents a collector station
type Station struct {
	Name      string
	Latitude  float64
	Longitude float64
	Elevation float64
}

// Point3D represents a 3D coordinate
type Point3D struct {
	X, Y, Z float64
}

// NoiseProfile defines noise characteristics
type NoiseProfile struct {
	GaussianNoise float64 // Gaussian noise level
	ImpulseNoise  float64 // Impulse noise probability
	ImpulseLevel  float64 // Impulse noise amplitude
	PhaseDrift    float64 // Phase drift rate (rad/sec)
	DCOffset      float64 // DC offset
}

// calculateDistance3D calculates 3D distance between two lat/lon/elevation points
func calculateDistance3D(lat1, lon1, elev1, lat2, lon2, elev2 float64) float64 {
	// Convert to ECEF coordinates for accurate 3D distance
	const (
		a  = 6378137.0         // WGS84 semi-major axis (meters)
		f  = 1.0 / 298.257223563 // WGS84 flattening
		e2 = 2*f - f*f         // First eccentricity squared
	)
	
	// Convert both points to ECEF
	lat1Rad := lat1 * math.Pi / 180
	lon1Rad := lon1 * math.Pi / 180
	lat2Rad := lat2 * math.Pi / 180
	lon2Rad := lon2 * math.Pi / 180
	
	N1 := a / math.Sqrt(1 - e2*math.Sin(lat1Rad)*math.Sin(lat1Rad))
	N2 := a / math.Sqrt(1 - e2*math.Sin(lat2Rad)*math.Sin(lat2Rad))
	
	x1 := (N1 + elev1) * math.Cos(lat1Rad) * math.Cos(lon1Rad)
	y1 := (N1 + elev1) * math.Cos(lat1Rad) * math.Sin(lon1Rad)
	z1 := (N1*(1-e2) + elev1) * math.Sin(lat1Rad)
	
	x2 := (N2 + elev2) * math.Cos(lat2Rad) * math.Cos(lon2Rad)
	y2 := (N2 + elev2) * math.Cos(lat2Rad) * math.Sin(lon2Rad)
	z2 := (N2*(1-e2) + elev2) * math.Sin(lat2Rad)
	
	dx := x2 - x1
	dy := y2 - y1
	dz := z2 - z1
	
	return math.Sqrt(dx*dx + dy*dy + dz*dz)
}

// generateWeakSignal creates a realistic weak signal with various impairments
func generateWeakSignal(samples int, frequency, sampleRate, amplitude, phase float64, noiseProfile NoiseProfile) []complex64 {
	signal := make([]complex64, samples)
	
	phaseDrift := 0.0
	
	for i := 0; i < samples; i++ {
		t := float64(i) / sampleRate
		omega := 2 * math.Pi * frequency
		
		// Phase drift over time (simulates oscillator instability)
		phaseDrift += noiseProfile.PhaseDrift / sampleRate
		
		// Base signal with phase drift
		currentPhase := omega*t + phase + phaseDrift
		realPart := amplitude * math.Cos(currentPhase)
		imagPart := amplitude * math.Sin(currentPhase)
		
		// Add DC offset
		realPart += noiseProfile.DCOffset
		imagPart += noiseProfile.DCOffset
		
		// Add Gaussian noise
		if noiseProfile.GaussianNoise > 0 {
			realPart += noiseProfile.GaussianNoise * rand.NormFloat64()
			imagPart += noiseProfile.GaussianNoise * rand.NormFloat64()
		}
		
		// Add impulse noise (occasional spikes)
		if rand.Float64() < noiseProfile.ImpulseNoise {
			realPart += noiseProfile.ImpulseLevel * (2*rand.Float64() - 1)
			imagPart += noiseProfile.ImpulseLevel * (2*rand.Float64() - 1)
		}
		
		signal[i] = complex(float32(realPart), float32(imagPart))
	}
	
	return signal
}

// generateStrongSignal creates a strong, clean signal
func generateStrongSignal(samples int, frequency, sampleRate, amplitude, phase float64) []complex64 {
	signal := make([]complex64, samples)
	
	for i := 0; i < samples; i++ {
		t := float64(i) / sampleRate
		omega := 2 * math.Pi * frequency
		
		// Clean sinusoid
		realPart := amplitude * math.Cos(omega*t + phase)
		imagPart := amplitude * math.Sin(omega*t + phase)
		
		// Add minimal noise for realism (very low level)
		realPart += 0.001 * rand.NormFloat64()
		imagPart += 0.001 * rand.NormFloat64()
		
		signal[i] = complex(float32(realPart), float32(imagPart))
	}
	
	return signal
}

// simulateWeakSignalStation creates a simulated .dat file with weak reference, strong target
func simulateWeakSignalStation(sim *WeakSignalSimulator, stationName string, station Station) error {
	fmt.Printf("Simulating weak signal station: %s\n", stationName)
	
	// Calculate distance from station to transmitter
	distance := calculateDistance3D(
		station.Latitude, station.Longitude, station.Elevation,
		sim.Transmitter.Latitude, sim.Transmitter.Longitude, sim.Transmitter.Elevation)
	
	fmt.Printf("  Distance to transmitter: %.2f km\n", distance/1000)
	
	// Calculate signal travel time (time delay)
	const c = 299792458.0 // Speed of light (m/s)
	travelTime := distance / c
	refPhaseDelay := 2 * math.Pi * sim.RefFrequency * travelTime
	tgtPhaseDelay := 2 * math.Pi * sim.TargetFrequency * travelTime
	
	fmt.Printf("  Travel time: %.6f μs\n", travelTime*1e6)
	fmt.Printf("  REF phase delay: %.3f radians\n", refPhaseDelay)
	fmt.Printf("  TGT phase delay: %.3f radians\n", tgtPhaseDelay)
	
	// Calculate signal amplitudes (1/r falloff + power scaling)
	refAmplitude := sim.Transmitter.RefPower / distance * 0.1
	tgtAmplitude := sim.Transmitter.TgtPower / distance * 0.1
	
	fmt.Printf("  REF signal amplitude: %.9f (weak)\n", refAmplitude)
	fmt.Printf("  TGT signal amplitude: %.9f (strong)\n", tgtAmplitude)
	
	// Define noise profiles
	// Weak reference signal - lots of impairments
	weakNoiseProfile := NoiseProfile{
		GaussianNoise: refAmplitude * 0.8,  // 80% of signal level - very noisy
		ImpulseNoise:  0.001,               // 0.1% impulse noise probability
		ImpulseLevel:  refAmplitude * 5.0,  // Strong impulses
		PhaseDrift:    0.05,                // Phase drift: 0.05 rad/sec
		DCOffset:      refAmplitude * 0.1,  // DC offset
	}
	
	// Strong target signal - minimal noise
	strongNoiseProfile := NoiseProfile{
		GaussianNoise: tgtAmplitude * 0.02, // 2% noise level - very clean
		ImpulseNoise:  0.0001,              // Very rare impulses
		ImpulseLevel:  tgtAmplitude * 0.5,  // Weaker impulses
		PhaseDrift:    0.001,               // Minimal phase drift
		DCOffset:      tgtAmplitude * 0.01, // Small DC offset
	}
	
	// Generate the three blocks: ref, target, ref
	// Block 1: Weak reference signal with lots of noise
	block1 := generateWeakSignal(sim.BlockSamples, sim.RefFrequency, 
		sim.SampleRate, refAmplitude, refPhaseDelay, weakNoiseProfile)
	
	// Block 2: Strong target signal with minimal noise
	block2 := generateStrongSignal(sim.BlockSamples, sim.TargetFrequency, 
		sim.SampleRate, tgtAmplitude, tgtPhaseDelay)
	
	// Block 3: Weak reference signal again (same as block 1 characteristics)
	block3 := generateWeakSignal(sim.BlockSamples, sim.RefFrequency, 
		sim.SampleRate, refAmplitude, refPhaseDelay, weakNoiseProfile)
	
	// Combine blocks into single data stream
	totalData := make([]complex64, sim.TotalSamples)
	copy(totalData[0:sim.BlockSamples], block1)
	copy(totalData[sim.BlockSamples:2*sim.BlockSamples], block2)  
	copy(totalData[2*sim.BlockSamples:3*sim.BlockSamples], block3)
	
	// Convert to raw I/Q bytes (RTL-SDR format: unsigned 8-bit I, Q)
	rawData := make([]byte, len(totalData)*2)
	for i, sample := range totalData {
		// Convert from float32 to [0,255] uint8 centered at 127.5
		iVal := real(sample)*127.5 + 127.5
		qVal := imag(sample)*127.5 + 127.5
		
		// Clamp to valid range
		if iVal < 0 { iVal = 0 }
		if iVal > 255 { iVal = 255 }
		if qVal < 0 { qVal = 0 }
		if qVal > 255 { qVal = 255 }
		
		rawData[i*2] = byte(iVal)
		rawData[i*2+1] = byte(qVal)
	}
	
	// Write to file with same naming convention as collector
	timestamp := int64(1754950000) // Different timestamp from perfect simulator
	filename := fmt.Sprintf("weak-%s-%d.dat", stationName, timestamp)
	
	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("failed to create %s: %v", filename, err)
	}
	defer file.Close()
	
	_, err = file.Write(rawData)
	if err != nil {
		return fmt.Errorf("failed to write %s: %v", filename, err)
	}
	
	fmt.Printf("  Generated: %s (%.1f MB)\n", filename, float64(len(rawData))/1e6)
	
	// Calculate SNR estimates for verification
	refSNR := 20 * math.Log10(refAmplitude / weakNoiseProfile.GaussianNoise)
	tgtSNR := 20 * math.Log10(tgtAmplitude / strongNoiseProfile.GaussianNoise)
	fmt.Printf("  Estimated REF SNR: %.1f dB (weak)\n", refSNR)
	fmt.Printf("  Estimated TGT SNR: %.1f dB (strong)\n", tgtSNR)
	
	return nil
}

// loadStationsFromCSV loads station data from CSV file (simplified version)
func loadStationsFromCSV(csvPath string) (map[string]Station, error) {
	// For now, hard-code the stations from lat-lon-table.csv
	stations := make(map[string]Station)
	
	stations["162400000"] = Station{
		Name: "162400000",
		Latitude: 41.25703803095629,
		Longitude: -95.95512763589404,
		Elevation: 349.07,
	}
	stations["kx0u"] = Station{
		Name: "kx0u", 
		Latitude: 41.18660274289527,
		Longitude: -95.96064116595667,
		Elevation: 355.69,
	}
	stations["n3pay"] = Station{
		Name: "n3pay",
		Latitude: 41.24669616513154,
		Longitude: -96.08366304481238,
		Elevation: 329.0,
	}
	stations["kf0mtl"] = Station{
		Name: "kf0mtl",
		Latitude: 41.32916620016985,
		Longitude: -96.03513381562004, 
		Elevation: 373.18,
	}
	
	return stations, nil
}

func main() {
	// Initialize random number generator
	rand.Seed(time.Now().UnixNano())
	
	if len(os.Args) < 7 {
		fmt.Printf("Usage: %s <csv_file> <target_freq> <tx_lat> <tx_lon> <tx_elev> <ref_power> <tgt_power>\n", os.Args[0])
		fmt.Printf("Example: ./weak_signal_simulator lat-lon-table.csv 92300000 41.20 -96.00 400 10 1000\n")
		fmt.Printf("This creates weak reference signals with noise and strong target signals\n")
		fmt.Printf("ref_power: Reference signal power (low = weak signal)\n")
		fmt.Printf("tgt_power: Target signal power (high = strong signal)\n")
		os.Exit(1)
	}

	csvFile := os.Args[1]
	targetFreq, _ := strconv.ParseFloat(os.Args[2], 64)
	txLat, _ := strconv.ParseFloat(os.Args[3], 64)
	txLon, _ := strconv.ParseFloat(os.Args[4], 64) 
	txElev, _ := strconv.ParseFloat(os.Args[5], 64)
	refPower, _ := strconv.ParseFloat(os.Args[6], 64)
	tgtPower, _ := strconv.ParseFloat(os.Args[7], 64)

	fmt.Println("=== Weak Signal TDOA Simulator ===")
	fmt.Printf("Target frequency: %.3f MHz\n", targetFreq/1e6)
	fmt.Printf("Transmitter: %.6f°, %.6f°, %.1fm\n", txLat, txLon, txElev)
	fmt.Printf("Reference power: %.1f (weak signal)\n", refPower)
	fmt.Printf("Target power: %.1f (strong signal)\n", tgtPower)

	// Load stations
	stations, err := loadStationsFromCSV(csvFile)
	if err != nil {
		log.Fatalf("Failed to load stations: %v", err)
	}

	fmt.Printf("Loaded %d stations\n", len(stations))

	// Setup simulation config for realistic weak signal testing
	sim := &WeakSignalSimulator{
		SampleRate:   2000000.0, // 2 MHz
		TotalSamples: 60000000,  // 30 seconds total (3×10 second blocks)
		BlockSamples: 20000000,  // 10 seconds per block (increased for contamination dilution)
		RefFrequency: 162400000.0, // NOAA weather
		TargetFrequency: targetFreq,
		Stations: stations,
		Transmitter: SimulatedTransmitter{
			Latitude:  txLat,
			Longitude: txLon, 
			Elevation: txElev,
			RefPower:  refPower,  // Weak reference
			TgtPower:  tgtPower,  // Strong target
		},
	}

	// Generate simulated data for each collector station
	collectorStations := []string{"kx0u", "n3pay", "kf0mtl"}
	
	for _, stationName := range collectorStations {
		station, exists := stations[stationName]
		if !exists {
			fmt.Printf("Warning: Station %s not found in CSV\n", stationName)
			continue
		}
		
		err := simulateWeakSignalStation(sim, stationName, station)
		if err != nil {
			log.Printf("Failed to simulate %s: %v", stationName, err)
		}
	}

	fmt.Printf("\n=== Weak Signal Simulation Complete ===\n")
	fmt.Printf("Generated files: weak-kx0u-*.dat, weak-n3pay-*.dat, weak-kf0mtl-*.dat\n")
	fmt.Printf("Test with: ./processor 162400000 %.0f lat-lon-table.csv weak-*.dat\n", targetFreq)
	fmt.Printf("Expected location: %.6f°, %.6f°, %.1fm\n", txLat, txLon, txElev)
	fmt.Printf("\nCharacteristics:\n")
	fmt.Printf("- Reference signal: WEAK (%.0f power) with noise, phase drift, impulses\n", refPower)
	fmt.Printf("- Target signal: STRONG (%.0f power) with minimal noise\n", tgtPower)
	fmt.Printf("- Tests processor performance with realistic weak reference signals\n")
}