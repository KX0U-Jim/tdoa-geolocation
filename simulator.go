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

// SimulatedTransmitter represents a transmitter to locate
type SimulatedTransmitter struct {
	Latitude  float64 // degrees
	Longitude float64 // degrees
	Elevation float64 // meters
	Power     float64 // signal power (arbitrary units)
}

// SimulationConfig contains simulation parameters
type SimulationConfig struct {
	SampleRate      float64 // 2 MHz
	TotalSamples    int     // 60M samples total (30 seconds)
	BlockSamples    int     // 20M samples per block
	RefFrequency    float64 // 162.4 MHz (reference signal)
	TargetFrequency float64 // Variable (signal to locate)
	NoiseLevel      float64 // Noise floor level
	Stations        map[string]Station // From lat-lon-table.csv
	Transmitter     SimulatedTransmitter
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

// generatePerfectSignal creates a perfect reference signal 
func generatePerfectSignal(samples int, frequency, sampleRate, amplitude, phase float64) []complex64 {
	signal := make([]complex64, samples)
	
	for i := 0; i < samples; i++ {
		t := float64(i) / sampleRate
		omega := 2 * math.Pi * frequency
		
		// Pure sinusoid with specified phase
		realPart := float32(amplitude * math.Cos(omega*t + phase))
		imagPart := float32(amplitude * math.Sin(omega*t + phase))
		
		signal[i] = complex(realPart, imagPart)
	}
	
	return signal
}

// addNoise adds Gaussian noise to signal
func addNoise(signal []complex64, noiseLevel float64) []complex64 {
	noisy := make([]complex64, len(signal))
	
	for i, sample := range signal {
		// Add Gaussian noise to both I and Q
		noiseI := float32(noiseLevel * (2*rand.Float64() - 1))
		noiseQ := float32(noiseLevel * (2*rand.Float64() - 1))
		
		noisy[i] = complex(real(sample) + noiseI, imag(sample) + noiseQ)
	}
	
	return noisy
}

// simulateStation creates a simulated .dat file for one collector station
func simulateStation(config *SimulationConfig, stationName string, station Station) error {
	fmt.Printf("Simulating station: %s\n", stationName)
	
	// Calculate distance from station to transmitter
	distance := calculateDistance3D(
		station.Latitude, station.Longitude, station.Elevation,
		config.Transmitter.Latitude, config.Transmitter.Longitude, config.Transmitter.Elevation)
	
	fmt.Printf("  Distance to transmitter: %.2f km\n", distance/1000)
	
	// Calculate signal travel time (time delay)
	c := 299792458.0 // Speed of light (m/s)
	travelTime := distance / c
	phaseDelay := 2 * math.Pi * config.TargetFrequency * travelTime
	
	fmt.Printf("  Travel time: %.6f μs\n", travelTime*1e6)
	fmt.Printf("  Phase delay: %.3f radians\n", phaseDelay)
	
	// Calculate signal amplitude (1/r falloff + some variations)
	amplitude := config.Transmitter.Power / distance
	amplitude *= 0.1 // Scale factor to reasonable levels
	
	fmt.Printf("  Signal amplitude: %.9f\n", amplitude)
	
	// Generate the three blocks: ref, target, ref
	// Block 1: Reference signal (perfect, no delay)
	block1 := generatePerfectSignal(config.BlockSamples, config.RefFrequency, 
		config.SampleRate, 0.01, 0.0) // Small amplitude for reference
	block1 = addNoise(block1, config.NoiseLevel)
	
	// Block 2: Target signal (with distance-based delay and amplitude)
	block2 := generatePerfectSignal(config.BlockSamples, config.TargetFrequency, 
		config.SampleRate, amplitude, phaseDelay)
	block2 = addNoise(block2, config.NoiseLevel)
	
	// Block 3: Reference signal again (same as block 1)
	block3 := generatePerfectSignal(config.BlockSamples, config.RefFrequency, 
		config.SampleRate, 0.01, 0.0)
	block3 = addNoise(block3, config.NoiseLevel)
	
	// Combine blocks into single data stream
	totalData := make([]complex64, config.TotalSamples)
	copy(totalData[0:config.BlockSamples], block1)
	copy(totalData[config.BlockSamples:2*config.BlockSamples], block2)  
	copy(totalData[2*config.BlockSamples:3*config.BlockSamples], block3)
	
	// Convert to raw I/Q bytes (RTL-SDR format: unsigned 8-bit I, Q)
	rawData := make([]byte, len(totalData)*2)
	for i, sample := range totalData {
		// Convert from [-1,1] float32 to [0,255] uint8 centered at 127.5
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
	timestamp := int64(1754900000) // Simulated timestamp
	filename := fmt.Sprintf("sim-%s-%d.dat", stationName, timestamp)
	
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
	return nil
}

// Station represents a collector station (copy from processor.go to avoid import conflict)
type Station struct {
	Name      string
	Latitude  float64
	Longitude float64
	Elevation float64
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
		fmt.Printf("Usage: %s <csv_file> <target_freq> <tx_lat> <tx_lon> <tx_elev> <tx_power>\n", os.Args[0])
		fmt.Println("Example: ./simulator lat-lon-table.csv 101700000 41.20 -96.00 400 1000")
		fmt.Println("This creates perfect simulated .dat files for TDOA testing")
		os.Exit(1)
	}

	csvFile := os.Args[1]
	targetFreq, _ := strconv.ParseFloat(os.Args[2], 64)
	txLat, _ := strconv.ParseFloat(os.Args[3], 64)
	txLon, _ := strconv.ParseFloat(os.Args[4], 64) 
	txElev, _ := strconv.ParseFloat(os.Args[5], 64)
	txPower, _ := strconv.ParseFloat(os.Args[6], 64)

	fmt.Println("=== TDOA Signal Simulator ===")
	fmt.Printf("Target frequency: %.3f MHz\n", targetFreq/1e6)
	fmt.Printf("Transmitter: %.6f°, %.6f°, %.1fm\n", txLat, txLon, txElev)
	fmt.Printf("Transmitter power: %.1f\n", txPower)

	// Load stations
	stations, err := loadStationsFromCSV(csvFile)
	if err != nil {
		log.Fatalf("Failed to load stations: %v", err)
	}

	fmt.Printf("Loaded %d stations\n", len(stations))

	// Setup simulation config - very small for fast testing
	config := &SimulationConfig{
		SampleRate:   2000000.0, // 2 MHz
		TotalSamples: 60000000,  // 30 seconds total (3×10 second blocks)
		BlockSamples: 20000000,  // 10 seconds per block (increased for contamination dilution)
		RefFrequency: 162400000.0, // NOAA weather
		TargetFrequency: targetFreq,
		NoiseLevel: 0.01, // Low noise for perfect signals
		Stations: stations,
		Transmitter: SimulatedTransmitter{
			Latitude:  txLat,
			Longitude: txLon, 
			Elevation: txElev,
			Power:     txPower,
		},
	}

	// Generate simulated data for each collector station (not reference)
	collectorStations := []string{"kx0u", "n3pay", "kf0mtl"}
	
	for _, stationName := range collectorStations {
		station, exists := stations[stationName]
		if !exists {
			fmt.Printf("Warning: Station %s not found in CSV\n", stationName)
			continue
		}
		
		err := simulateStation(config, stationName, station)
		if err != nil {
			log.Printf("Failed to simulate %s: %v", stationName, err)
		}
	}

	fmt.Printf("\n=== Simulation Complete ===\n")
	fmt.Printf("Generated files: sim-kx0u-*.dat, sim-n3pay-*.dat, sim-kf0mtl-*.dat\n")
	fmt.Printf("Test with: ./processor 162400000 %.0f lat-lon-table.csv sim-*.dat\n", targetFreq)
	fmt.Printf("Expected location: %.6f°, %.6f°, %.1fm\n", txLat, txLon, txElev)
}