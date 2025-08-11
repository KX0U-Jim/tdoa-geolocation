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

// loadIQData loads I/Q data from .dat file (placeholder for now)
func (p *TDOAProcessor) loadIQData(filename string) ([]complex64, error) {
	// TODO: Implement actual I/Q data loading
	// For now, return empty slice to allow compilation
	fmt.Printf("Loading I/Q data from: %s\n", filename)
	return make([]complex64, 0), nil
}

// extractReferenceSignal extracts reference frequency samples from dual-freq data
func (p *TDOAProcessor) extractReferenceSignal(data []complex64) []complex64 {
	// TODO: Implement signal extraction from dual-frequency pattern
	// librtlsdr-2freq pattern: freq1, freq2, freq1 (blocks of samples)
	// Extract blocks 1 and 3 (reference frequency)
	fmt.Println("Extracting reference signal from dual-frequency data")
	return make([]complex64, 0)
}

// crossCorrelate performs cross-correlation between two signals
func (p *TDOAProcessor) crossCorrelate(signal1, signal2 []complex64) (int, float64) {
	// TODO: Implement cross-correlation algorithm
	// Returns: time delay (samples), correlation peak value
	fmt.Println("Performing cross-correlation")
	return 0, 0.0
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

		collectorData = append(collectorData, struct {
			Station Station
			Data    []complex64
			RefSig  []complex64
		}{
			Station: station,
			Data:    data,
			RefSig:  refSig,
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

	// Perform cross-correlation between reference signals
	fmt.Println("\nCross-correlation analysis:")
	var timeDifferences []float64
	
	for i := 0; i < len(collectorData); i++ {
		for j := i + 1; j < len(collectorData); j++ {
			delay, correlation := p.crossCorrelate(collectorData[i].RefSig, collectorData[j].RefSig)
			
			// Convert sample delay to time difference
			sampleRate := 2e6 // 2 Msps
			timeDiff := float64(delay) / sampleRate
			
			timeDifferences = append(timeDifferences, timeDiff)
			
			fmt.Printf("%s - %s: delay=%d samples (%.3f μs), correlation=%.3f\n", 
				collectorData[i].Station.Name, collectorData[j].Station.Name, 
				delay, timeDiff*1e6, correlation)
		}
	}

	// TODO: Implement TDOA triangulation using time differences
	fmt.Println("\nTDOA triangulation:")
	fmt.Println("TODO: Calculate target position from time differences")

	return nil
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