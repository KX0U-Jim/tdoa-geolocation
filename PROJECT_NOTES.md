# TDOA Project Key Information

## Signal Information

### 92.3 MHz - FM Broadcast Station (Reference Signal)
- **Type**: Strong commercial WBFM broadcaster 
- **Reception**: All collectors receive this signal with great strength
- **Purpose**: Reference signal for timing synchronization between collectors
- **Expected**: Should provide excellent correlation for TDOA timing measurements

### 162.4 MHz - NOAA Weather Radio (Target Signal) 
- **Type**: NOAA weather radio NBFM broadcast
- **Reception**: All collectors receive this well enough for close to full quieting or better
- **Purpose**: Target signal to locate via TDOA triangulation
- **Expected**: Should provide sufficient signal strength for correlation

## Hardware Information

### RTL-SDR Receivers
- **Oscillator**: All units equipped with TCXO (Temperature Compensated Crystal Oscillator)
- **Frequency Stability**: Excellent - no concerns about frequency drift between units
- **Implication**: Timing differences should be due to signal propagation delays only, not hardware drift

## Baseline Distances
- kx0u - n3pay: 12.29 km
- kx0u - kf0mtl: 17.02 km  
- n3pay - kf0mtl: 10.02 km

## Expected TDOA Delays
- **Maximum possible delay**: ~57 μs (17 km ÷ speed of light)
- **Realistic range**: 0-57 μs between any collector pair
- **Any delays >100 μs indicate correlation algorithm problems**

## Key Technical Notes
- Both signals should provide strong, reliable correlation
- TCXO ensures frequency stability - no need for frequency drift compensation
- Focus on correlation algorithm accuracy, not hardware timing issues