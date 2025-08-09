# Usage Guide

## Basic Operations

### Planning a Collection

1. **Choose frequencies:**
   - Reference: Strong FM broadcast station (e.g., 96.9 MHz)
   - Target: Signal to locate (e.g., 162.55 MHz NOAA weather radio)

2. **Coordinate timing:**
   - All stations must start simultaneously
   - Use epoch seconds for precision timing
   - Allow 30+ seconds for preparation

3. **Select duration:**
   - Start with 10-30 seconds for testing
   - Longer collections give better correlation but risk drift
   - Maximum 100 seconds due to memory constraints

### Running a Collection

#### Calculate Start Time
```bash
# Start 60 seconds from now
start_time=$(($(date +%s) + 60))
echo "Start time: $start_time ($(date -d @$start_time))"
```

#### Execute on All Stations

**Station kx0u:**
```bash
./collector --duration=30 96900000 162550000 $start_time kx0u
```

**Station n3pay:**
```bash
./collector --duration=30 96900000 162550000 $start_time n3pay
```

**Station kf0mtl:**
```bash
./collector --duration=30 96900000 162550000 $start_time kf0mtl
```

All commands must use:
- **Same frequencies** (reference and target)
- **Same start time** (epoch seconds)
- **Same duration**
- **Unique station IDs**

### Data Validation

After collection, validate each station's data:

```bash
# Check each station's output
./reader kx0u-$start_time.dat 30
./reader n3pay-$start_time.dat 30  
./reader kf0mtl-$start_time.dat 30
```

**Expected Results:**
```
✓ File size matches expected exactly
✓ Sample count follows 3×n pattern: 20000000 samples per frequency block
✓ I/Q components show good dynamic range
✓ DC bias appears normal (I: 127.2, Q: 128.1)
✓ Data contains non-zero values
```

## Advanced Usage

### Custom Sample Rates

The collector defaults to 2 Msps, but you can modify for specific needs:

```go
// In collector.go, change this constant:
const sampleRate = 2400000 // 2.4 Msps for wider bandwidth
```

Remember to rebuild: `go build -o collector collector.go`

### Extended Duration Collections

For longer collections, monitor system resources:

```bash
# Monitor during collection
htop

# Check for USB overruns in dmesg
dmesg | tail -20
```

### Frequency Selection Guidelines

**Reference Signal Selection:**
- Choose strongest local FM station
- Verify coverage at all collector locations
- Avoid stations with frequent dead air or automation

**Target Signal Selection:**
- Ensure signal is receivable at all stations
- Stationary targets work best initially
- Test with known locations first

## File Management

### Naming Convention

Files are automatically named: `{station_id}-{start_epoch}.dat`

Examples:
```
kx0u-1692123456.dat    # Station kx0u, epoch time 1692123456
n3pay-1692123456.dat   # Station n3pay, same collection
kf0mtl-1692123456.dat  # Station kf0mtl, same collection
```

### Storage Requirements

**Per collection:**
- Duration × 2 Msps × 2 bytes = File size
- 30 seconds = 120 MB per station
- 100 seconds = 400 MB per station

**Disk space planning:**
```bash
# Check available space
df -h

# Clean up old collections
rm *-$(date -d '1 week ago' +%s).dat
```

### Data Transfer

**Centralize collections:**
```bash
# Copy from remote stations
scp user@n3pay-ip:~/*-$start_time.dat ./
scp user@kf0mtl-ip:~/*-$start_time.dat ./

# Organize by collection time
mkdir collection-$start_time
mv *-$start_time.dat collection-$start_time/
```

## Troubleshooting

### Common Issues

**"Error executing rtl_sdr"**
- Check RTL-SDR connection: `lsusb | grep Realtek`
- Verify permissions: `groups $USER` should include `plugdev`
- Test manually: `rtl_test`

**"Sample count doesn't follow 3×n pattern"**
- USB overrun occurred during collection
- Reduce duration or check system load
- Ensure no other processes using RTL-SDR

**"DC bias may be off"**
- RTL-SDR hardware issue or poor connection
- Try different USB port
- Check antenna connection

**"Data appears to be all zeros"**
- RTL-SDR not receiving signal
- Check antenna connection and tuning
- Verify frequency is correct

### Performance Optimization

**Reduce sample loss:**
```bash
# Stop unnecessary services during collection
sudo systemctl stop bluetooth
sudo systemctl stop cups

# Increase USB buffer size
echo 'SUBSYSTEM=="usb", ATTRS{idVendor}=="0bda", ATTRS{idProduct}=="2838", ATTR{bConfigurationValue}="1"' | sudo tee -a /etc/udev/rules.d/20-rtlsdr.rules
```

**Monitor system performance:**
```bash
# Before collection
free -h
iostat 1 5

# During collection (different terminal)
top -p $(pgrep collector)
```

## Operational Procedures

### Pre-Collection Checklist

1. ✅ All stations have correct time synchronization
2. ✅ RTL-SDR devices connected and tested  
3. ✅ Reference and target frequencies verified
4. ✅ Start time calculated and communicated
5. ✅ Sufficient disk space available
6. ✅ No interfering processes running

### Collection Workflow

1. **Coordinate timing** across all team members
2. **Execute simultaneously** on all stations
3. **Monitor completion** - all should finish together
4. **Validate immediately** using reader
5. **Transfer data** to central location if needed
6. **Archive successfully** validated collections

### Post-Collection Analysis

1. **Data validation** on all stations
2. **Quality assessment** - SNR, timing accuracy
3. **Pattern matching** (when processor is ready)
4. **Results comparison** between stations
5. **Documentation** of collection parameters

## Automation

### Scripted Collections

Create collection scripts for repeated operations:

```bash
#!/bin/bash
# collect.sh - Automated collection script

REFERENCE_FREQ=96900000
TARGET_FREQ=162550000
DURATION=30
STATION_ID=$(hostname)

# Calculate start time
START_TIME=$(($(date +%s) + 60))

echo "Collection starting at $(date -d @$START_TIME)"
echo "Reference: $REFERENCE_FREQ Hz"
echo "Target: $TARGET_FREQ Hz"
echo "Duration: $DURATION seconds"

# Execute collection
./collector --duration=$DURATION $REFERENCE_FREQ $TARGET_FREQ $START_TIME $STATION_ID

# Validate results
./reader ${STATION_ID}-${START_TIME}.dat $DURATION
```

### Scheduled Collections

Use cron for regular collections:

```bash
# Edit crontab
crontab -e

# Add scheduled collection (example: every day at 2 AM)
0 2 * * * cd /home/user/tdoa && ./collect.sh >> collection.log 2>&1
```

## Safety and Legal

### RF Safety
- Use appropriate antennas for frequency ranges
- Follow local RF exposure guidelines
- Ensure proper grounding of equipment

### Legal Compliance
- Verify local regulations for RF monitoring
- Respect privacy and communication laws  
- Use only for amateur radio experimentation
- Coordinate with authorities if required

### Best Practices
- Document all collection parameters
- Maintain equipment logs and calibration
- Share results with amateur radio community
- Contribute improvements back to project