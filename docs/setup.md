# Setup Guide

## System Requirements

### Hardware
- **3 RTL-SDR receivers** (V4/V5 with TCXO recommended)
- **3 collector stations** (Raspberry Pi 4+ or x86_64 Linux)
- **Internet connectivity** for synchronization
- **Adequate separation** between collectors (10-50 km for good geometry)

### Software
- Linux-based OS (Ubuntu 20.04+, Raspberry Pi OS)
- Go 1.19+ (installed by deployment script)
- Build tools (gcc, cmake, make)
- USB access permissions

## Installation

### Automated Deployment

The easiest way to set up each collector station:

```bash
# Copy files to each station
scp collector.go reader.go deploy.sh user@station-ip:~/

# SSH to each station and run deployment
ssh user@station-ip
cd ~/
chmod +x deploy.sh
./deploy.sh

# Reboot for USB permissions to take effect
sudo reboot
```

### Manual Installation

If you prefer manual setup or need to troubleshoot:

#### 1. Install Dependencies

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install -y build-essential cmake git libusb-1.0-0-dev pkg-config golang-go
```

**Raspberry Pi:**
```bash
sudo apt update
sudo apt install -y build-essential cmake git libusb-1.0-0-dev pkg-config golang-go libraspberrypi-dev
```

#### 2. Build librtlsdr-2freq

```bash
git clone https://github.com/DC9ST/librtlsdr-2freq.git
cd librtlsdr-2freq
mkdir build && cd build
cmake ..
make -j$(nproc)
```

Verify the build:
```bash
./src/rtl_sdr --help
# Should show -f and -h dual-frequency parameters
```

#### 3. Build Collector Programs

```bash
cd ~/  # or your project directory
go mod init tdoa-collector
go build -o collector collector.go
go build -o reader reader.go
```

#### 4. Set Up Permissions

```bash
# Add user to plugdev group
sudo usermod -a -G plugdev $USER

# Install udev rules
sudo cp librtlsdr-2freq/rtl-sdr.rules /etc/udev/rules.d/
sudo udevadm control --reload-rules
sudo udevadm trigger

# Log out and back in (or reboot) for group changes
```

## Configuration

### Collector Station Coordinates

Record the precise GPS coordinates of each collector station:

```bash
# Example coordinates (replace with your actual locations)
kx0u:   41.18660274289527, -95.96064116595667
n3pay:  41.24669616513154, -96.08366304481238  
kf0mtl: 41.32916620016985, -96.03513381562004
```

### RTL-SDR Testing

Verify each RTL-SDR is working:

```bash
# Test basic RTL-SDR functionality
rtl_test -t

# Test dual-frequency capability (will fail without antenna, but should start)
timeout 5 ./librtlsdr-2freq/build/src/rtl_sdr -f 96900000 -h 162550000 -n 1000 test.dat
```

### Network Time Synchronization

Ensure all stations have synchronized clocks:

```bash
# Install NTP
sudo apt install ntp

# Check time synchronization
timedatectl status

# Should show "NTP service: active" and "synchronized: yes"
```

## Station-Specific Configuration

### Update Collector Paths

If you moved files or have different directory structure:

```go
// In collector.go, update this line if needed:
rtlSdrPath := "/path/to/your/librtlsdr-2freq/build/src/rtl_sdr"
```

### Performance Tuning

**For Raspberry Pi:**
```bash
# Increase GPU memory split for better performance
echo "gpu_mem=128" | sudo tee -a /boot/config.txt

# Use fast SD card (Class 10 or better)
# Consider USB 3.0 SSD for intensive operations
```

**For all stations:**
```bash
# Increase USB buffer sizes if experiencing sample loss
echo 'SUBSYSTEM=="usb", ATTRS{idVendor}=="0bda", ATTRS{idProduct}=="2838", MODE="0666"' | sudo tee /etc/udev/rules.d/20-rtlsdr.rules
```

## Validation

### Test Collection

Run a short test collection:

```bash
# Calculate start time 30 seconds in future
start_time=$(($(date +%s) + 30))

# Run 5-second test collection
./collector --duration=5 96900000 162550000 $start_time test-$(hostname)

# Validate the collected data
./reader test-$(hostname)-$start_time.dat 5
```

Expected output:
```
✓ File size matches expected exactly
✓ File size is valid multiple of sample size  
✓ Sample count follows 3×n pattern
✓ I/Q components show good dynamic range
✓ DC bias appears normal
✓ Data contains non-zero values
```

### Troubleshooting

**USB Permission Errors:**
```bash
# Check if user is in plugdev group
groups $USER

# If not listed, add and reboot
sudo usermod -a -G plugdev $USER
sudo reboot
```

**RTL-SDR Not Detected:**
```bash
# Check USB connection
lsusb | grep Realtek

# Should show something like:
# Bus 001 Device 004: ID 0bda:2838 Realtek Semiconductor Corp. RTL2838
```

**Sample Loss Issues:**
```bash
# Check system load during collection
htop

# Ensure no other processes using RTL-SDR
sudo killall rtl_sdr rtl_tcp rtl_fm
```

**Clock Synchronization:**
```bash
# Check NTP status
ntpq -p

# Force time sync if needed
sudo ntpdate -s time.nist.gov
```

## Network Setup

### Firewall Configuration

If using firewall, allow SSH and any custom ports:

```bash
# UFW example
sudo ufw allow ssh
sudo ufw allow from trusted.network.ip
sudo ufw enable
```

### Remote Access

Set up SSH keys for easy access to all stations:

```bash
# Generate key pair (if not already done)
ssh-keygen -t rsa -b 4096

# Copy to each station
ssh-copy-id user@station-ip
```

## Next Steps

After successful installation:

1. **Coordinate Test:** Run simultaneous collections on all stations
2. **Data Analysis:** Use reader to validate all collections
3. **Pattern Development:** Begin working on FM audio pattern matching
4. **Field Testing:** Start with known target locations for validation

See [Usage Guide](usage.md) for operational procedures.