#!/bin/bash

# TDOA Collector Deployment Script
# Deploys librtlsdr-2freq and collector to Raspberry Pi stations

set -e  # Exit on any error

echo "=== TDOA Collector Deployment Script ==="
echo "This script will:"
echo "1. Install build dependencies"
echo "2. Build librtlsdr-2freq with dual-frequency support"
echo "3. Build and deploy the collector"
echo "4. Set up proper permissions"
echo ""

# Check if we're on a Raspberry Pi
if grep -q "Raspberry Pi" /proc/cpuinfo 2>/dev/null; then
    echo "✓ Detected Raspberry Pi hardware"
    PI_DETECTED=true
else
    echo "⚠ Not detected as Raspberry Pi - continuing anyway"
    PI_DETECTED=false
fi

# Function to detect available CPU cores for parallel builds
detect_cores() {
    CORES=$(nproc)
    if [ "$CORES" -gt 4 ]; then
        BUILD_JOBS=4
    elif [ "$CORES" -gt 2 ]; then
        BUILD_JOBS=$CORES
    else
        BUILD_JOBS=2  # Conservative for older Pi models
    fi
    echo "Using $BUILD_JOBS parallel build jobs"
}

# Step 1: Install dependencies
echo ""
echo "=== Step 1: Installing dependencies ==="
sudo apt update
sudo apt install -y build-essential cmake git libusb-1.0-0-dev pkg-config golang-go

# Optional: Install additional Pi-specific optimizations
if [ "$PI_DETECTED" = true ]; then
    echo "Installing Raspberry Pi specific packages..."
    sudo apt install -y libraspberrypi-dev
fi

# Step 2: Build librtlsdr-2freq
echo ""
echo "=== Step 2: Building librtlsdr-2freq ==="

# Clean up any existing build
if [ -d "librtlsdr-2freq" ]; then
    echo "Removing existing librtlsdr-2freq directory..."
    rm -rf librtlsdr-2freq
fi

# Clone and build
git clone https://github.com/DC9ST/librtlsdr-2freq.git
cd librtlsdr-2freq

mkdir build
cd build

echo "Configuring build with cmake (enabling automatic kernel driver detaching)..."
cmake .. -DDETACH_KERNEL_DRIVER=ON

detect_cores
echo "Building with make -j$BUILD_JOBS..."
make -j$BUILD_JOBS

# Verify build
if [ ! -f "src/rtl_sdr" ]; then
    echo "✗ Error: rtl_sdr binary not found after build"
    exit 1
fi

echo "✓ librtlsdr-2freq built successfully"

# Test dual-frequency support
echo "Testing dual-frequency parameters..."
if ./src/rtl_sdr --help 2>&1 | grep -q "\-f.*\-h"; then
    echo "✓ Dual-frequency support confirmed"
else
    echo "✗ Warning: Dual-frequency parameters not detected in help output"
fi

cd ../..

# Step 3: Build collector
echo ""
echo "=== Step 3: Building collector ==="

# Check if Go files exist
if [ ! -f "collector.go" ]; then
    echo "✗ Error: collector.go not found in current directory"
    echo "Please run this script from the directory containing collector.go"
    exit 1
fi

# Initialize Go module if needed
if [ ! -f "go.mod" ]; then
    echo "Initializing Go module..."
    go mod init tdoa-collector
fi

# Build collector
echo "Building collector..."
go build -o collector collector.go

if [ ! -f "collector" ]; then
    echo "✗ Error: collector binary not created"
    exit 1
fi

echo "✓ Collector built successfully"

# Build reader
if [ -f "reader.go" ]; then
    echo "Building reader..."
    go build -o reader reader.go
    echo "✓ Reader built successfully"
fi

# Step 4: Set up permissions
echo ""
echo "=== Step 4: Setting up permissions ==="

# Add user to plugdev group for USB device access
echo "Adding user to plugdev group..."
sudo usermod -a -G plugdev $USER

# Install udev rules for RTL-SDR
echo "Installing RTL-SDR udev rules..."
if [ -f "librtlsdr-2freq/rtl-sdr.rules" ]; then
    sudo cp librtlsdr-2freq/rtl-sdr.rules /etc/udev/rules.d/
    sudo udevadm control --reload-rules
    sudo udevadm trigger
    echo "✓ Udev rules installed"
fi

# Make binaries executable
chmod +x collector
chmod +x librtlsdr-2freq/build/src/rtl_sdr
if [ -f "reader" ]; then
    chmod +x reader
fi

# Step 5: Update collector configuration
echo ""
echo "=== Step 5: Updating collector configuration ==="

# Create a local version of collector with correct rtl_sdr path
RTL_SDR_PATH="$(pwd)/librtlsdr-2freq/build/src/rtl_sdr"
echo "RTL-SDR binary location: $RTL_SDR_PATH"

# Check if we need to update the path in collector.go
if grep -q 'rtlSdrPath := "src/rtl_sdr"' collector.go; then
    echo "⚠ Note: collector.go still has relative path 'src/rtl_sdr'"
    echo "   Consider updating to absolute path: $RTL_SDR_PATH"
    echo "   Or ensure collector is run from the correct directory"
fi

# Step 6: Final validation
echo ""
echo "=== Step 6: Final validation ==="

echo "Testing collector help..."
if ./collector 2>&1 | grep -q "Usage:"; then
    echo "✓ Collector executable and shows help"
else
    echo "✗ Error: Collector not responding properly"
fi

echo "Testing RTL-SDR detection..."
if timeout 5 ./librtlsdr-2freq/build/src/rtl_sdr --help >/dev/null 2>&1; then
    echo "✓ RTL-SDR binary responds"
    
    # Test kernel driver detaching if RTL-SDR hardware is present
    echo "Testing automatic kernel driver detaching..."
    if timeout 3 ./librtlsdr-2freq/build/src/rtl_sdr -f 100000000 -h 101000000 -g 20 -n 100 /tmp/test.dat >/dev/null 2>&1; then
        echo "✓ Automatic kernel driver detaching works"
        rm -f /tmp/test.dat
    else
        echo "⚠ Could not test kernel driver detaching (RTL-SDR may not be connected)"
    fi
else
    echo "⚠ RTL-SDR binary test timed out (normal if no hardware connected)"
fi

# Display final information
echo ""
echo "=== Deployment Complete ==="
echo ""
echo "Binaries created:"
echo "  - ./collector (TDOA data collector)"
if [ -f "reader" ]; then
echo "  - ./reader (data validation tool)"
fi
echo "  - ./librtlsdr-2freq/build/src/rtl_sdr (dual-frequency RTL-SDR)"
echo ""
echo "Usage example:"
echo "  future_time=\$((\$(date +%s) + 30))"
echo "  ./collector --duration=10 96900000 162550000 \$future_time kx0u"
echo "  ./reader kx0u-\$future_time.dat 10"
echo ""
echo "⚠ IMPORTANT: You may need to log out and back in for group changes to take effect"
echo "⚠ If USB permissions fail, try: sudo reboot"
echo ""

if [ "$PI_DETECTED" = true ]; then
    echo "Raspberry Pi specific notes:"
    echo "- Build completed with $BUILD_JOBS parallel jobs"
    echo "- Consider increasing GPU memory split if needed: sudo raspi-config"
    echo "- For best performance, use a fast SD card (Class 10 or better)"
fi

echo "Deployment script completed successfully!"