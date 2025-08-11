#!/bin/bash

# RTL-SDR DVB Module Killer
echo "=== RTL-SDR DVB Module Removal ==="

echo "Before - loaded modules:"
lsmod | grep -i "rtl\|dvb" || echo "No modules found"
echo ""

echo "Removing DVB and RTL kernel modules in correct order..."

# Remove in dependency order (most dependent first)
sudo modprobe -r rtl2832_sdr 2>/dev/null && echo "✓ Removed rtl2832_sdr" || echo "✗ rtl2832_sdr removal failed"
sudo modprobe -r dvb_usb_rtl28xxu 2>/dev/null && echo "✓ Removed dvb_usb_rtl28xxu" || echo "✗ dvb_usb_rtl28xxu removal failed"
sudo modprobe -r rtl2832 2>/dev/null && echo "✓ Removed rtl2832" || echo "✗ rtl2832 removal failed"
sudo modprobe -r dvb_usb_v2 2>/dev/null && echo "✓ Removed dvb_usb_v2" || echo "✗ dvb_usb_v2 removal failed"
sudo modprobe -r dvb_core 2>/dev/null && echo "✓ Removed dvb_core" || echo "✗ dvb_core removal failed"

# Also remove any video/media modules that might interfere
sudo modprobe -r videobuf2_vmalloc 2>/dev/null && echo "✓ Removed videobuf2_vmalloc" || echo "✗ videobuf2_vmalloc removal failed"
sudo modprobe -r i2c_mux 2>/dev/null && echo "✓ Removed i2c_mux" || echo "✗ i2c_mux removal failed"

echo ""
echo "After - remaining modules:"
lsmod | grep -i "rtl\|dvb" || echo "No RTL/DVB modules remain (GOOD!)"
echo ""

echo "Testing RTL-SDR access..."
if [ -f "./librtlsdr-2freq/build/src/rtl_sdr" ]; then
    timeout 5 ./librtlsdr-2freq/build/src/rtl_sdr -f 100000000 -h 101000000 -g 20 -n 100 /tmp/rtl_test.dat 2>&1
    RESULT=$?
    rm -f /tmp/rtl_test.dat 2>/dev/null
    
    if [ $RESULT -eq 0 ]; then
        echo ""
        echo "✅ SUCCESS: RTL-SDR is now accessible!"
        echo "You can now run collector operations."
    else
        echo ""
        echo "❌ FAILED: RTL-SDR still not accessible (exit code: $RESULT)"
        echo ""
        echo "USB device check:"
        lsusb | grep -i "rtl\|realtek"
        echo ""
        echo "Remaining troubleshooting:"
        echo "1. Try a different USB port"
        echo "2. Reboot the system"
        echo "3. Check if user is in plugdev group: groups \$USER"
    fi
else
    echo "❌ ERROR: librtlsdr-2freq binary not found at ./librtlsdr-2freq/build/src/rtl_sdr"
fi