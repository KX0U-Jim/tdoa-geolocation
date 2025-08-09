# FM Audio Pattern Matching for TDOA

## The Innovation

Traditional TDOA systems avoid FM broadcast signals because their modulation is "unpredictable." This project takes the opposite approach: **use the complexity of FM modulation as the solution**.

## Technical Approach

### The Problem with FM Signals
- Constantly changing audio content
- Dynamic compression and processing
- Variable modulation depth
- No predictable timing markers

### Our Solution: Pre-recorded Audio Prediction

Instead of trying to work around FM's variability, we embrace it:

1. **Record the FM audio** from the reference station
2. **Predict the RF signal** that audio would generate
3. **Search for those patterns** in the collected RF samples
4. **Extract timing differences** from pattern matches

## Implementation Details

### Phase 1: Audio Recording and Analysis

```bash
# Record FM audio content (example using external tools)
# This gives us the baseband audio that will modulate the carrier

# Example audio characteristics:
# - Sample rate: 44.1 kHz or 48 kHz
# - Duration: Match your collection window
# - Format: WAV, uncompressed
```

### Phase 2: RF Pattern Prediction

The audio modulates the FM carrier according to:
```
f_instantaneous = f_carrier + k_f * audio_sample
```

Where:
- `f_carrier` = Reference station frequency (e.g., 96.9 MHz)
- `k_f` = Frequency deviation constant (~75 kHz max for WBFM)
- `audio_sample` = Instantaneous audio amplitude

**FFT-based Approach:**
```
For each audio frame:
1. FFT the audio to get frequency components
2. Map audio frequencies to FM deviation
3. Generate expected RF spectrum at that moment
4. Create searchable "fingerprint" pattern
```

### Phase 3: Pattern Matching in RF Samples

```go
// Pseudocode for pattern matching
func FindAudioPattern(rfSamples []complex64, expectedPattern []complex64) float64 {
    // Cross-correlate expected pattern with received RF
    correlation := CrossCorrelate(rfSamples, expectedPattern)
    
    // Find peak correlation
    peakIndex := FindMaxPeak(correlation)
    
    // Return timing offset
    return float64(peakIndex) / sampleRate
}
```

## Advantages Over Traditional Methods

### Compared to GPS Disciplined Oscillators
- **Cost:** $50 RTL-SDR vs $500+ GPSDO
- **Availability:** FM stations everywhere vs GPS timing requirements
- **Simplicity:** No external timing infrastructure needed

### Compared to Other Reference Signals
- **Coverage:** FM stations have wide, known coverage areas
- **Power:** High-power transmitters ensure good SNR at all collectors
- **Stability:** Professional broadcast equipment with stable carriers

### Compared to Avoiding FM Entirely
- **Target Availability:** Many targets of interest use FM (amateur, public safety)
- **Reference Availability:** Abundant FM broadcast stations with known locations

## Pattern Matching Algorithms

### Cross-Correlation Method
```
correlation[n] = Σ(received[m] * conj(expected[m-n]))
```

### Sliding Window Approach
- Process RF samples in overlapping windows
- Match predicted patterns against each window
- Track correlation peaks over time

### Multi-Scale Analysis
- Match patterns at different time scales
- Short patterns: Individual audio transients (syllables, notes)
- Long patterns: Phrase-level or sentence-level content

## Implementation Challenges

### Synchronization
- **Problem:** Don't know exact timing relationship between audio recording and RF collection
- **Solution:** Search for pattern matches across entire collection window

### Doppler Effects
- **Problem:** Mobile targets may have frequency shifts
- **Solution:** Search with frequency offsets, or limit to stationary targets

### Multipath
- **Problem:** Reflected signals create false correlation peaks
- **Solution:** Use strongest correlation peak, or multiple peak analysis

### Audio Processing Variations
- **Problem:** Station audio processing may change
- **Solution:** Use robust pattern features, or update audio recordings

## Performance Characteristics

### Timing Precision
- **Theoretical:** Limited by sample rate (e.g., 500ns at 2 Msps)
- **Practical:** Depends on pattern sharpness and SNR
- **Target:** Sub-microsecond accuracy for meter-level geolocation

### Processing Requirements
- **FFT Operations:** Real-time capable on modern hardware
- **Memory Usage:** Pattern libraries scale with audio content
- **Correlation:** Parallelizable across multiple cores

## Future Enhancements

### Machine Learning
- Train models to recognize optimal correlation patterns
- Automatic pattern selection from audio content
- Adaptive filtering for varying RF conditions

### Real-Time Implementation
- Streaming audio capture from FM receivers
- Live pattern generation and matching
- Dynamic pattern library updates

### Multi-Station Pattern Libraries
- Share successful patterns across team members
- Build database of effective reference signals
- Collaborative pattern validation

## Validation Methods

### Simulation Testing
- Generate synthetic FM signals from known audio
- Test pattern matching accuracy in controlled environment
- Validate before field deployment

### Cross-Validation
- Use multiple FM stations as references
- Compare timing results between different patterns
- Statistical validation of consistency

### Ground Truth Testing
- Use targets at known locations for validation
- Compare computed location with actual position
- Measure and improve accuracy over time

---

This innovation transforms FM broadcast signals from a problem into a powerful solution for TDOA geolocation, making the system practical and cost-effective for amateur radio applications.