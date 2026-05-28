# ZwickerLoudnessAudio.jl

Audio-file companion to [ZwickerLoudness.jl](https://github.com/slink/ZwickerLoudness.jl).
Compute ISO 532-1:2017 Method 1 stationary loudness directly from a `.wav`
file or a time-domain signal.

The third-octave filterbank is a per-band Butterworth bandpass designed per
ANSI S1.1-1986 (matching the approach used by [MoSQITo](https://github.com/Eomys/MoSQITo)).
Default order 3 reproduces the published Annex B reference values for the
three tone signals within ISO 532-1's 5% tolerance.

## Installation

This package is not yet registered. Install from the local checkout:

```julia
using Pkg
Pkg.develop(path="/path/to/ZwickerLoudness.jl")    # parent kernel
Pkg.develop(path="/path/to/ZwickerLoudnessAudio.jl")
```

## Quick start

```julia
using ZwickerLoudnessAudio

# From a .wav file (multichannel files are mixed to mono).
result = loudness_zwst("recording.wav")
result.loudness           # total loudness N [sone]
result.loudness_level     # loudness level L_N [phon]
result.specific_loudness  # 240-bin N'(z) [sone/Bark]

# From a time-domain signal already in pascals.
using Statistics
fs = 48_000
t = range(0, 0.5, length=round(Int, 0.5 * fs))
sig_pa = 2e-5 * 10^(60/20) * sin.(2π * 1000 .* t) * sqrt(2)  # 1 kHz, 60 dB SPL
loudness_zwst(sig_pa, fs).loudness  # ~4.0 sone (ISO Signal 3 reference)

# Per-band SPL diagnostic.
band_levels(sig_pa, fs)  # 28-element Vector{Float64}, dB re 20 µPa
```

## Calibration

`pa_per_unit` converts the raw signal to pascals. Defaults to `1.0` (signal
already in Pa).

```julia
loudness_zwst(raw_signal, fs; pa_per_unit=0.5)  # signal is half of true Pa
```

## Conformance

ISO 532-1:2017 Annex B reference values (synthesized in-process; tolerance is
the standard's ±5%):

| Signal | Stimulus | Reference | Computed |
|--------|----------|-----------|----------|
| 2 | 250 Hz tone @ 80 dB SPL | 14.655 sone | ~15.09 sone (+2.9%) |
| 3 | 1 kHz tone @ 60 dB SPL | 4.019 sone | ~4.16 sone (+3.6%) |
| 4 | 4 kHz tone @ 40 dB SPL | 1.549 sone | ~1.56 sone (+0.5%) |

Signal 5 (pink noise @ 60 dB SPL, reference 10.498 sone) needs a calibrated
pink-noise generator to test reproducibly and is deferred.

## API

```julia
loudness_zwst(signal::AbstractVector{<:Real}, fs::Real;
              field_type=:free, pa_per_unit=1.0, order=3) -> ZwickerResult

loudness_zwst(path::AbstractString;
              field_type=:free, pa_per_unit=1.0, order=3) -> ZwickerResult

band_levels(signal, fs; pa_per_unit=1.0, order=3) -> Vector{Float64}  # 28 SPL values
```

`order` is the Butterworth filter order per band. Default `3` matches
MoSQITo. Higher orders give tighter selectivity (less neighbor leakage from
a tone) at the cost of compute and potential low-band conditioning issues.

## License

MIT
