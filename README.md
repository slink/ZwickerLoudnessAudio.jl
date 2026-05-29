# ZwickerLoudnessAudio.jl

[![CI](https://github.com/slink/ZwickerLoudnessAudio.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/slink/ZwickerLoudnessAudio.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/slink/ZwickerLoudnessAudio.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/slink/ZwickerLoudnessAudio.jl)

Audio-file companion to [ZwickerLoudness.jl](https://github.com/slink/ZwickerLoudness.jl).
Compute ISO 532-1:2017 Method 1 stationary loudness directly from a `.wav`
file or a time-domain signal.

The third-octave filterbank is a per-band Butterworth bandpass designed per
ANSI S1.1-1986 (matching the approach used by [MoSQITo](https://github.com/Eomys/MoSQITo)).
Default order 3 reproduces all four ISO 532-1 Annex B reference signals
(two tones, one broadband tone, and pink noise) within the standard's ±5%
tolerance.

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

# From a .wav file. Multichannel files are averaged to mono by default;
# pass `channel=n` to pick one 1-based channel instead.
result = loudness_zwst("recording.wav")
result.loudness           # total loudness N [sone]
result.loudness_level     # loudness level L_N [phon]
result.specific_loudness  # 240-bin N'(z) [sone/Bark]

loudness_zwst("stereo.wav"; channel=1)  # left channel only

# From a time-domain signal already in pascals.
using Statistics
fs = 48_000
t = range(0, 0.5, length=round(Int, 0.5 * fs))
sig_pa = 2e-5 * 10^(60/20) * sin.(2π * 1000 .* t) * sqrt(2)  # 1 kHz, 60 dB SPL
loudness_zwst(sig_pa, fs).loudness  # ~4.0 sone (ISO Signal 3 reference)

# From a samples × channels matrix.
loudness_zwst(stereo_matrix, fs; channel=2)

# Per-band SPL diagnostic.
band_levels(sig_pa, fs)  # 28-element Vector{Float64}, dB re 20 µPa
```

## Calibration

`pa_per_unit` converts the raw signal to pascals. Defaults to `1.0` (signal
already in Pa).

```julia
loudness_zwst(raw_signal, fs; pa_per_unit=0.5)  # signal is half of true Pa
```

## Sample rate

Inputs at 32 kHz or above pass through to the filterbank untouched. Lower
rates are polyphase-resampled to 48 kHz via `DSP.resample`, with a `@warn`
that bands above the original Nyquist will be near-silent (upsampling can't
recover bandwidth the signal never had). High rates up to at least 384 kHz
are also supported directly — the SOS-form filter cascade stays
well-conditioned for the 25 Hz band even at very small normalized edges.

## Conformance

ISO 532-1:2017 Annex B reference values (tone signals synthesized
in-process; pink noise generated via seeded FFT 1/√f shaping. Tolerance is
the standard's ±5%):

| Signal | Stimulus | Reference | Computed |
|--------|----------|-----------|----------|
| 2 | 250 Hz tone @ 80 dB SPL | 14.655 sone | ~15.09 sone (+2.9%) |
| 3 | 1 kHz tone @ 60 dB SPL | 4.019 sone | ~4.16 sone (+3.6%) |
| 4 | 4 kHz tone @ 40 dB SPL | 1.549 sone | ~1.56 sone (+0.5%) |
| 5 | pink noise @ 60 dB SPL | 10.498 sone | ~10.77 sone (+2.6%) |

## API

```julia
loudness_zwst(signal::AbstractVector{<:Real}, fs::Real;
              field_type=:free, pa_per_unit=1.0, order=3) -> ZwickerResult

loudness_zwst(samples::AbstractMatrix{<:Real}, fs::Real;
              channel=:mono, field_type=:free, pa_per_unit=1.0,
              order=3) -> ZwickerResult

loudness_zwst(path::AbstractString;
              channel=:mono, field_type=:free, pa_per_unit=1.0,
              order=3) -> ZwickerResult

band_levels(signal, fs; pa_per_unit=1.0, order=3) -> Vector{Float64}  # 28 SPL values
```

`channel=:mono` (default) averages all channels into one signal; an integer
picks one 1-based channel.

`order` is the Butterworth filter order per band. Default `3` matches
MoSQITo. Higher orders give tighter selectivity (less neighbor leakage from
a tone) at the cost of compute and potential low-band conditioning issues.

## License

MIT
