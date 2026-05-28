# ZwickerLoudnessAudio.jl

Audio-file companion to [ZwickerLoudness.jl](https://github.com/slink/ZwickerLoudness.jl).
Compute ISO 532-1:2017 Method 1 stationary loudness directly from a `.wav`
file or a time-domain signal.

> **Status: smoke scaffold.** The third-octave filterbank is currently a
> simple FFT-power integrator and is **not** IEC 61260 compliant. Numerical
> conformance against ISO 532-1 Annex B Signals 3 and 4 is deferred to a
> follow-up that lands a proper IIR filterbank. Use this package now for
> wiring, exploration, and rough estimates; not yet for certified loudness
> reporting.

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
loudness_zwst(sig_pa, fs).loudness  # ~3.5 sone (smoke filterbank)

# Per-band SPL diagnostic.
band_levels(sig_pa, fs)  # 28-element Vector{Float64}, dB re 20 µPa
```

## Calibration

`pa_per_unit` converts the raw signal to pascals. Defaults to `1.0` (signal
already in Pa).

```julia
loudness_zwst(raw_signal, fs; pa_per_unit=0.5)  # signal is half of true Pa
```

## Why the IIR filterbank matters

ISO 532-1 Annex B Signals 2-5 are `.wav` files; the standard validates
implementations on the full pipeline `wav → IEC 61260 filterbank → band SPL →
Method 1 kernel`. The smoke FFT integrator here gets band SPL roughly right
for stationary tones, but it misses real-filterbank effects (passband ripple,
finite skirt rejection, time-frequency tradeoff) that materially affect mid-
and high-frequency tones. See `docs/specs/2026-05-28-scaffold-design.md` for
the planned upgrade path.

## API

```julia
loudness_zwst(signal::AbstractVector{<:Real}, fs::Real;
              field_type=:free, pa_per_unit=1.0) -> ZwickerResult

loudness_zwst(path::AbstractString;
              field_type=:free, pa_per_unit=1.0) -> ZwickerResult

band_levels(signal, fs; pa_per_unit=1.0) -> Vector{Float64}  # 28 SPL values
```

## License

MIT
