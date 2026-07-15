# ZwickerLoudnessAudio.jl

[![CI](https://github.com/slink/ZwickerLoudnessAudio.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/slink/ZwickerLoudnessAudio.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/slink/ZwickerLoudnessAudio.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/slink/ZwickerLoudnessAudio.jl)
[![DOI](https://zenodo.org/badge/1253784080.svg)](https://zenodo.org/badge/latestdoi/1253784080)

Audio-file companion to [ZwickerLoudness.jl](https://github.com/slink/ZwickerLoudness.jl).
Compute ISO 532-1:2017 loudness directly from a `.wav` file or a
time-domain signal: the stationary method (clause 5; `loudness_zwst`) and
the time-varying method (clause 6; `loudness_zwtv`).

For the stationary method, the third-octave filterbank is a per-band
Butterworth bandpass designed per ANSI S1.1-1986 (matching the approach
used by [MoSQITo](https://github.com/Eomys/MoSQITo)).
Default order 3 reproduces all four ISO 532-1 Annex B reference signals
(two tones, one broadband tone, and pink noise) within the standard's ±5%
tolerance. The time-varying method instead uses the standard's own Annex A
filter design (see below).

## Installation

```julia
using Pkg
Pkg.add("ZwickerLoudnessAudio")
```

The parent kernel [ZwickerLoudness.jl](https://github.com/slink/ZwickerLoudness.jl)
is installed automatically as a dependency.

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

## Time-Varying Loudness (ISO 532-1:2017 clause 6)

`loudness_zwtv` computes time-varying loudness per ISO 532-1:2017
clause 6, "Method for time-varying sounds" (colloquially "Method 2" —
this project's own shorthand, not the standard's terminology; **not to be
confused with ISO 532-2**, the unrelated Moore–Glasberg standard). It
returns a `ZwickerTimeVaryingResult` from
[ZwickerLoudness.jl](https://github.com/slink/ZwickerLoudness.jl): the
loudness time series N(t) on a 2 ms axis, the percentile loudness N5/N10
(loudness reached or exceeded 5 % / 10 % of the time, ISO/DIN definition,
Hyndman–Fan type-7 quantiles), and the specific loudness N'(z, t).

```julia
using ZwickerLoudnessAudio
fs = 48_000
t = range(0, 0.3, length=round(Int, 0.3 * fs))
env = @. 0.5 * (1 - cos(2π * 4 * t))            # 4 Hz raised-cosine envelope
sig_pa = @. 2e-5 * 10^(70/20) * sqrt(2) * env * sin(2π * 1000 * t)
result = loudness_zwtv(sig_pa, fs)
length(result.loudness_over_time)  # 150 — N(t) [sone] on the 2 ms axis
result.N5                          # 6.733 [sone]
result.N10                         # 6.629 [sone]
size(result.specific_loudness)     # (240, 150) — N'(z, t) [sone/Bark]

loudness_zwtv("recording.wav")               # wav entry, mirrors loudness_zwst
loudness_zwtv(stereo_matrix, fs; channel=2)  # matrix entry, ditto
```

The band stage is a faithful transcription of MoSQITo's
`_third_octave_levels`/`_square_and_smooth` at the same pinned commit the
kernel transcribes (ISO 532-1:2017 Annex A three-section filter cascade,
squaring, triple low-pass smoothing, 0.5 ms band-level series) — it is
deliberately **not** the ANSI Butterworth bank `loudness_zwst` uses,
because the time-varying method specifies its own filter design.

**48 kHz only:** `loudness_zwtv` throws an `ArgumentError` for any other
sample rate (the reference front end hard-requires 48 kHz; no resampler
is reimplemented). Resample first, e.g. `DSP.resample(signal, 48_000/fs)`,
and retry. There is no auto-resample fallback like `loudness_zwst`'s.

Cross-checked end-to-end against MoSQITo `loudness_zwtv` on identical
48 kHz signals (steady tone, tone burst, AM tone): N(t) agrees to ~3e-16
relative, N5/N10 to ~2e-15 absolute. An ISO 532-1 Annex B.4 end-to-end
conformance test runs whenever the (never vendored, ISO-copyrighted)
reference recording is available locally and self-skips otherwise — see
`test/fixtures/NOTICE.md` for attribution and licensing posture.

## Calibration

`pa_per_unit` converts the raw signal to pascals. Defaults to `1.0` (signal
already in Pa).

```julia
loudness_zwst(raw_signal, fs; pa_per_unit=0.5)  # signal is half of true Pa
```

## Sample rate

For `loudness_zwst`: inputs at 32 kHz or above pass through to the
filterbank untouched. Lower rates are polyphase-resampled to 48 kHz via
`DSP.resample`, with a `@warn` that bands above the original Nyquist will
be near-silent (upsampling can't recover bandwidth the signal never had).
High rates up to at least 384 kHz are also supported directly — the
SOS-form filter cascade stays well-conditioned for the 25 Hz band even at
very small normalized edges.

For `loudness_zwtv`: exactly 48 kHz, no exceptions — see the time-varying
section above.

## Conformance

Stationary method — ISO 532-1:2017 Annex B reference values (tone signals
synthesized in-process; pink noise generated via seeded FFT 1/√f shaping.
Tolerance is the standard's ±5%):

| Signal | Stimulus | Reference | Computed |
|--------|----------|-----------|----------|
| 2 | 250 Hz tone @ 80 dB SPL | 14.655 sone | ~15.09 sone (+2.9%) |
| 3 | 1 kHz tone @ 60 dB SPL | 4.019 sone | ~4.16 sone (+3.6%) |
| 4 | 4 kHz tone @ 40 dB SPL | 1.549 sone | ~1.56 sone (+0.5%) |
| 5 | pink noise @ 60 dB SPL | 10.498 sone | ~10.77 sone (+2.6%) |

Time-varying method — MoSQITo crosscheck fixtures on synthesized signals
(vendored, always tested) plus an Annex B.4 end-to-end conformance test
that self-skips when the local (never vendored) ISO reference recording is
absent; see the time-varying section above for measured agreement.

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

loudness_zwtv(signal::AbstractVector{<:Real}, fs::Real;
              field_type=:free, pa_per_unit=1.0) -> ZwickerTimeVaryingResult

loudness_zwtv(samples::AbstractMatrix{<:Real}, fs::Real;
              channel=:mono, field_type=:free,
              pa_per_unit=1.0) -> ZwickerTimeVaryingResult

loudness_zwtv(path::AbstractString;
              channel=:mono, field_type=:free,
              pa_per_unit=1.0) -> ZwickerTimeVaryingResult
```

`channel=:mono` (default) averages all channels into one signal; an integer
picks one 1-based channel.

`order` (`loudness_zwst`/`band_levels` only) is the Butterworth filter
order per band. Default `3` matches MoSQITo. Higher orders give tighter
selectivity (less neighbor leakage from a tone) at the cost of compute and
potential low-band conditioning issues. `loudness_zwtv`'s filterbank is
fixed by the standard and takes no `order`.

## License

MIT. The time-varying band-stage filter tables and algorithm in
`src/zwtv_frontend.jl` are transcribed from
[MoSQITo](https://github.com/Eomys/MoSQITo) (Apache License 2.0) — see
`test/fixtures/NOTICE.md` for full attribution.
