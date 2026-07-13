# =========================================================================== #
#  ISO 532-1:2017 clause 6, "Method for time-varying sounds"
#  (colloquially "Method 2" -- NOT to be confused with ISO 532-2, the
#  unrelated Moore-Glasberg standard).
#
#  Front-end stage: 48 kHz time-domain signal -> 28-band, 2000 Hz (0.5 ms)
#  third-octave level time series, the input contract of
#  ZwickerLoudness.jl's `zwicker_loudness_time_varying` kernel (v0.3.0).
#
#  Faithful transcription of MoSQITo's `_third_octave_levels.py` and
#  `_square_and_smooth.py` @ d990c33f94f1 (Apache-2.0) -- see
#  test/fixtures/NOTICE.md. Deliberately NOT the ANSI S1.1-1986
#  Butterworth bank used by this package's Method 1 `band_levels`: the
#  reference algorithm's Method 2 band stage is a different filter design
#  (ISO 532-1:2017 Annex A 6th-order per-band SOS bank) at a different
#  output rate (2 ms internal write-out here, 0.5 ms actual before
#  decimation), and fidelity to the pinned reference governs.
# =========================================================================== #

export loudness_zwtv

const _ZWTV_FS = 48000.0
const _ZWTV_DEC_FACTOR = 24            # Int(fs / 2000) at fs = 48 kHz.
const _ZWTV_N_BANDS = 28
const _ZWTV_TINY_VALUE = 1e-12
const _ZWTV_I_REF = 4e-10
# Kernel-side decimation from the 0.5 ms band rate to the 2 ms output axis
# (loudness_zwtv.py:128-131 @ d990c33f94f1) -- a SEPARATE constant from
# `_ZWTV_DEC_FACTOR` above (front-end 48 kHz -> 2 kHz), used here only to
# decimate MoSQITo's own linspace time axis to match the kernel's output
# length (see the `time_axis` note in `loudness_zwtv`'s docstring).
const _ZWTV_KERNEL_DEC_FACTOR = 4

# ISO 532-1:2017 Annex A one-third-octave filter design: three cascaded
# second-order sections (b0,b1,b2,a0,a1,a2 per row) shared by every band,
# perturbed per band by `ZWTV_THIRD_OCTAVE_FILTER_DIFF` below.
# Transcribed from `_third_octave_levels.py:41-43` @ d990c33f94f1.
const ZWTV_THIRD_OCTAVE_FILTER_REF = Float64[
  1 2 1 1 -2 1;
  1 0 -1 1 -2 1;
  1 -2 1 1 -2 1
]

# Per-band correction to `ZWTV_THIRD_OCTAVE_FILTER_REF` (only the a1/a2
# columns differ from zero): `coeff = REF - DIFF[i]` gives the actual
# per-band SOS coefficients. 28 bands x 3 sections x 6 coefficients.
# Transcribed from `_third_octave_levels.py:47-190` @ d990c33f94f1
# (ISO 532-1:2017 Table A.2).
const ZWTV_THIRD_OCTAVE_FILTER_DIFF = Matrix{Float64}[
    [
      0.0 0.0 0.0 0.0 -0.00067026 0.000659453;
      0.0 0.0 0.0 0.0 -0.000375071 0.000361926;
      0.0 0.0 0.0 0.0 -0.000306523 0.000297634
    ],
    [
      0.0 0.0 0.0 0.0 -0.000847258 0.000830131;
      0.0 0.0 0.0 0.0 -0.000476448 0.000455616;
      0.0 0.0 0.0 0.0 -0.000388773 0.000374685
    ],
    [
      0.0 0.0 0.0 0.0 -0.0010721 0.00104496;
      0.0 0.0 0.0 0.0 -0.000606567 0.000573553;
      0.0 0.0 0.0 0.0 -0.000494004 0.000471677
    ],
    [
      0.0 0.0 0.0 0.0 -0.00135836 0.00131535;
      0.0 0.0 0.0 0.0 -0.000774327 0.000722007;
      0.0 0.0 0.0 0.0 -0.000629154 0.000593771
    ],
    [
      0.0 0.0 0.0 0.0 -0.0017238 0.00165564;
      0.0 0.0 0.0 0.0 -0.00099178 0.000908866;
      0.0 0.0 0.0 0.0 -0.000803529 0.000747455
    ],
    [
      0.0 0.0 0.0 0.0 -0.00219188 0.00208388;
      0.0 0.0 0.0 0.0 -0.00127545 0.00114406;
      0.0 0.0 0.0 0.0 -0.00102976 0.0009409
    ],
    [
      0.0 0.0 0.0 0.0 -0.00279386 0.00262274;
      0.0 0.0 0.0 0.0 -0.00164828 0.00144006;
      0.0 0.0 0.0 0.0 -0.0013252 0.00118438
    ],
    [
      0.0 0.0 0.0 0.0 -0.00357182 0.00330071;
      0.0 0.0 0.0 0.0 -0.00214252 0.00181258;
      0.0 0.0 0.0 0.0 -0.00171397 0.00149082
    ],
    [
      0.0 0.0 0.0 0.0 -0.00458305 0.00415355;
      0.0 0.0 0.0 0.0 -0.00280413 0.00228135;
      0.0 0.0 0.0 0.0 -0.00223006 0.00187646
    ],
    [
      0.0 0.0 0.0 0.0 -0.00590655 0.00522622;
      0.0 0.0 0.0 0.0 -0.00369947 0.00287118;
      0.0 0.0 0.0 0.0 -0.00292205 0.00236178
    ],
    [
      0.0 0.0 0.0 0.0 -0.00765243 0.00657493;
      0.0 0.0 0.0 0.0 -0.0049254 0.00361318;
      0.0 0.0 0.0 0.0 -0.00386007 0.0029724
    ],
    [
      0.0 0.0 0.0 0.0 -0.0100023 0.0082961;
      0.0 0.0 0.0 0.0 -0.00663788 0.00455999;
      0.0 0.0 0.0 0.0 -0.00515982 0.00375306
    ],
    [
      0.0 0.0 0.0 0.0 -0.013123 0.010422;
      0.0 0.0 0.0 0.0 -0.00902274 0.00573132;
      0.0 0.0 0.0 0.0 -0.00694543 0.00471734
    ],
    [
      0.0 0.0 0.0 0.0 -0.0173693 0.0130947;
      0.0 0.0 0.0 0.0 -0.0124176 0.00720526;
      0.0 0.0 0.0 0.0 -0.00946002 0.00593145
    ],
    [
      0.0 0.0 0.0 0.0 -0.0231934 0.0164308;
      0.0 0.0 0.0 0.0 -0.0173009 0.00904761;
      0.0 0.0 0.0 0.0 -0.0130358 0.00744926
    ],
    [
      0.0 0.0 0.0 0.0 -0.0313292 0.020637;
      0.0 0.0 0.0 0.0 -0.0244342 0.0113731;
      0.0 0.0 0.0 0.0 -0.0182108 0.00936778
    ],
    [
      0.0 0.0 0.0 0.0 -0.0428261 0.0259325;
      0.0 0.0 0.0 0.0 -0.0349619 0.0143046;
      0.0 0.0 0.0 0.0 -0.0257855 0.0117912
    ],
    [
      0.0 0.0 0.0 0.0 -0.0591733 0.0325054;
      0.0 0.0 0.0 0.0 -0.0506072 0.0179513;
      0.0 0.0 0.0 0.0 -0.0369401 0.0148094
    ],
    [
      0.0 0.0 0.0 0.0 -0.0826348 0.0405894;
      0.0 0.0 0.0 0.0 -0.0740348 0.0224476;
      0.0 0.0 0.0 0.0 -0.0534977 0.0185371
    ],
    [
      0.0 0.0 0.0 0.0 -0.117018 0.0508116;
      0.0 0.0 0.0 0.0 -0.109516 0.0281387;
      0.0 0.0 0.0 0.0 -0.0785097 0.0232872
    ],
    [
      0.0 0.0 0.0 0.0 -0.167714 0.0637872;
      0.0 0.0 0.0 0.0 -0.163378 0.0353729;
      0.0 0.0 0.0 0.0 -0.116419 0.0293723
    ],
    [
      0.0 0.0 0.0 0.0 -0.242528 0.0798576;
      0.0 0.0 0.0 0.0 -0.245161 0.044337;
      0.0 0.0 0.0 0.0 -0.173972 0.0370015
    ],
    [
      0.0 0.0 0.0 0.0 -0.353142 0.099633;
      0.0 0.0 0.0 0.0 -0.369163 0.0553535;
      0.0 0.0 0.0 0.0 -0.261399 0.0465428
    ],
    [
      0.0 0.0 0.0 0.0 -0.516316 0.124177;
      0.0 0.0 0.0 0.0 -0.555473 0.0689403;
      0.0 0.0 0.0 0.0 -0.393998 0.0586715
    ],
    [
      0.0 0.0 0.0 0.0 -0.756635 0.155023;
      0.0 0.0 0.0 0.0 -0.834281 0.0858123;
      0.0 0.0 0.0 0.0 -0.594547 0.074396
    ],
    [
      0.0 0.0 0.0 0.0 -1.10165 0.191713;
      0.0 0.0 0.0 0.0 -1.23939 0.105243;
      0.0 0.0 0.0 0.0 -0.891666 0.0940354
    ],
    [
      0.0 0.0 0.0 0.0 -1.58477 0.239049;
      0.0 0.0 0.0 0.0 -1.80505 0.128794;
      0.0 0.0 0.0 0.0 -1.325 0.121333
    ],
    [
      0.0 0.0 0.0 0.0 -2.5063 0.142308;
      0.0 0.0 0.0 0.0 -2.19464 0.27647;
      0.0 0.0 0.0 0.0 -1.90231 0.147304
    ],
]

# Overall gain applied once per band, after the 3-section cascade (not
# folded into any individual section) -- matches `filter_gain[i_bands] *
# sosfilt(coeff, sig)` in `_third_octave_levels.py:269` @ d990c33f94f1
# (ISO 532-1:2017 Table A.2).
const ZWTV_FILTER_GAIN = Float64[
    4.30764e-11,
    8.5934e-11,
    1.71424e-10,
    3.41944e-10,
    6.82035e-10,
    1.36026e-09,
    2.71261e-09,
    5.4087e-09,
    1.07826e-08,
    2.1491e-08,
    4.28228e-08,
    8.54316e-08,
    1.70009e-07,
    3.38215e-07,
    6.7199e-07,
    1.33531e-06,
    2.65172e-06,
    5.25477e-06,
    1.0378e-05,
    2.0487e-05,
    4.05198e-05,
    7.97914e-05,
    0.000156511,
    0.000304954,
    0.000599157,
    0.00116544,
    0.00227488,
    0.00391006,
]

# Build the fixed 3-section SOS filter for one-based band `band_idx`
# (1..28). `coeff = REF - DIFF[band_idx]` per band, cascaded via
# `DSP.SecondOrderSections` with the band's overall gain applied once at
# the end -- mirrors `sosfilt(coeff, sig)` then `filter_gain[i] * ...`.
function _zwtv_band_filter(band_idx::Int)
    coeff = ZWTV_THIRD_OCTAVE_FILTER_REF .- ZWTV_THIRD_OCTAVE_FILTER_DIFF[band_idx]
    biquads = [Biquad(coeff[s, 1], coeff[s, 2], coeff[s, 3], coeff[s, 4], coeff[s, 5], coeff[s, 6])
               for s in 1:3]
    return SecondOrderSections{:z}(biquads, ZWTV_FILTER_GAIN[band_idx])
end

# Computed once at load time -- fs is fixed at 48 kHz for this stage (per
# the "resample first" contract below), so there's nothing to key a cache
# on, unlike Method 1's `_band_filters`.
const ZWTV_BAND_FILTERS = [_zwtv_band_filter(i) for i in 1:_ZWTV_N_BANDS]

"""
    _third_octave_levels(sig, fs) -> (band_levels, band_time_axis)

Time-domain one-third-octave band levels at the ISO 532-1:2017 clause 6
front-end rate (2000 Hz, 0.5 ms). `sig` must be sampled at 48 kHz. Returns
`band_levels` (28 x nblocks, dB SPL re 20 µPa) -- the exact input contract
of `ZwickerLoudness.zwicker_loudness_time_varying` -- and `band_time_axis`
(nblocks,), MoSQITo's own `linspace(0, length(sig)/fs, nblocks)` time
axis (see the `loudness_zwtv` docstring for why this, not an idealized
2 ms grid, is what gets reported downstream).

Transcribed from `_third_octave_levels.py` @ d990c33f94f1 (Apache-2.0),
citing ISO 532-1:2017 section 6.3 and Annex A.
"""
function _third_octave_levels(sig::Vector{Float64}, fs::Float64)
    n = length(sig)
    n_time = length(1:_ZWTV_DEC_FACTOR:n)
    band_time_axis = collect(range(0.0, n / fs; length=n_time))

    band_levels = Matrix{Float64}(undef, _ZWTV_N_BANDS, n_time)
    for i in 1:_ZWTV_N_BANDS
        sig_filt = filt(ZWTV_BAND_FILTERS[i], sig)
        # ISO preferred center frequency for this band (i is 1-based here;
        # MoSQITo's 0-based `i_bands = i - 1`, so `(i_bands - 16) == (i -
        # 17)`); coincides with this package's `BAND_CENTERS_HZ[i]`.
        center_freq = 10.0^((i - 17) / 10.0) * 1000.0
        sig_filt = _square_and_smooth(sig_filt, center_freq, fs)
        decimated = @view sig_filt[1:_ZWTV_DEC_FACTOR:end]
        @. band_levels[i, :] = 10.0 * log10((decimated + _ZWTV_TINY_VALUE) / _ZWTV_I_REF)
    end
    return band_levels, band_time_axis
end

"""
    _square_and_smooth(sig, center_freq, fs) -> Vector{Float64}

Square `sig` and pass it through three cascaded first-order low-pass
filters with a frequency-dependent time constant (`2/(3*center_freq)` up
to 1000 Hz, `2/3000` above), each with zero initial state.

Transcribed from `_square_and_smooth.py` @ d990c33f94f1 (Apache-2.0),
citing ISO 532-1:2017 section 6.3.
"""
function _square_and_smooth(sig::AbstractVector{Float64}, center_freq::Float64, fs::Float64)
    tau = center_freq <= 1000.0 ? 2.0 / (3.0 * center_freq) : 2.0 / (3.0 * 1000.0)
    a1 = exp(-1.0 / (fs * tau))
    b0 = 1.0 - a1
    y = sig .^ 2
    for _ in 1:3
        y = _zwtv_onepole_lowpass(y, b0, a1)
    end
    return y
end

# y[n] = b0*x[n] + a1*y[n-1], y[-1] = 0 -- the difference equation for
# scipy's `lfilter([b0], [1, -a1], x)` (zero initial conditions).
function _zwtv_onepole_lowpass(x::AbstractVector{Float64}, b0::Float64, a1::Float64)
    y = similar(x)
    acc = 0.0
    @inbounds for i in eachindex(x)
        acc = muladd(a1, acc, b0 * x[i])
        y[i] = acc
    end
    return y
end

"""
    loudness_zwtv(signal, fs; field_type=:free, pa_per_unit=1.0) -> ZwickerTimeVaryingResult

Compute ISO 532-1:2017 clause 6 ("method for time-varying sounds",
colloquially "Method 2" -- not `ISO 532-2`, the unrelated Moore-Glasberg
standard) time-varying loudness from a time-domain signal. Internally:
[`ZwickerLoudnessAudio._third_octave_levels`](@ref) (a faithful
transcription of MoSQITo's front-end stage, distinct from
[`band_levels`](@ref)'s ANSI Butterworth bank used for Method 1) →
`ZwickerLoudness.zwicker_loudness_time_varying`.

# Keyword arguments
- `field_type::Symbol`: `:free` (default) or `:diffuse` sound field.
- `pa_per_unit::Real`: scales `signal` to pascals (default `1.0` assumes
  the signal is already in Pa).

# Returns
`ZwickerTimeVaryingResult` (from ZwickerLoudness.jl): `loudness_over_time`
(N(t), sone, 2 ms axis), `N5`/`N10` (sone), `specific_loudness` (240 x
n_out, sone/Bark), `time_axis` (s).

# Throws
`ArgumentError` if `fs != 48000` -- this front-end stage is a faithful
transcription of MoSQITo's `_third_octave_levels`, which itself hard-
requires 48 kHz (no resampler reimplementation); resample the signal to
48 kHz first (e.g. `DSP.resample`) and retry.

# `time_axis` note
The kernel (`ZwickerLoudness.zwicker_loudness_time_varying`) never sees
the original signal or `fs`, so it synthesizes an idealized 2 ms axis
(`0.0, 0.002, 0.004, ...`). MoSQITo's own front end instead builds
`time_axis` as `linspace(0, length(signal)/fs, nblocks)`, i.e. it divides
the ORIGINAL signal duration into `nblocks - 1` equal intervals -- a step
about 0.5 % longer than the ideal 2 ms (`_third_octave_levels.py:259` @
d990c33f94f1). **This function reports MoSQITo's linspace-based axis**
(decimated by 4, matching MoSQITo's public `loudness_zwtv`'s `time_axis`
output), not the kernel's idealized one, so that `loudness_zwtv`'s output
matches MoSQITo's own end-to-end result field-for-field (crosscheck
parity was the deciding factor; see the package's fixture generation
scripts for the measured near-bit-exact match).

# Example
```julia
using ZwickerLoudnessAudio
fs = 48_000
t = range(0, 0.2, length=round(Int, 0.2 * fs))
sig_pa = 2e-5 * 10^(70/20) * sin.(2π * 1000 .* t) * sqrt(2)  # 1 kHz, 70 dB SPL
result = loudness_zwtv(sig_pa, fs)
result.loudness_over_time  # N(t) [sone]
result.N5                  # loudness exceeded 5% of the time [sone]
```
"""
function loudness_zwtv(signal::AbstractVector{<:Real}, fs::Real;
                       field_type::Symbol=:free, pa_per_unit::Real=1.0)
    isempty(signal) && throw(ArgumentError("signal is empty"))
    Float64(fs) == _ZWTV_FS || throw(ArgumentError(
        "loudness_zwtv requires fs == 48000 Hz (ISO 532-1:2017 clause 6 " *
        "front-end stage is a faithful transcription of MoSQITo's " *
        "_third_octave_levels, which itself hard-requires 48 kHz); " *
        "resample the signal to 48 kHz first (e.g. DSP.resample) and " *
        "retry. Got fs=$(fs)."))

    x = Float64.(signal) .* Float64(pa_per_unit)
    band_levels, band_time_axis = _third_octave_levels(x, _ZWTV_FS)
    kernel_result = zwicker_loudness_time_varying(band_levels; field_type=field_type)

    # See the `time_axis` note in this function's docstring: report
    # MoSQITo's own linspace-based axis, decimated by 4, in place of the
    # kernel's idealized one.
    frontend_time_axis = band_time_axis[1:_ZWTV_KERNEL_DEC_FACTOR:end]
    return ZwickerTimeVaryingResult(kernel_result.loudness_over_time, kernel_result.N5,
                                    kernel_result.N10, kernel_result.specific_loudness,
                                    frontend_time_axis)
end

"""
    loudness_zwtv(samples::AbstractMatrix, fs; channel=:mono, field_type=:free,
                 pa_per_unit=1.0) -> ZwickerTimeVaryingResult

Compute time-varying loudness from a samples × channels matrix. `channel=:mono`
(the default) averages all channels into one signal; an integer picks that
1-based channel. Mirrors [`loudness_zwst`](@ref)'s matrix entry.
"""
function loudness_zwtv(samples::AbstractMatrix{<:Real}, fs::Real;
                       channel::Union{Symbol,Integer}=:mono,
                       field_type::Symbol=:free, pa_per_unit::Real=1.0)
    signal = _select_channel(samples, channel)
    return loudness_zwtv(signal, Float64(fs); field_type=field_type, pa_per_unit=pa_per_unit)
end

"""
    loudness_zwtv(path::AbstractString; channel=:mono, field_type=:free,
                 pa_per_unit=1.0) -> ZwickerTimeVaryingResult

Read a `.wav` file and compute ISO 532-1:2017 clause 6 time-varying
loudness. For multichannel files, `channel=:mono` averages all channels;
an integer picks one 1-based channel. Mirrors [`loudness_zwst`](@ref)'s
wav entry. The file must be sampled at 48 kHz (see the vector method's
`ArgumentError` for the resample-first contract).
"""
function loudness_zwtv(path::AbstractString;
                       channel::Union{Symbol,Integer}=:mono,
                       field_type::Symbol=:free, pa_per_unit::Real=1.0)
    raw, fs = wavread(path)
    return loudness_zwtv(raw, Float64(fs);
                         channel=channel, field_type=field_type, pa_per_unit=pa_per_unit)
end
