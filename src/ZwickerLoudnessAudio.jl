module ZwickerLoudnessAudio

using FFTW
using Statistics
using WAV
using ZwickerLoudness

export band_levels, loudness_zwst

# ANSI/IEC standard one-third-octave band center frequencies (Hz),
# 25 Hz to 12.5 kHz, 28 bands. Matches the input shape ZwickerLoudness expects.
const BAND_CENTERS_HZ = [
    25.0,    31.5,    40.0,    50.0,    63.0,    80.0,   100.0,
   125.0,   160.0,   200.0,   250.0,   315.0,   400.0,   500.0,
   630.0,   800.0,  1000.0,  1250.0,  1600.0,  2000.0,  2500.0,
  3150.0,  4000.0,  5000.0,  6300.0,  8000.0, 10000.0, 12500.0,
]

const _ONE_SIXTH_OCTAVE = 2.0^(1.0/6.0)
const _P_REF = 2e-5  # Reference acoustic pressure in pascals (20 µPa).

"""
    band_levels(signal, fs; pa_per_unit=1.0) -> Vector{Float64}

Compute one-third-octave band SPL [dB re 20 µPa] for `signal` sampled at `fs` Hz.
Returns 28 values for the standard band centers 25 Hz – 12.5 kHz.

`pa_per_unit` scales the raw signal to pascals (default `1.0` assumes the signal
is already in Pa).

# Implementation note (smoke)

This is the **smoke-level** filterbank: it FFTs the whole signal (no window,
no segmentation) and sums power for FFT bins falling inside each band's
brick-wall edges `[fc / 2^(1/6), fc · 2^(1/6)]`. It is not IEC 61260
compliant — neighboring-band rejection is overstated and tone leakage is
under-modeled. Replace with an IIR filterbank for ISO 532-1 Annex B
conformance against Signals 2-5.
"""
function band_levels(signal::AbstractVector{<:Real}, fs::Real; pa_per_unit::Real=1.0)
    n = length(signal)
    n == 0 && throw(ArgumentError("signal is empty"))
    fs > 0 || throw(ArgumentError("fs must be positive"))

    # Convert to pascals.
    x = Float64.(signal) .* Float64(pa_per_unit)

    # Real FFT → power per positive-frequency bin. Parseval-consistent scaling:
    # one-sided power spectrum of a real signal.
    X = rfft(x)
    nbins = length(X)
    # |X|^2 / n gives power per bin in (Pa^2 · samples). Divide by fs to get
    # power spectral density (Pa^2/Hz); equivalently sum |X[k]|^2 * 2 / n^2 to
    # get mean-square per bin in Pa^2 (one-sided, doubled except DC/Nyquist).
    psq = abs2.(X)
    # Convert to mean-square contribution per bin (Pa^2).
    # For RMS reconstruction over the full positive spectrum (one-sided):
    #   sum_k contribution_k == mean(x^2)
    # Use the standard one-sided normalization with DC + Nyquist not doubled.
    contrib = psq ./ (Float64(n)^2)
    contrib[2:end] .*= 2.0
    # Nyquist bin (for even n) shouldn't be doubled — undo:
    if iseven(n)
        contrib[end] /= 2.0
    end

    # Bin frequencies.
    freqs = (0:nbins-1) .* (Float64(fs) / Float64(n))

    # Integrate power across each band, then convert to dB SPL.
    levels = fill(-Inf, length(BAND_CENTERS_HZ))
    for (i, fc) in enumerate(BAND_CENTERS_HZ)
        fl = fc / _ONE_SIXTH_OCTAVE
        fu = fc * _ONE_SIXTH_OCTAVE
        # Avoid bands extending past Nyquist.
        fu > fs / 2 && (fu = fs / 2)
        fl >= fu && continue
        in_band = (freqs .>= fl) .& (freqs .< fu)
        p2 = sum(contrib[in_band])
        if p2 > 0
            levels[i] = 10.0 * log10(p2 / _P_REF^2)
        end
    end

    # Replace -Inf placeholders with -60 dB so the downstream kernel treats
    # them as effectively silent without producing NaN/Inf math.
    return [isfinite(L) ? L : -60.0 for L in levels]
end

"""
    loudness_zwst(signal, fs; field_type=:free, pa_per_unit=1.0) -> ZwickerResult

Compute ISO 532-1:2017 Method 1 stationary loudness from a time-domain signal.
Internally: [`band_levels`](@ref) → `ZwickerLoudness.zwicker_loudness`.

See `band_levels` for the smoke-level caveats on the filterbank.
"""
function loudness_zwst(signal::AbstractVector{<:Real}, fs::Real;
                       field_type::Symbol=:free, pa_per_unit::Real=1.0)
    levels = band_levels(signal, fs; pa_per_unit=pa_per_unit)
    return zwicker_loudness(levels; field_type=field_type)
end

"""
    loudness_zwst(path::AbstractString; field_type=:free, pa_per_unit=1.0) -> ZwickerResult

Read a `.wav` file and compute ISO 532-1:2017 Method 1 stationary loudness.
Multichannel files are mixed to mono by averaging channels.
"""
function loudness_zwst(path::AbstractString;
                       field_type::Symbol=:free, pa_per_unit::Real=1.0)
    raw, fs = wavread(path)
    # Mix to mono if multichannel.
    signal = size(raw, 2) == 1 ? vec(raw) : vec(mean(raw; dims=2))
    return loudness_zwst(signal, Float64(fs); field_type=field_type, pa_per_unit=pa_per_unit)
end

end # module
