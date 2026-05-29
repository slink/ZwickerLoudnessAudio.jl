module ZwickerLoudnessAudio

using DSP
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

const _P_REF = 2e-5            # Reference acoustic pressure in pascals (20 µPa).
const _DEFAULT_ORDER = 3       # Filter order matches MoSQITo / ANSI S1.1-1986 default.
const _SILENCE_DB = -60.0      # Placeholder for bands above Nyquist or numerically silent.

# Memoize the per-band filterbank by (fs, order). digitalfilter+SOS conversion
# is the bulk of band_levels' work; repeated invocations on the same sample
# rate (e.g. processing many WAV files) reuse the cached designs. Entries are
# Vectors of length 28; bands above Nyquist hold `nothing`.
const _BAND_FILTER_CACHE = Dict{Tuple{Float64,Int},Vector{Any}}()
const _BAND_FILTER_CACHE_LOCK = ReentrantLock()

function _band_filters(fs::Real, order::Int)
    key = (Float64(fs), order)
    lock(_BAND_FILTER_CACHE_LOCK) do
        get!(_BAND_FILTER_CACHE, key) do
            α = _ansi_alpha(order)
            nyq = Float64(fs) / 2.0
            filters = Vector{Any}(undef, length(BAND_CENTERS_HZ))
            for (i, fc) in enumerate(BAND_CENTERS_HZ)
                w1 = fc / nyq / α
                w2 = fc / nyq * α
                if w1 <= 0.0 || w2 >= 1.0
                    filters[i] = nothing
                else
                    filters[i] = convert(SecondOrderSections,
                        digitalfilter(Bandpass(w1, w2), Butterworth(order)))
                end
            end
            return filters
        end
    end
end

_band_filter_cache_size() = length(_BAND_FILTER_CACHE)

function _clear_band_filter_cache!()
    lock(_BAND_FILTER_CACHE_LOCK) do
        empty!(_BAND_FILTER_CACHE)
    end
    return nothing
end

# ANSI S1.1-1986 design-bandwidth correction. The nominal third-octave edges
# `[fc / 2^(1/6), fc * 2^(1/6)]` give the *ideal* brick-wall band. For an
# Order-N Butterworth approximation, the design passband must be widened so
# that the noise-equivalent bandwidth equals the ideal third-octave bandwidth.
# Returns the multiplier `α` such that the design edges are `[fc/α, fc·α]`.
function _ansi_alpha(order::Int, bands_per_octave::Int=3)
    b = 1.0 / bands_per_octave
    # Reference (ideal-band) quotient: fc / (f2 - f1) for f1 = fc/2^(b/2), f2 = fc·2^(b/2).
    qr = 1.0 / (2.0^(b/2) - 2.0^(-b/2))
    # Design quotient from ANSI S1.1-1986 eq.9.
    qd = (π / 2 / order) / sin(π / 2 / order) * qr
    return (1.0 + sqrt(1.0 + 4.0 * qd^2)) / (2.0 * qd)
end

"""
    band_levels(signal, fs; pa_per_unit=1.0, order=3) -> Vector{Float64}

Compute one-third-octave band SPL [dB re 20 µPa] for `signal` sampled at `fs`
Hz. Returns 28 values for the standard band centers 25 Hz – 12.5 kHz.

Each band is filtered with a Butterworth bandpass designed per
ANSI S1.1-1986 (the same approach MoSQITo uses): the design edges
`[fc/α, fc·α]` are widened from the nominal third-octave edges by the
Order-N bandwidth correction `α`, which makes the noise-equivalent
bandwidth match the ideal band. Default order 3 gives Class-3 ANSI
behavior; bump to 6 or 8 for tighter selectivity.

`pa_per_unit` scales the raw signal to pascals (default `1.0` assumes Pa).
Bands whose upper edge exceeds Nyquist are clamped to a silence placeholder.

Filter designs are memoized by `(fs, order)`, so repeated calls at the same
sample rate (e.g. processing many WAV files) reuse the cached SOS objects.
"""
function band_levels(signal::AbstractVector{<:Real}, fs::Real;
                     pa_per_unit::Real=1.0, order::Int=_DEFAULT_ORDER)
    isempty(signal) && throw(ArgumentError("signal is empty"))
    fs > 0 || throw(ArgumentError("fs must be positive"))
    order >= 1 || throw(ArgumentError("order must be >= 1"))

    x = Float64.(signal) .* Float64(pa_per_unit)
    filters = _band_filters(fs, order)

    levels = fill(_SILENCE_DB, length(BAND_CENTERS_HZ))
    for (i, sos) in enumerate(filters)
        sos === nothing && continue
        y = filt(sos, x)
        rms = sqrt(mean(abs2, y))
        rms > 0 && (levels[i] = 20.0 * log10(rms / _P_REF))
    end
    return levels
end

"""
    loudness_zwst(signal, fs; field_type=:free, pa_per_unit=1.0, order=3) -> ZwickerResult

Compute ISO 532-1:2017 Method 1 stationary loudness from a time-domain signal.
Internally: [`band_levels`](@ref) → `ZwickerLoudness.zwicker_loudness`.
"""
function loudness_zwst(signal::AbstractVector{<:Real}, fs::Real;
                       field_type::Symbol=:free, pa_per_unit::Real=1.0,
                       order::Int=_DEFAULT_ORDER)
    levels = band_levels(signal, fs; pa_per_unit=pa_per_unit, order=order)
    return zwicker_loudness(levels; field_type=field_type)
end

"""
    loudness_zwst(samples::AbstractMatrix, fs; channel=:mono, field_type=:free,
                  pa_per_unit=1.0, order=3) -> ZwickerResult

Compute stationary loudness from a samples × channels matrix. `channel=:mono`
(the default) averages all channels into one signal; an integer picks that
1-based channel.
"""
function loudness_zwst(samples::AbstractMatrix{<:Real}, fs::Real;
                       channel::Union{Symbol,Integer}=:mono,
                       field_type::Symbol=:free, pa_per_unit::Real=1.0,
                       order::Int=_DEFAULT_ORDER)
    signal = _select_channel(samples, channel)
    return loudness_zwst(signal, Float64(fs);
                         field_type=field_type, pa_per_unit=pa_per_unit, order=order)
end

function _select_channel(samples::AbstractMatrix{<:Real},
                         channel::Union{Symbol,Integer})
    nch = size(samples, 2)
    if channel === :mono
        return nch == 1 ? vec(samples) : vec(mean(samples; dims=2))
    elseif channel isa Integer
        (1 <= channel <= nch) || throw(ArgumentError(
            "channel $channel out of range for $nch-channel input"))
        return Vector(samples[:, channel])
    else
        throw(ArgumentError("channel must be :mono or a 1-based Integer, got $channel"))
    end
end

"""
    loudness_zwst(path::AbstractString; channel=:mono, field_type=:free,
                  pa_per_unit=1.0, order=3) -> ZwickerResult

Read a `.wav` file and compute ISO 532-1:2017 Method 1 stationary loudness.
For multichannel files, `channel=:mono` averages all channels; an integer
picks one 1-based channel.
"""
function loudness_zwst(path::AbstractString;
                       channel::Union{Symbol,Integer}=:mono,
                       field_type::Symbol=:free, pa_per_unit::Real=1.0,
                       order::Int=_DEFAULT_ORDER)
    raw, fs = wavread(path)
    return loudness_zwst(raw, Float64(fs);
                         channel=channel, field_type=field_type,
                         pa_per_unit=pa_per_unit, order=order)
end

end # module
