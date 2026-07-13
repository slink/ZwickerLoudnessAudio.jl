# Deterministic time-varying test signal generators, shared between
# scripts/generate_mosqito_zwtv_frontend_crosscheck.jl (fixture generation)
# and test/test_zwtv_frontend.jl (regenerating the same signals to feed
# this package's own `loudness_zwtv`/`_third_octave_levels`). Identical to
# ZwickerLoudness.jl's kernel-repo rig generator (confirmed byte-identical
# via matching band_levels_sha256 across both repos' fixtures) -- kept in
# one place here per the PsychoacousticMetrics.jl `test/support/` convention.
using Statistics: std

# Calibrated pure tone, `duration` seconds at `fc` Hz, `spl` dB re 20 uPa.
# Amplitude calibration transcribed from MoSQITo utils/am_sine_generator.py
# @ d990c33f94f1 (std-normalization, not RMS -- ddof=0).
function zwtv_tone(fs::Real, fc::Real, duration::Real, spl::Real)
    n = round(Int, duration * fs)
    t = (0:n-1) ./ fs
    y = sin.(2π * fc .* t)
    A = 2e-5 * 10.0^(spl / 20)
    return y .* (A / std(y; corrected=false))
end

# Tone burst: `burst_duration` s of a `spl`-calibrated tone (level defined
# over the burst itself, not the padded silence) starting at `burst_start`
# s within a `total_duration` s silent buffer.
function zwtv_tone_burst(fs::Real, fc::Real, total_duration::Real, burst_duration::Real,
                          burst_start::Real, spl::Real)
    n_total = round(Int, total_duration * fs)
    n_burst = round(Int, burst_duration * fs)
    i0 = round(Int, burst_start * fs)
    burst = zwtv_tone(fs, fc, burst_duration, spl)
    sig = zeros(n_total)
    sig[i0+1:i0+n_burst] .= burst
    return sig
end

# Amplitude-modulated tone: (1 + sin(2*pi*fm*t)) * sin(2*pi*fc*t), calibrated
# to `spl` dB by std-normalization (same convention as `zwtv_tone`).
function zwtv_am_tone(fs::Real, fc::Real, fm::Real, duration::Real, spl::Real)
    n = round(Int, duration * fs)
    t = (0:n-1) ./ fs
    xmod = sin.(2π * fm .* t)
    xc = sin.(2π * fc .* t)
    y = (1 .+ xmod) .* xc
    A = 2e-5 * 10.0^(spl / 20)
    return y .* (A / std(y; corrected=false))
end

# (name, fs, field_type, signal) -- the SAME 3 cases/params as
# ZwickerLoudness.jl's kernel rig, so band_levels computed on them is
# directly comparable to that repo's already-vendored kernel fixture.
zwtv_generator_cases() = [
    ("steady_tone_1k_60db", 48000.0, :free,
        zwtv_tone(48000, 1000.0, 0.10, 60.0)),
    ("tone_burst_1k_70db_30ms", 48000.0, :free,
        zwtv_tone_burst(48000, 1000.0, 0.15, 0.03, 0.06, 70.0)),
    ("am_tone_1k_10hz_65db", 48000.0, :free,
        zwtv_am_tone(48000, 1000.0, 10.0, 0.20, 65.0)),
]
