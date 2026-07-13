# One-time generator for test/fixtures/zwtv_frontend_fixtures.jl.
# Synthesizes the SAME deterministic time-varying signals as
# ZwickerLoudness.jl's scripts/generate_mosqito_zwtv_crosscheck.jl (same
# tone/tone_burst/am_tone recipe and case parameters -- see that repo's
# .superpowers/sdd/zwtv-pins.md for the verified kernel/front-end
# boundary), hands them to MoSQITo (via uv, both the private front-end
# stage `_third_octave_levels` and the public `loudness_zwtv` API), and
# vendors the front-end-level fixtures Task 4 TDDs against: per-case
# (band_levels, band_time_axis, N_t, time_axis, N5/N10).
# Requires uv and network on first run.
# Usage: julia scripts/generate_mosqito_zwtv_frontend_crosscheck.jl
include(joinpath(@__DIR__, "..", "test", "support", "zwtv_generators.jl"))

# (name, fs, field_type, signal) -- same 3 cases/params as
# ZwickerLoudness.jl's kernel rig, so band_levels here is directly
# comparable to that repo's already-vendored ZWTV_KERNEL_FIXTURE_CASES.
# `test/support/zwtv_generators.jl` is also included directly by
# test/test_zwtv_frontend.jl, so both the fixture and the test suite
# regenerate byte-identical signals from one source.
cases = [(name, fs, string(field_type), sig) for (name, fs, field_type, sig) in zwtv_generator_cases()]

jsonvec(v) = string("[", join(string.(v), ","), "]")
tmp = tempname() * ".json"
open(tmp, "w") do io
    entries = String[]
    for (name, fs, field_type, sig) in cases
        push!(entries,
            """{"name": "$name", "fs": $fs, "field_type": "$field_type", "signal": $(jsonvec(sig))}""")
    end
    print(io, "[", join(entries, ","), "]")
end

out = joinpath(@__DIR__, "..", "test", "fixtures", "zwtv_frontend_fixtures.jl")

# Optional: one ISO 532-1 Annex B time-varying signal, end-to-end, DERIVED
# NUMBERS ONLY (see crosscheck_zwtv_frontend.py's Annex B block for the
# licensing reasoning, mirrored from ZwickerLoudness.jl). Points at a
# *local* MoSQITo checkout at the pin -- never committed, and the path is
# only used transiently at generation time. Same Annex B.4 signal as
# ZwickerLoudness.jl's kernel rig for cross-repo consistency.
annexb_wav = joinpath(
    "/tmp", "mosqito-pinned", "validations", "sq_metrics", "loudness_zwtv",
    "input", "ISO_532-1", "Annex B.4", "Test signal 6 (tone 250 Hz 30 dB - 80 dB).wav",
)
annexb_name = "annexb_test_signal_6_tone_250hz"
annexb_provenance = "MoSQITo validations/sq_metrics/loudness_zwtv/input/ISO_532-1/Annex B.4/" *
    "Test signal 6 (tone 250 Hz 30 dB - 80 dB).wav @ d990c33f94f1 (Apache-2.0); " *
    "underlying signal is ISO 532-1:2017 Annex B.4 Test signal 6, ISO-copyrighted."

if isfile(annexb_wav)
    run(`uv run --with mosqito==1.2.1 --with numpy --with matplotlib python $(joinpath(@__DIR__, "crosscheck_zwtv_frontend.py")) $tmp $out $annexb_wav $annexb_name $annexb_provenance`)
else
    @warn "Annex B wav not found locally at $annexb_wav; skipping Annex B derived-numbers fixture (clone MoSQITo @ d990c33f94f1 to /tmp/mosqito-pinned to include it)."
    run(`uv run --with mosqito==1.2.1 --with numpy --with matplotlib python $(joinpath(@__DIR__, "crosscheck_zwtv_frontend.py")) $tmp $out`)
end
