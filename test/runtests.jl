using Test
using Statistics
using ZwickerLoudness: ZwickerResult
using ZwickerLoudnessAudio

# Generate a sinusoidal stimulus at frequency `f` Hz, sampled at `fs`, of
# duration `dur` seconds, scaled so its RMS pressure corresponds to `level_db`
# SPL re 20 µPa.
function make_tone(f::Real, fs::Real, dur::Real, level_db::Real)
    t = range(0, dur, length=round(Int, dur * fs))
    raw = sin.(2π * f .* t)
    rms_target = 2e-5 * 10.0^(level_db / 20.0)
    rms_now = sqrt(mean(raw .^ 2))
    return raw .* (rms_target / rms_now)
end

@testset "ZwickerLoudnessAudio" begin

    @testset "band_levels shape and basic correctness" begin
        fs = 48_000
        sig = make_tone(1000, fs, 0.5, 60.0)
        L = band_levels(sig, fs)
        @test length(L) == 28
        @test all(isfinite, L)
        # The 1 kHz band (index 17) should dominate.
        @test argmax(L) == 17
        # 1 kHz tone at 60 dB SPL: band-17 level should be in the 55-65 dB ballpark.
        @test 55 < L[17] < 65
    end

    @testset "loudness_zwst(signal, fs) returns a ZwickerResult" begin
        fs = 48_000
        sig = make_tone(1000, fs, 0.5, 60.0)
        r = loudness_zwst(sig, fs)
        @test r isa ZwickerResult
        @test isfinite(r.loudness)
        @test isfinite(r.loudness_level)
        @test length(r.specific_loudness) == 240
        # 1 kHz tone at 60 dB SPL: ISO Signal 3 reference is 4.019 sone.
        @test isapprox(r.loudness, 4.019; rtol=0.05)
    end

    @testset "broadband noise produces non-zero loudness" begin
        fs = 48_000
        n = round(Int, 0.5 * fs)
        # Pink-ish noise via cumulative sum (rough; smoke check only).
        raw = randn(n)
        rms_now = sqrt(mean(raw .^ 2))
        target = 2e-5 * 10.0^(60.0 / 20.0)
        sig = raw .* (target / rms_now)
        r = loudness_zwst(sig, fs)
        @test r.loudness > 0
    end

    @testset "field_type kwarg is honored" begin
        fs = 48_000
        sig = make_tone(1000, fs, 0.5, 60.0)
        rf = loudness_zwst(sig, fs; field_type=:free)
        rd = loudness_zwst(sig, fs; field_type=:diffuse)
        @test rf.loudness != rd.loudness
    end

    @testset "invalid input throws" begin
        @test_throws ArgumentError band_levels(Float64[], 48_000)
        @test_throws ArgumentError band_levels([1.0], 0.0)
    end

    # ============================================================ #
    #  Annex B Signals 2-4 conformance.
    #
    #  Reference values come from the .wav files in ISO 532-1:2017
    #  Annex B (also distributed by MoSQITo). All three pass within
    #  the ISO 5% tolerance with the order-3 Butterworth filterbank
    #  designed per ANSI S1.1-1986 (matching MoSQITo's noct_spectrum).
    #
    #  Signal 5 (pink noise) needs a calibrated pink-noise generator
    #  to test reproducibly; deferred.
    # ============================================================ #
    @testset "Annex B Signal 2: 250 Hz @ 80 dB SPL -> ~14.655 sone" begin
        fs = 48_000
        sig = make_tone(250, fs, 0.5, 80.0)
        r = loudness_zwst(sig, fs)
        @test isapprox(r.loudness, 14.655; rtol=0.05)
    end

    @testset "Annex B Signal 3: 1 kHz @ 60 dB SPL -> ~4.019 sone" begin
        fs = 48_000
        sig = make_tone(1000, fs, 0.5, 60.0)
        r = loudness_zwst(sig, fs)
        @test isapprox(r.loudness, 4.019; rtol=0.05)
    end

    @testset "Annex B Signal 4: 4 kHz @ 40 dB SPL -> ~1.549 sone" begin
        fs = 48_000
        sig = make_tone(4000, fs, 0.5, 40.0)
        r = loudness_zwst(sig, fs)
        @test isapprox(r.loudness, 1.549; rtol=0.05)
    end

    @testset "order kwarg is plumbed through" begin
        fs = 48_000
        sig = make_tone(1000, fs, 0.5, 60.0)
        r3 = loudness_zwst(sig, fs; order=3)
        r6 = loudness_zwst(sig, fs; order=6)
        # Higher-order filter has tighter selectivity -> less neighbor leakage
        # from a tone -> lower spreading contribution -> smaller total loudness.
        @test r6.loudness < r3.loudness
    end

end
