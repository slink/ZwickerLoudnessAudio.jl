using Test
using FFTW
using Logging
using Random
using Statistics
using WAV: WAVE_FORMAT_IEEE_FLOAT, wavwrite
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

# Generate deterministic pink noise by shaping Gaussian white noise in the
# frequency domain to an exact 1/sqrt(f) magnitude spectrum (i.e. 1/f power
# spectral density). Output is scaled so its overall RMS pressure corresponds
# to `level_db` SPL re 20 µPa.
#
# Why this approach over a Paul Kellet 6-pole IIR pink generator: Kellet's
# approximation overweights the 25 Hz – 12.5 kHz band region (because its
# lowest pole sits at ~9 Hz, so less of the signal's total energy escapes
# below the Zwicker band range), and consistently overestimates Signal 5
# loudness by ~7% across all seeds. The FFT-shaped generator centers on the
# ISO 10.498 sone reference within ±5% for typical seeds.
function make_pink_noise(fs::Real, dur::Real, level_db::Real; seed::Integer=42)
    rng = MersenneTwister(seed)
    n = round(Int, dur * fs)
    X = rfft(randn(rng, n))
    df = fs / n
    @inbounds for k in 2:length(X)
        X[k] /= sqrt((k - 1) * df)
    end
    X[1] = 0     # zero DC
    pink = irfft(X, n)
    rms_target = 2e-5 * 10.0^(level_db / 20.0)
    return pink .* (rms_target / sqrt(mean(pink .^ 2)))
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

    @testset "Annex B Signal 5: pink noise @ 60 dB SPL -> ~10.498 sone" begin
        fs = 48_000
        sig = make_pink_noise(fs, 4.0, 60.0; seed=20260528)
        r = loudness_zwst(sig, fs)
        @test isapprox(r.loudness, 10.498; rtol=0.05)
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

    @testset "filter cache: design reused across calls" begin
        ZwickerLoudnessAudio._clear_band_filter_cache!()
        @test ZwickerLoudnessAudio._band_filter_cache_size() == 0

        sig = make_tone(1000, 48_000, 0.1, 60.0)

        # First call populates one cache entry.
        band_levels(sig, 48_000)
        @test ZwickerLoudnessAudio._band_filter_cache_size() == 1

        # Repeat call with the same (fs, order) does not add a new entry.
        band_levels(sig, 48_000)
        @test ZwickerLoudnessAudio._band_filter_cache_size() == 1

        # Different order keys a fresh entry.
        band_levels(sig, 48_000; order=6)
        @test ZwickerLoudnessAudio._band_filter_cache_size() == 2

        # Different fs keys a fresh entry.
        sig2 = make_tone(1000, 44_100, 0.1, 60.0)
        band_levels(sig2, 44_100)
        @test ZwickerLoudnessAudio._band_filter_cache_size() == 3

        # Different fs types collapse on the same Float64 key.
        band_levels(sig2, 44_100.0)
        @test ZwickerLoudnessAudio._band_filter_cache_size() == 3
    end

    @testset "filter cache: cached results match uncached" begin
        sig = make_tone(1000, 48_000, 0.2, 60.0)
        ZwickerLoudnessAudio._clear_band_filter_cache!()
        L_first = band_levels(sig, 48_000)
        L_cached = band_levels(sig, 48_000)
        @test L_first == L_cached
    end

    @testset "multichannel matrix input" begin
        fs = 48_000
        L = make_tone(1000, fs, 0.5, 60.0)
        R = make_tone(4000, fs, 0.5, 40.0)
        samples = hcat(L, R)

        # :mono default mixes both channels into a mono signal.
        r_mono = loudness_zwst(samples, fs)
        @test r_mono isa ZwickerResult

        # channel=1 picks the 1 kHz channel only (≈ ISO Signal 3).
        r1 = loudness_zwst(samples, fs; channel=1)
        @test isapprox(r1.loudness, 4.019; rtol=0.05)

        # channel=2 picks the 4 kHz channel only (≈ ISO Signal 4).
        r2 = loudness_zwst(samples, fs; channel=2)
        @test isapprox(r2.loudness, 1.549; rtol=0.05)

        # Out-of-range channels are rejected.
        @test_throws ArgumentError loudness_zwst(samples, fs; channel=3)
        @test_throws ArgumentError loudness_zwst(samples, fs; channel=0)

        # Unknown Symbols (anything other than :mono) are rejected.
        @test_throws ArgumentError loudness_zwst(samples, fs; channel=:left)
    end

    @testset "bands above input Nyquist resolve to nothing in the cache" begin
        # At fs=10 kHz, the top several bands' design edges exceed Nyquist
        # and `_band_filters` stores `nothing` for them. The public API
        # auto-resamples before this code runs, but the internal contract
        # still applies if callers invoke `_band_filters` directly.
        filters = ZwickerLoudnessAudio._band_filters(10_000.0, 3)
        @test filters[end] === nothing
        @test any(x -> x === nothing, filters)
        @test any(x -> x !== nothing, filters)
    end

    @testset "high sample rates stay numerically stable" begin
        # At 192 kHz the 25 Hz band's normalized edge is ~1.3e-4 and at 384 kHz
        # it's ~6.5e-5 — small enough that direct-form filter coefficients
        # would lose precision. The SOS cascade representation keeps the
        # filterbank well-conditioned, and Annex B Signals 2/3/4 round-trip
        # within tolerance at all these rates. Lock that in.
        for fs in (96_000, 192_000, 384_000)
            r2 = loudness_zwst(make_tone(250, fs, 0.5, 80.0), fs)
            @test isapprox(r2.loudness, 14.655; rtol=0.05)
            r3 = loudness_zwst(make_tone(1000, fs, 0.5, 60.0), fs)
            @test isapprox(r3.loudness, 4.019; rtol=0.05)
            r4 = loudness_zwst(make_tone(4000, fs, 0.5, 40.0), fs)
            @test isapprox(r4.loudness, 1.549; rtol=0.05)
        end
    end

    @testset "low fs is auto-resampled" begin
        # A 16 kHz signal can't carry the 12.5 kHz band's design edge (~14 kHz
        # for order=3), so band_levels at 16 kHz would silently drop the top
        # band. With auto-resampling, the result should recover the ISO
        # Signal 3 reference within tolerance after upsampling to 48 kHz.
        sig = make_tone(1000, 16_000, 0.5, 60.0)
        r = @test_logs (:warn, r"resampling"i) loudness_zwst(sig, 16_000)
        @test isapprox(r.loudness, 4.019; rtol=0.05)
    end

    @testset "supported fs passes through without resampling" begin
        # Common CD/DVD rates (≥ 32 kHz) leave the input untouched.
        sig44 = make_tone(1000, 44_100, 0.5, 60.0)
        r = @test_logs min_level=Logging.Warn loudness_zwst(sig44, 44_100)
        @test isapprox(r.loudness, 4.019; rtol=0.05)
        sig48 = make_tone(1000, 48_000, 0.5, 60.0)
        r2 = @test_logs min_level=Logging.Warn loudness_zwst(sig48, 48_000)
        @test isapprox(r2.loudness, 4.019; rtol=0.05)
    end

    @testset "multichannel WAV file via channel kwarg" begin
        fs = 48_000
        L = make_tone(1000, fs, 0.5, 60.0)
        R = make_tone(4000, fs, 0.5, 40.0)
        path = tempname() * ".wav"
        try
            wavwrite(hcat(L, R), path; Fs=fs, nbits=32, compression=WAVE_FORMAT_IEEE_FLOAT)
            r_mono = loudness_zwst(path)
            @test r_mono isa ZwickerResult
            r1 = loudness_zwst(path; channel=1)
            @test isapprox(r1.loudness, 4.019; rtol=0.05)
            r2 = loudness_zwst(path; channel=2)
            @test isapprox(r2.loudness, 1.549; rtol=0.05)
        finally
            isfile(path) && rm(path)
        end
    end

end
