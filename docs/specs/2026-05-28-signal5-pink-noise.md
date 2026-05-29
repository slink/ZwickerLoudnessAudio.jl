# Annex B Signal 5 (Pink Noise) Conformance

**Date:** 2026-05-28
**Status:** Implemented
**Builds on:** [`2026-05-28-iir-filterbank.md`](2026-05-28-iir-filterbank.md)

## Goal

Add a deterministic pink-noise test stimulus and an Annex B Signal 5
conformance testset (10.498 sone @ 60 dB SPL, free field) — the only
remaining Annex B signal that wasn't covered when the IIR filterbank
landed.

## Generator design

A test helper `make_pink_noise(fs, dur, level_db; seed)` in `test/runtests.jl`:

1. Draw `n = round(dur*fs)` samples of Gaussian white noise from a seeded
   `MersenneTwister`.
2. Take the real FFT, scale bin `k` (k≥1) by `1/sqrt(k·fs/n)`, zero DC.
   This is an *exact* `1/sqrt(f)` magnitude shaping → `1/f` power spectral
   density → pink noise.
3. Inverse-FFT and rescale so overall RMS pressure corresponds to
   `level_db` SPL re 20 µPa.

## Why not Paul Kellet's 6-pole IIR

The Kellet refined pink filter is pure-Julia and seemingly the obvious
zero-dep choice. It produces a passable `1/f` spectrum (±0.05 dB from
~9 Hz to ~22 kHz), but on Signal 5 it consistently overshoots:

| Generator | 8-seed mean | 8-seed stddev |
|-----------|-------------|---------------|
| Kellet IIR, 4 s | 11.286 sone (+7.5%) | 0.096 sone (0.85%) |
| FFT 1/√f, 4 s | 10.531 sone (+0.3%) | 0.257 sone (2.4%) |

The Kellet miss is systematic — every seed and duration we tried lands
between +5.9% and +9.0%. The likely cause is that Kellet's lowest pole
(~9 Hz) leaves less of the total energy below 25 Hz than ideal `1/f`
would, so a larger fraction of the signal's energy lands inside the 28
Zwicker bands. Per-band SPL averages ~44 dB with Kellet vs ~42.4 dB with
the FFT generator at the same 60 dB overall SPL.

The FFT generator centers on the ISO reference and the chosen
`seed=20260528, dur=4.0s` lands at +2.6%, comfortably inside the ISO
±5% tolerance.

## Cost

`FFTW` is added to `[extras]` (test-only). It's already pulled in as a
transitive dep of `DSP`, so no new artifact downloads. The package's
runtime dependency set is unchanged.

## Verification

All 18 testsets pass on Julia 1.12 (macOS arm64). The new testset:

```
Annex B Signal 5: pink noise @ 60 dB SPL -> ~10.498 sone
```

asserts `isapprox(r.loudness, 10.498; rtol=0.05)`.

## Remaining gaps

None from the original IIR filterbank spec — Signal 5 was the last open
item. Filter caching, multichannel input handling, low-fs auto-resampling,
and high-fs numerical stability all landed in follow-up commits on `main`.
