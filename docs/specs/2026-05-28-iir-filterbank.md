# IIR Filterbank — Annex B Conformance

**Date:** 2026-05-28
**Status:** Implemented
**Supersedes the smoke-filterbank caveats in:** [`2026-05-28-scaffold-design.md`](2026-05-28-scaffold-design.md)

## Goal

Replace the smoke FFT-power integrator in `band_levels` with a per-band IIR
filterbank that produces 28-band SPL values from a time signal accurately
enough to reproduce ISO 532-1:2017 Annex B Signals 2, 3, and 4 within the
standard's ±5% tolerance.

## Why the smoke filterbank wasn't enough

The smoke version summed FFT power inside brick-wall third-octave edges
`[fc/2^(1/6), fc·2^(1/6)]` for each band. For a stationary tone this
produces a single dominant band at the tone's center frequency and
near-silent neighbors — `(band_n=60 dB, others≈-∞)`. The Zwicker kernel
applied to that idealized distribution gives ~3.49 sone for a 1 kHz tone at
60 dB SPL, vs the ISO reference 4.019 sone (a −13% miss).

A real third-octave filter has finite skirt rolloff. For a 1 kHz tone, the
800 Hz and 1.25 kHz neighbor filters still pass roughly 41 dB of the tone
(at order 3) through their attenuated skirts. Those non-trivial neighbor
levels feed the Zwicker spreading function and lift total loudness up to
~4 sone, matching ISO.

## Design

Per band:
1. Compute design edges using the **ANSI S1.1-1986 Order-N bandwidth
   correction**: `[fc/α, fc·α]` where `α` is chosen so the
   noise-equivalent bandwidth of the Order-N Butterworth filter equals the
   ideal third-octave bandwidth. This is the same correction MoSQITo uses.
2. Normalize edges by Nyquist (`fs/2`).
3. Design a Butterworth bandpass via `DSP.digitalfilter(Bandpass(w1, w2),
   Butterworth(order))`.
4. Convert to Second-Order Sections (`SecondOrderSections`) for numerical
   stability at low normalized frequencies.
5. Filter the signal: `y = filt(sos, x)`.
6. Compute RMS, convert to dB SPL re 20 µPa.

Bands whose upper design edge exceeds Nyquist are clamped to the kernel's
silence placeholder.

Default `order = 3` matches MoSQITo and ANSI S1.1-1986's default Order-3 /
Type 0–C specification. An `order` kwarg is plumbed through `band_levels`
and `loudness_zwst`; higher orders give tighter selectivity at the cost of
filter-design conditioning for very low bands.

## Verification

| Signal | Reference | Computed | Δ |
|--------|-----------|----------|---|
| 2 (250 Hz @ 80 dB) | 14.655 sone | 15.087 sone | +2.9% |
| 3 (1 kHz @ 60 dB) | 4.019 sone | 4.162 sone | +3.6% |
| 4 (4 kHz @ 40 dB) | 1.549 sone | 1.557 sone | +0.5% |

All within the ISO 532-1 ±5% tolerance. The previously `@test_broken`
Annex B testsets are now plain `@test`.

## Not addressed in this change

- **Signal 5** (pink noise @ 60 dB SPL, reference 10.498 sone): needs a
  calibrated pink-noise generator (Voss-McCartney or a 1/f filter on white
  noise). The generator implementation, not the filterbank, is the gap.
- **Multichannel signals.** Still mixed to mono.
- **Sample-rate resampling.** Callers must provide `fs ≥ 48 kHz` (or at
  least high enough that 12.5 kHz · α < Nyquist). No automatic resampling.
- **Low-band decimation** (à la MoSQITo's `_n_oct_time_filter` when
  `fc < fs/200`). The current per-call SOS design appears stable across
  all 28 bands at the test sample rates, but very low rates may eventually
  warrant decimation for the lowest bands.
- **Performance.** Filters are designed on every call. A future optimization
  is to cache `(fs, order) → Vector{SOS}` so batch processing reuses them.

## Dependencies

- **Added:** `DSP.jl` (filter design + SOS filtering).
- **Removed:** direct `FFTW.jl` dep — no longer used in the source. DSP
  still pulls FFTW in transitively for unrelated functionality.
