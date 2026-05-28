# ZwickerLoudnessAudio.jl — Initial Scaffold Design

**Date:** 2026-05-28
**Status:** Approved (user-confirmed in brainstorming session)

## Purpose

Companion package to [ZwickerLoudness.jl](https://github.com/slink/ZwickerLoudness.jl).
Wraps the kernel with audio-file ingest and a third-octave filterbank so users
can compute Zwicker loudness directly from a `.wav` file or a time-domain
signal, instead of supplying pre-computed 28-band SPL values.

This satisfies the full ISO 532-1:2017 Method 1 pipeline shape:
`wav → third-octave bands → core+spreading → N (sone), L_N (phon), N'(z)`.

## Scope of this scaffold

**Goal:** end-to-end smoke. The pipeline runs from `.wav` to a `ZwickerResult`
and returns plausible numbers. **Numerical conformance against ISO 532-1
Annex B Signals 2–5 is explicitly deferred** to a follow-up that lands the
proper IEC 61260 IIR filterbank.

This is intentional. The filterbank is the hard part; shipping a smoke
pipeline first proves the wiring and gives a single place to swap the
filterbank later without touching the public API.

## Package identity

- **Name:** `ZwickerLoudnessAudio`
- **License:** MIT (match parent)
- **Initial version:** `0.1.0`
- **Repo:** local-only initially (no GitHub remote, no push) — matches solo
  workflow preference; remote can be added later.

## Dependencies

| Package | Purpose |
|---|---|
| `ZwickerLoudness` | Core kernel; declared via `Pkg.develop(path="../ZwickerLoudness.jl")` during dev so parent edits flow through immediately |
| `WAV` | Read `.wav` files |
| `FFTW` | FFT for band-energy estimation (smoke path) |
| `Statistics` (stdlib) | RMS / power calculations |

## Public API

```julia
loudness_zwst(signal::AbstractVector{<:Real}, fs::Real; field_type=:free, pa_per_unit=1.0) -> ZwickerResult
loudness_zwst(path::AbstractString; field_type=:free, pa_per_unit=1.0) -> ZwickerResult
band_levels(signal::AbstractVector{<:Real}, fs::Real; pa_per_unit=1.0) -> Vector{Float64}
```

- `signal` is a 1-D time series. Multichannel handling is out of scope; callers
  mix to mono themselves for now.
- `fs` is sample rate in Hz.
- `pa_per_unit` is the calibration factor mapping signal units to pascals. By
  convention, leave at `1.0` if the signal is already in Pa.
- `band_levels` returns 28 SPL values (dB re 20 µPa) for the standard
  one-third-octave centers from 25 Hz to 12.5 kHz.
- `loudness_zwst` returns the same `ZwickerResult` struct from `ZwickerLoudness`.

## Filterbank (smoke implementation)

1. Compute FFT of the full signal (no windowing for the smoke path; deliberately
   crude).
2. Convert to per-bin power.
3. For each of 28 third-octave bands (centers 25 Hz to 12.5 kHz per ANSI/IEC),
   sum power for all FFT bins whose frequency falls within the band edges
   `[fc / 2^(1/6), fc · 2^(1/6)]`.
4. Convert to RMS pressure, then to dB SPL (re 20 µPa).

Limitations of the smoke path (documented in code):
- No window function — leakage will redistribute energy across bands.
- No filter rolloff modeling — sharp brick-wall band edges instead of IEC 61260
  slope masks.
- Steady-state assumption: entire signal averaged as one segment.

These are exactly the things the IIR upgrade will fix.

## Tests (smoke level)

1. **End-to-end runs:** generate a 1 kHz tone at ~60 dB SPL, run through
   `loudness_zwst`, assert the result is finite and in a loose plausible range
   (e.g., 2–6 sone). Tight conformance values come later.
2. **band_levels shape:** returns 28 finite values.
3. **Broadband:** white noise at known RMS → loudness > 0.
4. **Annex B Signals 2–5 placeholder:** scaffolded as `@test_broken` with the
   target loudness values (14.655, 4.019, 1.549, 10.498) and a comment
   pointing at the IIR upgrade as the gating work.

## Infrastructure

Mirrors `ZwickerLoudness.jl`:
- `.github/workflows/CI.yml` — Julia 1.10 + 1, ubuntu/macos/windows
- `.github/dependabot.yml` — github-actions weekly
- `codecov.yml` + codecov-action step (token wiring deferred until a remote
  exists)
- `.gitignore` — Manifest.toml, coverage outputs, `.worktrees/`
- `README.md` — install, usage, **explicit scope/limitation note**

## Out of scope

Deferred to follow-up work (each is its own design conversation when we get
there):
- IEC 61260-compliant IIR filterbank
- Annex B Signals 2–5 numerical conformance
- Multichannel signal support
- Sample-rate resampling (will require `fs ≥ 48 kHz` initially)
- Calibration handling beyond the `pa_per_unit` scaling parameter
- Time-varying loudness (Method 2)
- GitHub remote + Codecov token wiring
