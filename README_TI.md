# Transversely isotropic (TI) CPS wrappers

The `*_TI` functions in `functions/` are transversely isotropic (radially
anisotropic) versions of the isotropic CPS wrappers. Each one keeps the
calling sequence of its isotropic counterpart. Instead of surf96 / srfker96 /
sprep96 / sdisp96 / sregn96 / slegn96, they call the CPS TI programs
`tprep96`, `tdisp96`, `tregn96` and `tlegn96`, all of which are already in
`bin_v3.30/`.

| Isotropic | TI | Notes |
|---|---|---|
| `dispR_surf96` | `dispR_surf96_TI` | same arguments (`vec_T, model, nmode, disp_type, fref`) |
| `dispL_surf96` | `dispL_surf96_TI` | same arguments, plus optional `fref` |
| `calc_kernel96` | `calc_kernel96_TI` | returns `[kernel, dispersion]` structures (like `calc_anelastic_kernel96`) |
| `calc_anelastic_kernel96` | `calc_anelastic_kernel96_TI` | same fields as the isotropic version, plus TI kernels |
| `calc_eigenfunctions96` | `calc_eigenfunctions96_TI` | same fields, plus phv/grv/gamma and phase-velocity kernels |
| `calc_amplification_ellipticity_gamma` | `calc_amplification_ellipticity_gamma_TI` | same fields |
| `calc_admittance96_Ruan14` | `calc_admittance96_Ruan14_TI` | same fields |
| `writemod_surf96` | `writemod_surf96_TI` | writes a `TRANSVERSE ISOTROPIC` MODEL96 file |
| `readmod_surf96` | `readmod_surf96_TI` | also reads isotropic files (converted to TI) |
| `read_eig96` | `read_der96_TI` | reads the binary `t?egn96.der` file (full precision) |

New helpers:
- `convert_iso2TI`: isotropic model to TI model. Pass `xi`, `phi` and `eta` to add anisotropy while keeping the Voigt-average Vp and Vs.
- `convert_TI2iso`: TI model to its Voigt-average isotropic model (also returns xi, phi, eta).
- `calc_der96_TI`: the engine shared by the wrappers. It runs the CPS TI programs and rejects spurious roots (see below).
- `check_halfspace_eta_TI`, `cleanup_cps96_TI`: a warning check and temporary-file cleanup.

## TI model format

```
model_TI = [H Vpv Vsv rho Vph Vsh eta]            % elastic
model_TI = [H Vpv Vsv rho Vph Vsh eta Qp Qs]      % with Q
```

Columns 1–4 mean the same as in the isotropic `[H Vp Vs rho (Qp Qs)]` model.
`eta = F/(A-2L)`, which is 1 for an isotropic layer. To convert an isotropic
model to an equivalent TI model:

```matlab
model_TI = convert_iso2TI(model);                    % xi = phi = eta = 1
model_TI = convert_iso2TI(model, 1.1, 0.95, 0.92);   % add radial anisotropy
```

CPS stores `VPF = sqrt(F/rho)` in the model file, so layers with
`eta*(Vph^2 - 2*Vsv^2) < 0` cannot be written. This happens for Vp/Vs < sqrt(2)
when eta = 1.

## Kernels

`calc_kernel96_TI` returns layer kernels with respect to Vpv, Vph, Vsv, Vsh,
eta and rho, each with the others held fixed (the CPS parameterization). It
also returns the layer-thickness kernels (`dcdh`, `dudh`) and the
isotropic-equivalent kernels `dcdvs = dcdvsv + dcdvsh` and
`dcdvp = dcdvpv + dcdvph` (and the dU equivalents). The isotropic-equivalent
kernels can be compared directly with `calc_kernel96`.

CPS has no TI version of srfker96, so the code computes:
- **Group-velocity kernels** from the phase-velocity kernels (Rodi et al.,
  1975): `dU/dm = (U/c)(2-U/c) dc/dm + (U^2 w/c^2) d(dc/dm)/dw`, with a
  centered finite difference in frequency.
- **Anelastic kernels** (`calc_anelastic_kernel96_TI`) by applying the
  srfker96 causal-Q correction (`getgam` and `cslmod` in `srfker96.f`) to the
  elastic TI kernels.

> **Changes to the isotropic `calc_kernel96` and `calc_anelastic_kernel96`
> (November 2026).**
> - Outputs 6–7 (`dcdrho`, `dudrho`) used to be columns 7–8 of the srfker96
>   output, which are dc/dh and dU/dh (layer thickness), not density. They
>   are now real density kernels, from `sregn96`/`slegn96 -DER`
>   (`calc_density_kernel96`). srfker96's thickness kernels are returned as
>   new outputs 8–9 (`dcdh`, `dudh`), followed by the elastic phase and group
>   velocities.
> - srfker96's phase-velocity kernels dc/dVs and dc/dVp are now multiplied
>   by the earth-flattening factor. Because `calc_kernel96` requests only
>   group velocity, `srfdrr96`/`srfdrl96` map only the group-velocity
>   kernels back to the spherical model. Without the factor, the kernels
>   were too small by up to ~2.4% at 150 km depth.
> - `calc_anelastic_kernel96` now computes γ, the causal kernels and the Q
>   kernels from these corrected elastic kernels, using the srfker96
>   formulas. It uses `fref` and supports Love waves.
>
> The TI functions use the same conventions.

## Inversion for radial anisotropy

`test_kernel_inversion_vsvh_vpvh_rho_eta.m` is the TI version of
`test_kernel_inversion_vs_vp_rho.m`. It is a synthetic test that uses
`run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta.m` to invert Rayleigh and/or Love
phase velocities (`data_types = 'RL'`, `'R'` or `'L'`) for Vsv, Vsh, Vpv,
Vph, ρ and η. The setup follows the MINEOS Rayleigh-wave inversion:
- damping and smoothing of Vsv and Vsh
- Vpv/Vsv, Vph/Vsh and ρ/Vsv tied to the starting model
- ξ and η damped toward the starting model (`eps_xi`, `eps_eta`)
- ξ held fixed where Vsv < `vsv_thresh_xi` (crust and sediments)
- water layer fixed

With a single wave type, ξ is not constrained by the data, so use a large
`eps_xi`.

In the synthetic test, the true anisotropy (ξ = 1.08 at 30–150 km, with Vph
raised alongside Vsh, so φ = 1/ξ) is recovered well with 10–40 s R + L data
(ξ ≈ 1.06–1.09 at 50–120 km) provided that:
- **the norm damping is weak** (`eps_H` = 0.005). Stronger damping toward
  the starting model leaves the long-period Love data underpredicted
  (about +1.2σ) and ξ underestimated (≈ 1.04).
- **the smoothing is broken at the model's discontinuities** (water bottom,
  sediments, Moho). Without the breaks, the smoothing across the Moho
  prevents fitting the short-period data, and weak damping then puts the
  misfit into an oscillating ξ (1.16 at 80 km, 0.97 at 50 km).

> **Fix (November 2026): `break_constraint.m`.** For the second-derivative
> smoothing, it removed three equations at each break, including one that
> lies entirely below the break. The first layer below each break (e.g. the
> uppermost mantle) then had no smoothing and could take a large spurious
> velocity jump. Now only the two equations that cross the break are
> removed. This affects all inversions that break the smoothing
> (isotropic and TI).

### Oceanic LVZ tests (10–150 s)

- `test_kernel_inversion_vsvh_vpvh_rho_eta_LVZ.m`: true model = the
  starting card plus a fast lid, a strong Gaussian LVZ (Vsv ≈ 4.0 km/s at
  ~135 km) and ξ = 1.10 tapering to 1 at 220 km.
- `test_kernel_inversion_vsvh_vpvh_rho_eta_PA5.m`: true model = PA5
  (Gaherty et al., 1999; `MINEOS_CARD/pa5_5km.card`), starting model with
  a smooth isotropic mantle (no lid, LVZ or 68/166 km discontinuities).

Each script plots the R and L kernels in the starting and true models and
compares three inversions: kernels recalculated every iteration, kernels
fixed at the starting model, and a control with the same ξ but no LVZ.
Results (noise-free data, 1% uncertainties):
- **Synthetic LVZ:** with updated kernels, ξ is recovered to within 0.01 at
  22–260 km, as well as without the LVZ. Fixed kernels move the anisotropy
  too deep (ξ too low in the lid, +0.02 to +0.03 at 140–220 km). The LVZ
  mainly changes the Rayleigh sensitivity to ξ, which grows and moves into
  the LVZ.
- **PA5:** the thin crust and the fast lid make the 10–30 s mode-0 Love
  wave an LVZ-trapped mode (the crust/lid Love wave is on the mode-1
  branch). The mode-0 Love data then constrain the lid poorly. ξ in the lid
  is not recovered (1.01 vs 1.08 at 12–68 km) and the anisotropy is moved
  below 166 km (1.05 vs 1.00), whether or not the kernels are updated.
  The no-LVZ control recovers the lid ξ (1.06).
- **PA5, two stages** (`test_kernel_inversion_vsvh_vpvh_rho_eta_PA5_twostage.m`):
  first invert the Rayleigh data alone for an isotropic model (ξ held at
  1), which recovers the LVZ in Vsv, then use it as the starting model for
  the R + L inversion. The lid ξ improves from 1.01 to 1.05 (true 1.08), and
  the spurious ξ below 166 km drops from 1.05 to 1.02. This is nearly as
  good as starting from the true Vsv (lid ξ 1.06). The remaining error is
  the vertical smoothing of a sharp ξ profile.
- **Which Love branch is observed in PA5?** For a source and a receiver near
  the surface, the surface amplitude of the mode-1 Love wave (crust and lid)
  exceeds that of mode 0 (LVZ-trapped) by ×800 at 10 s, ×17 at 20 s and ×2
  at 40 s; they are comparable at 50–60 s and mode 0 dominates above
  ~70 s. Their group velocities are similar (4.45–4.56 km/s at 15–60 s), so
  they cannot be separated by time windowing. In such a model, a measured
  "fundamental-mode" Love wave below ~40 s is the mode-1 branch of the
  CPS numbering. The forward problem must follow the branch that is
  measured.

Extending the periods to 100 s improves ξ below ~100 km.

**Speed.** Each iteration computes the kernels and phase velocities
together in one CPS calculation. The kernels are phase-velocity only, so
two frequencies per period are computed instead of three, and only the
modes needed are computed (`nmodes`, default `nmode+1`).
`nit_recalc_G` (default `nit_recalc_c`) sets how often the kernels are
recalculated, independently of the phase velocities (`nit_recalc_c`). When
only the phase velocities are recalculated, the linearized prediction is
used to pick the root, so one frequency per period is computed. For the
synthetic test this cut the run time by about 3.6× with bit-identical
results. The final figure compares the true, starting and recovered Vsv,
Vsh, Vpv, Vph, ρ, η and ξ with formal 1σ uncertainties.
`card2mod_TI.m` converts a MINEOS card to the TI layered model format.

## Synthetic seismograms (mode summation)

`calc_synth96_TI.m` computes Z/R/T synthetic seismograms for a TI model by
surface-wave mode summation (CPS `tprep96`, `tdisp96`, `tlegn96`/`tregn96`,
`tpulse96`, `fmech96`, `f96tosac`):

```matlab
par.hs = 15; par.hr = 5; par.dt = 1; par.npts = 4096; par.fmax = 0.2; par.nmod = 20;
par.modesel = {'all', 0, 1};      % all modes, mode 0 only, mode 1 only
par.pulse = 't'; par.L = 3;       % triangle of base 2*L*dt
src.strike = 30; src.dip = 60; src.rake = 45; src.Mw = 6; src.az = 60;
syn = calc_synth96_TI(model, [1000 2000], par, src);   % syn.Z, syn.R, syn.T [npts x ndist x nsel]
```

The source is a double couple (strike/dip/rake with Mw or M0), a moment
tensor (`src.mt`) or a point force (`src.force`, e.g. to model noise
cross-correlations). A receiver exactly at the sea floor is moved 0.01 km
into the solid; at the water bottom, CPS gives zero Love-wave motion. Notes:
- **Frequency content:** the highest frequency is 1/(2·dt), or `par.fmax`
  if lower. The pulse, `par.nmod` (complete to ~5 s for 20 Love modes in
  PA5), Q and the layering also limit it. The record repeats every
  npts·dt seconds, so make it longer than the slowest arrival.
- **Mode numbers** are those of `tdisp96`, ordered by phase velocity.
  `tpulse96` does not reject the occasional spurious root (see
  `calc_der96_TI.m`). Body waves are not included.
- **Programs to build:** `tpulse96` (VOLIII), `fmech96` and `f96tosac`
  (VOLV) are **not in `bin_v3.30/`**. Build them from `src.zip`
  (PROGRAMS.330).

`test_synth_seismogram_TI.m` compares the all-modes, mode-0 and mode-1
synthetics for PA5 (Mw 6, 10 km below the sea floor, sea-floor station at
1000 and 2000 km). At 2000 km, mode 1 dominates the Love wave at 10–30 s
and mode 0 is small. At 30–60 s both are large and interfere, and mode 0
dominates above 60 s. On the vertical component, mode 1 also dominates at
10–15 s; mode 0 there is the slow wave guided by the water and sediment.
Mode 0 dominates the vertical component above 15 s.

## Benchmarks

```matlab
benchmark_TI_run_all   % or run benchmark_TI_dispersion, _kernels, _eigenfunctions, _anelastic_kernels
```

Each script compares the isotropic wrapper with the TI wrapper for a TI
model converted from an isotropic oceanic model (water, sediment, crust and
mantle). The TI kernels are also checked against finite differences for a
model with real anisotropy. A PASS/FAIL line is printed for each comparison.
The summary and the figures are written to `benchmark_TI_figs/`. Each
script plots the layered model(s) it uses (`plot_model_TI`) and keeps its
figures in its own numbering range, so `benchmark_TI_run_all` leaves every
figure open: dispersion 101–105, kernels 201–204, eigenfunctions 301–304,
anelastic kernels 401–403. To use CPS
binaries other than `bin_v3.30/`, set the environment variable `CPS_BIN`.

### Expected iso vs TI differences (spherical-earth models)

`writemod_surf96` and `writemod_surf96_TI` both write `SPHERICAL EARTH`
models. The CPS TI programs flatten the earth slightly differently from the
isotropic programs:
- **Earth radius.** The TI programs use 6371 km; the isotropic programs use 6370 km.
- **The modulus F.** For Rayleigh waves, `tdisp96`/`tregn96` do not flatten
  F, so a flattened isotropic layer ends up with eta slightly different from 1.

As a result, phase and group velocities differ by about 1e-4 (Rayleigh) and
1e-5 (Love) in relative terms. With `FLAT EARTH` models they agree to the
output precision.

Because F is not flattened, the Rayleigh-wave Vp-type kernels (dc/dVp,
dU/dVp and the Qp kernels) differ the most: up to ~1.5% at 50 s for the
benchmark model. The TI kernels are still consistent with TI dispersion
(finite-difference checks pass).

The kernels agree to about 1e-4 to 1e-3, which is the precision of the
srfker96 text output. That holds now that `calc_kernel96` applies the
spherical-earth mapping to srfker96's phase-velocity kernels (see above);
before the fix they differed by a factor of about `a/(a-z)`. The TI and
corrected isotropic kernels both agree with finite differences of the
spherical-earth phase velocity.

## Known CPS TI issues (handled or flagged by the wrappers)

1. **Root finding in `tdisp96` is less robust than in `sdisp96`.** For models
   with a water layer over slow sediments, `tdisp96` sometimes:
   - reports roots that are not modes of the model (nearly non-dispersive,
     near the slowest S velocity),
   - reports the same root twice, or
   - misses a root.

   Any of these renumbers the modes above it. The isotropic `sdisp96` and
   `surf96` did not do this for the same models. `calc_der96_TI`, which all
   the TI wrappers use, handles it in two ways:
   - **One period per run.** When given a list of periods, `tdisp96` caps the
     number of roots it accepts at each period at the number found at the
     previous (shorter) period. One spurious root then removes a genuine mode
     at every longer period. `tprep96`/`tdisp96`/`t?egn96` are therefore run
     separately for every period.
   - **Branch validation.** Every root is also computed at frequencies 2%
     higher and 2% lower. Two roots at adjacent frequencies belong to the
     same branch if the secant group velocity dω/dk between them matches the
     mean of their energy-integral group velocities within 1.2%. Genuine
     roots match to ~1e-6 to 3e-3; spurious roots miss by several percent.
     Roots exactly at the (flattened) half-space S velocity, the edge of the
     `tdisp96` search range, are discarded. A branch is kept only if it is
     found consistently at ≥ 2 of the 3 frequencies, and duplicates are
     removed. If a kept branch was missed at the requested frequency, it is
     rebuilt from the two neighbours (second-order accurate). The requested
     mode is then taken from the ordered list of kept branches. A
     `TI:spuriousRoot` warning reports any roots that were discarded or
     rebuilt.

   In tests on randomly perturbed isotropic-equivalent models, the TI
   phase velocities matched the isotropic `sdisp96` values to < 3e-4 at
   every period after this processing.
2. **Love waves with Vsh > Vsv.** `tdisp96` limits its phase-velocity
   search to the largest SV velocity in the model, for Love waves as well,
   so Love modes faster than that are never found. With ξ = 1.1 in the
   mantle, for example, the Love fundamental mode stopped at ~55 s. For
   Love waves `calc_der96_TI` passes `-cmax` set to the largest
   (flattened) SH velocity.
3. **eta ≠ 1 in the half-space.** This can produce a spurious Rayleigh root
   at low phase velocity (e.g. ~1.4 km/s with an unphysical group velocity).
   `check_halfspace_eta_TI` warns about it. Use eta = 1 in the half-space.
4. **Run time.** `tdisp96` (Rayleigh) is about 4× slower than `sdisp96`
   and takes nearly all of the time; `tprep96`, `tregn96`, reading the
   output and launching the programs are negligible. By default,
   `calc_der96_TI` computes three frequencies per period (two with
   `opts.nfreq = 2` or `calc_kernel96_TI`'s `opts.phase_only`, one with
   `opts.c_guess`). It also computes only the modes needed (`opts.nmodes`,
   default `nmode+1`), and recomputes a period once with 2 more modes if no
   valid root is found there (see `calc_der96_TI.m`).
5. **Root precision.** `tdisp96` refines each root with only 15 bisection
   steps, so phase velocities are good to ~1e-5 to 1e-4 km/s. Use finite
   differences of `dispR_surf96_TI`/`dispL_surf96_TI` with care.
6. **Q.** `tregn96`/`tlegn96` apply the causal-Q correction to the phase
   velocity but not to the group velocity. `calc_der96_TI` therefore always
   runs CPS elastically and computes γ and the causal phase velocity itself,
   using the formulas from `tregn96` (subroutine `gammap`). The results agree
   with the CPS values to ~1e-7. As in CPS, group velocity, eigenfunctions
   and kernels are the elastic values.

## Requirements

The CPS TI programs `tprep96`, `tdisp96`, `tregn96` and `tlegn96` must be on
the PATH, along with `gtimeout` (on mac, `brew install coreutils`; on Linux,
link `gtimeout` to `timeout`). The functions run in MATLAB and in Octave
(tested with Octave 8.4).
