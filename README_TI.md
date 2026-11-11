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

> **Note on the isotropic `calc_kernel96`.** Its 6th and 7th outputs
> (`dcdrho`, `dudrho`), and the `dcdrho`/`dudrho` fields of
> `calc_anelastic_kernel96`, are read from columns 7–8 of the srfker96 output.
> Those columns are **dc/dh and dU/dh (layer thickness)**, not density. The TI
> functions return true density kernels as `dcdrho`/`dudrho` and thickness
> kernels as `dcdh`/`dudh`.

## Benchmarks

```matlab
benchmark_TI_run_all   % or run benchmark_TI_dispersion, _kernels, _eigenfunctions, _anelastic_kernels
```

Each script compares the isotropic wrapper with the TI wrapper for a TI
model converted from an isotropic oceanic model (water, sediment, crust and
mantle). The TI kernels are also checked against finite differences for a
model with real anisotropy. A PASS/FAIL line is printed for each comparison.
The summary and the figures are written to `benchmark_TI_figs/`. To use CPS
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

The velocity kernels differ by a depth-dependent factor of about `a/(a-z)`
(about 2% at 150 km depth). `tregn96` applies the earth-flattening chain
rule to the kernels and srfker96 does not. Finite differences of the
spherical-earth phase velocity agree with the **TI** kernels, so for
spherical-earth models the isotropic srfker96 velocity kernels are slightly
too small at depth. With flat-earth models the TI kernels, including the
computed group-velocity kernels, match srfker96 to about 1e-4 to 1e-3.

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
4. **Run time.** `tdisp96` (Rayleigh) is several times slower than
   `sdisp96`, and three frequencies are computed for each requested
   period. Expect a few seconds per call for ~40 layers and ~10–20 periods.
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
(tested with Octave 8.4). Note that the isotropic `read_eig96`, which
`calc_eigenfunctions96` uses, calls `contains`, which Octave 8 does not have.
For that reason `benchmark_TI_eigenfunctions` needs MATLAB, or a `contains`
shim, in Octave 8.
