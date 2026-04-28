# pC_rate_sim.C — Code Explanation

This document explains the structure, physics, and simulation strategy of the
ROOT macro `pC_rate_sim.C`, which simulates the BNL Booster recoil-carbon
polarimeter for the polarized proton and ³He ion beams.

Reference paper: Shmakova et al., *"An absolute recoil-carbon polarimeter for
polarized light-ions at the BNL Booster"* (April 2026 draft).

---

## 1. What the simulation computes

For each beam energy point, the macro produces:

1. **Rate predictions** at each detector ring (recoil and forward).
2. **Hit position maps** on each silicon pad station, with realistic charge
   sharing across pixels (cluster physics).
3. **Kinematic distributions**: T_C vs θ_C elastic locus, TOF, projectile and
   recoil angle spectra.
4. **Smearing diagnostics**: how much the carbon recoil's hit position is
   displaced by multiple scattering inside the carbon target fiber.

The default driver is `run_all(n_events, pad_pitch_um, smearing_on)`, which
runs four energy points: pC at 500/850/1200 MeV, ³HeC at 443 MeV.

---

## 2. Physics ingredients

### 2.1 Cross sections

- **Bonin parametrization** (`namespace Bonin`) — analytic dσ/dΩ and analyzing
  power Ay for inclusive p+¹²C in the proton lab frame, valid 1.0–1.9 GeV/c.
  Coefficients from paper Tables 2 and 3.
- **Kamiya placeholder** (`namespace Kamiya`) — for ³HeC at 443 MeV. Uses
  paper Table 13's acceptance-averaged value (30 mb/sr on the forward cone)
  rather than the full digitized angular dependence. Sufficient for rate.

### 2.2 Two-body kinematics

`MakeCMKin(T, m1, m2)` builds the CM frame for elastic scattering:
- `p_cm`, `E1_cm`, `E2_cm` from the standard relativistic invariants.
- `gamma_cm`, `beta_cm` for the lab boost.

### 2.3 Beam parameters

`MakeBeamTarget(species, T)` returns the rate-formula inputs from paper Table 1:
- `N_cyc` — protons or helions per cycle.
- `f_rev` — revolution frequency, interpolated linearly in β between paper's
  injection/extraction values (paper Table 1 only quotes endpoints).
- `n_t = d_t · N_A / A_C` — target areal atomic density (Eq. 46).
- `f_ov` — fraction of beam intercepted by the thin fiber target (Eq. 53).

### 2.4 Detector geometry (`namespace Det`)

Six identical stations per ring, equally spaced in azimuth (φ = 0°, 60°, ...,
300°). All from paper Table 6:

- **Recoil ring**: 40×40 mm² Si pads, r=100 mm, z=17.4 mm, face tilted
  perpendicular to recoil direction at θ_C = 80°.
- **Forward ring**: 98×98 mm² Si pads, z=750 mm, perpendicular to beam,
  inner edge at r=100 mm.

### 2.5 Target smearing (`ApplyTargetSmearing`)

The carbon target is a vacuum-deposited C ribbon (~50 nm × ~10 μm), modeled
as a round fiber of equal cross-sectional area: deff = 800 nm diameter
(Eq. 48). Implemented effects:

1. **Vertex** sampled uniformly inside the fiber cross-section (disc of
   radius deff/2 in the x–z plane; fiber axis along +y).
2. **Path length** to fiber exit found by analytic ray–cylinder intersection.
3. **Multiple scattering** (Highland formula) applied as Gaussian kicks in
   the two axes perpendicular to the recoil direction. Z = 6 for carbon ion.
4. **Energy loss** linear dE/dx ≈ 0.10 MeV/μm (negligible for 800 nm).
5. **New origin** for pad ray-tracing = exit point of the fiber.

### 2.6 Stopping power and charge sharing in Si (`namespace CinSi`)

SRIM-2013 lookup table (paper Table 9) for ¹²C in Si:
- `dEdx_MeV_per_um(T)` — interpolates dE/dx in log(T).
- `Range_um(T)` — interpolates range in log(T).
- `Edep_MeV(T, thickness)` — total deposited energy (full T_C if stopped,
  else dE/dx × thickness if punch-through).

The pad-pixel charge distribution is computed analytically: for each event,
the Gaussian charge cloud (effective σ = 25 μm, lumping drift diffusion +
δ-rays + capacitive coupling) is integrated over each pixel using `erf`.
Pixels above a configurable threshold (3 keV, Timepix3-like) "fire"; this
gives realistic 5–10 pixel clusters (matching e.g. Bergmann et al., DOI
10.1016/j.radmeas.2024.107086).

---

## 3. Monte Carlo strategy

`RunSim` generates events flat in cos(θ_CM) ∈ [-1, +1] and φ ∈ [0, 2π].
Each event carries a weight that, when summed and scaled, gives the cross
section (and from there, the rate).

For each event:

1. **Boost projectile to lab** — get θ_proj,lab, p_lab.
2. **Compute Jacobian** dΩ_lab/dΩ_CM exactly:
   ```
   dcos_lab/dcos_cm = γ p*² (p* + β E* cos_cm) / p_lab³
   ```
   (Note: this differs slightly from paper Eq. 34 by factors of γ(1+r·cos);
   for these kinematics the difference is <1% since γ_CM ≈ 1.0–1.01.)
3. **Evaluate dσ/dΩ_lab** (Bonin or Kamiya).
4. **Compute event weight** w = dσ/dΩ_lab × (dΩ_lab/dΩ_CM).
5. **Boost carbon recoil to lab** — get θ_C,lab, T_C,lab, recoil direction.
6. **Apply target smearing** to the carbon recoil (if enabled): MS, dE,
   exit origin.
7. **Pad intersection** for both nominal (no smearing) and smeared
   trajectories, in the nearest of the 6 stations.
8. **Charge sharing**: spread the deposited energy across pixels using
   analytic erf integrals over the 2D Gaussian charge cloud.
9. **Fill histograms**.

### 3.1 Weight integration

The MC integrates dσ over the full CM solid angle:

```
σ(region) = (4π / N_events) × Σ (weight in region)
```

with weight in mb·sr⁻¹ × sr (dimensionless Jacobian). This gives σ in mb,
which converts to rate via:

```
Rate[Hz] = N_cyc × f_rev × n_t × f_ov × σ[cm²]
```

### 3.2 Two acceptance definitions per ring

To separate physics from geometry, every event is classified by:

- **2π-azimuth θ-only**: θ_C ∈ [68.7°, 91.3°] for recoil, etc.
  Gives the paper's `<dσ/dΩ> × Ω_ring` value (ideal full-azimuth integral).
- **Physical 6-pad acceptance**: the carbon trajectory actually intersects
  one of the 6 physical 40×40 mm² pad faces. Gives the realistic 6-station
  rate as the detector sees it.

The ratio between these tells you the geometric overlap factor (~6×7.3°/360°
for the forward ring, ~6×22.6°/360° for the recoil ring).

### 3.3 Smearing run as A/B comparison

When `smearing_on = kTRUE`, every event is processed twice: once with the
nominal recoil direction (origin = IP), once with the smeared direction
(origin = fiber exit). Stored in parallel histograms so the per-event
displacement (Δu, Δv) can be plotted directly. This is what
`h_recoil_pad_diff` and the `h_TC_du`, `h_TC_dv` plots show.

---

## 4. Histogram outputs

All histograms are scaled to **rate** (counts/s) at the end. The conversion
factor is:

```
scale_counts_to_rate = N_cyc × f_rev × n_t × f_ov × (4π / N_events) × mb→cm²
```

For 1D rate-density histograms (vs angle, energy, TOF), the bin width is
divided out so units are e.g. Hz/deg, Hz/MeV, Hz/ns.

### 4.1 Kinematic histograms (always filled)

- `h_theta_C` — rate vs θ_C,lab [deg]
- `h_TC` — rate vs T_C,lab [MeV]
- `h_theta_p` — rate vs θ_proj,lab [deg]
- `h_tof` — rate vs TOF [ns]
- `h_TC_th` — 2D elastic locus T_C vs θ_C
- `h_ToF_vs_Edep` — 2D deposited energy vs TOF (particle ID handle)

### 4.2 Pad-level histograms

For each of the recoil and forward rings:
- `h_*_pad` — **hit-rate** map: event rate per spatial bin [Hz/pad] at the
  smeared hit position. Does not account for cluster spreading.
- `h_*_pad_charge` — charge density map [keV·Hz / pixel]: pixel firing rate
  map after spreading the deposited energy over pixels with the analytic
  Gaussian erf integrals (pixels above threshold counted).
- `h_*_pad_nosm` — same as `h_*_pad` but without target smearing.
- `h_*_pad_diff` — per-event displacement (smeared − nominal).

### 4.3 Cluster physics

- `h_clusterSize` — distribution of N_pixels per carbon hit.
- `h_clusterSize_vs_TC` — cluster size vs T_C (cluster size grows with
  deposited energy, since more pixels in the diffusion tails clear threshold).

---

## 5. How to extend / modify

### Adding a new energy
Just call `pC_rate_sim(T, n_events)` or `He3C_rate_sim(T, n_events)`. The
rest is automatic.

### Tuning detector geometry
All recoil/forward ring numbers are in `namespace Det`. Change them and
the pad intersection updates automatically.

### Tuning smearing physics
`SmearParams` struct (file scope, above `RunSim`). Defaults: ρ=1.5 g/cm³
(amorphous C), X₀=42.7 g/cm² (graphite), Z=6 (carbon ion), dE/dx=0.1 MeV/μm.
Instantiated inside `RunSim` and populated from `BeamTarget::deff_nm`.

### Tuning charge-sharing
`sigma_diff_mm`, `Q_thr_keV` constants near the top of the per-event
charge-sharing block. Should be tuned to reproduce measured cluster-size
distributions for the chosen sensor (Timepix3, Timepix4, etc.).

### Replacing Kamiya with full digitization
Replace `Kamiya::dSig_mb` to interpolate digitized data points from
Kamiya et al. Fig. 5 (after converting cm→lab via the Jacobian Eq. 34).

---

## 6. Known limitations

1. **Bonin validity range**: 1.0–1.9 GeV/c (T_p ≈ 433–1181 MeV). At 200 MeV
   (Linac injection energy) the parametrization is extrapolated; the macro
   warns but proceeds.
2. **Paper Table 13 inconsistency**: at pC energies, Bonin's integrated
   `<dσ/dΩ>` over the recoil ring is a factor ~13 below paper Table 13.
   ³HeC matches the paper exactly. Discrepancy is on the paper side.
3. **dE/dx in Si** uses log-T linear interpolation of SRIM points
   (~5% accuracy). Sufficient for cluster-size estimation.
4. **Charge-cloud σ is empirical** — drift diffusion alone gives ~8 μm,
   but real Timepix3 measurements show effective σ ≈ 25–35 μm including
   δ-ray range and other broadening. Tuning required against data.
5. **Forward-ring projectile MS** is not simulated (~20 μm displacement
   at z=750 mm, much smaller than pad pitch).
6. **No edge effects**: events landing on pad edges are accepted as hits;
   no dead-region around pad borders is modeled.

---

## 7. File layout

```
pC_rate_sim.C
├── PC namespace               — physical constants
├── Bonin namespace            — pC cross section + Ay
├── Kamiya namespace           — 3HeC placeholder cross section
├── CinSi namespace            — SRIM stopping power lookup
├── CMKin struct + MakeCMKin   — relativistic kinematics
├── BeamTarget struct + funcs  — Table 1 beam parameters
├── Det namespace              — detector geometry constants
├── SmearParams struct         — fiber smearing configuration (file scope)
├── ApplyTargetSmearing        — fiber MS + dE
├── HitRecoilPad / HitFwdPad   — ray–pad intersection (with origin)
├── SimResult struct           — all output histograms + summary numbers
├── RunSim                     — main MC loop (this is the work)
├── MakeKinPlot / MakePadPlot / MakeDiffPlot / MakeChargePlot  — drawing
└── pC_rate_sim / He3C_rate_sim / run_all   — top-level driver functions
```

---

## 8. Build & run

Requires ROOT (any version with TVector3 — i.e. ≥5).

```
root -l
.L pC_rate_sim.C+
run_all()                               // 4 energies, 5M events each (default)
run_all(1000000, 55.0, kTRUE)           // 1M events, 55 um pitch, smearing ON
run_all(500000, 55.0, kFALSE)           // disable smearing globally
pC_rate_sim(500., 500000)     // single proton energy
He3C_rate_sim(443., 500000)             // single 3He energy
```

Each call produces four output files per energy point:
- `kin_<species>_<T>MeV[_smear|_nosm].png/.pdf`   — kinematic distributions
- `pads_<species>_<T>MeV[_smear|_nosm].png/.pdf`  — pad rate maps + displacement
- `diff_<species>_<T>MeV[_smear|_nosm].png/.pdf`  — smeared vs nominal projections
- `charge_<species>_<T>MeV[_smear|_nosm].png/.pdf` — charge/cluster distributions

and prints a summary table comparing rates against paper Table 13.