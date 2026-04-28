# toy_MC_pp_pC_CNI_RHIC.ipynb — Code Description

Jupyter notebook simulating **proton-proton (pp) and proton-carbon (pC) elastic
scattering in the Coulomb-Nuclear Interference (CNI) region** for the RHIC HJET
polarimeter. A 100 GeV kinetic-energy proton beam hits a target at rest (either
a proton jet or a thin carbon fiber); recoil particles are tracked through the
HJET dipole magnet to silicon strip detectors at R = 220 mm.

Reference paper: Shmakova et al., *"An absolute recoil-carbon polarimeter for
polarized light-ions at the BNL Booster"* (April 2026 draft).

---

## What the notebook computes

1. **CNI event generation** — flat MC sample in momentum transfer *t*, boosted
   to the lab frame; *t*-cut selects the CNI region (0.001 < |t| < 0.03 GeV²).
2. **Analyzing power A_N(t)** — CNI model with electromagnetic–hadronic
   interference and proton form factor (Z=1, pp case).
3. **Kinematic distributions** — θ, T of scattered and recoil particles in both
   CM and lab frames; TOF to the detector plane at R = 22 cm.
4. **Magnetic-field propagation** — Boris-push integrator tracks recoil protons
   through Bx, By, Bz field maps read from CSV files; hit positions compared
   with and without each field component.
5. **Silicon detector response** — PSTAR stopping-power table for protons in Si;
   energy deposition in a 300 μm sensor; TOF vs E_dep correlation.
6. **Hit-position analysis** — left/right asymmetry masks; T vs z projections;
   charge calculation (w = 3.65 eV/pair, Si capacitance C_SSD = 17 pF).
7. **Detector digitization** — hit positions binned onto pixel (0.2 mm pitch)
   and strip (64 strips / 100 mm) grids; 2D occupancy maps.
8. **3D geometry visualization** — Python translation of `xzy_sketch_3D_det250.m`
   (MATLAB); draws Bx/By/Bz coil geometry and two detector planes at r = 220 mm.

---

## Physics ingredients

### Beam and target

| Parameter | Value |
|---|---|
| Beam kinetic energy | 100 GeV (proton) |
| Target | proton (pp) or ¹²C (pC), at rest |
| Detector radius R | 22 cm (220 mm) |
| CNI *t*-window | 0.001 – 0.03 GeV² |

### Analyzing power model

```python
def A_N_model(t):
    # CNI: electromagnetic-hadronic interference
    # G = 1.7928 (proton anomalous moment)
    # alpha = 1/137, sig_tot = 35 mb (pp), Z = 1
    t0 = 8*pi*alpha*Z / sig_tot
    A_N = G * t0 * t * sqrt(t) / (mp * (t**2 + t0**2))
```

Valid for |t| ≲ 0.1 GeV²; used for event weighting in asymmetry analysis.

### B-field propagation (Boris push)

```python
def boris_push(r, v, dt, use_B=True):
    # Lorentz force integration using the Boris algorithm
    # q = proton charge, m = proton mass
```

Field maps loaded from `Bfield_map_B{x,y,z}_Pair` CSV files (3D grids).
Four scenarios compared: no B, BX, BY, BZ.

### Energy loss in silicon

PSTAR tabulated stopping powers (MeV cm²/g) interpolated in log(T).
Deposited energy in a 300 μm Si layer:
- particle stops → full T deposited
- punch-through → dE/dx × thickness

---

## Notebook structure

| Section | Cells | Description |
|---|---|---|
| Setup | 0–1 | Imports, matplotlib backend, constants |
| CNI generation | 2–4 | Event generation, lab boost, TOF |
| A_N model | 6–9 | Analyzing power function, *t* selections |
| Lab/CM plots | 10–16 | θ, T distributions, elastic locus, 3D hit map |
| B-field tracking | 17, 23 | Boris pusher, field map loading |
| Si detector | 18–22 | Geometry, PSTAR dE/dx tables |
| Analysis | 25–41 | CSV input, TOF/E_dep 2D histograms, charge, asymmetry |
| 3D geometry | 42–58 | HJET magnet coils + detector visualization |
| Hit analysis | 59–62 | Left/right masks, T vs hit-position plots |
| Digitization | 63–77 | Pixel/strip binning, occupancy maps |
| Extras | 78–81 | 3D voxelization (commented), resolution smearing (commented) |

---

## Input files

| File pattern | Contents |
|---|---|
| `FieldMap_B{X,Y,Z}_tgt_orientation.csv` | 3D B-field maps for the HJET magnet |
| `Bfield_map_B{x,y,z}_Pair` (CSV) | Alternative field map set |
| `proton_drift_comparison_Vertex_1cm_noB_R22.0.csv` | Pre-generated pp hit data, no B |
| `proton_drift_comparison_Vertex_1cm_withB_Bfield_map_B{x,y,z}_Pair_R22.0.csv` | Pre-generated hit data with each B component |

All CSV files must be present in the working directory before running the
Analysis section (cell 26 onward).

---

## Key tunable parameters

```python
# Cell 3 — kinematics
R        = 22.0     # cm  — distance from IP to detector plane
t_cut    = 0.02     # GeV² — upper t limit for generation
N_steps  = 5000     # Boris-push steps per track
time_step = 1e-11   # s   — Boris time step

# Cell 9 — CNI window
selection = (t < 0.03) & (t > 0.001)  # GeV²

# Cell 45 — detector geometry
r_det    = 220      # mm — detector radius
detSize  = 100      # mm — square pad side length

# Cell 66 — digitization
bin_mm   = 0.20     # mm — pixel pitch

# Cell 72 — strip detector
n_stripes_per_side = 64
total_width_mm     = 100
```

---

## Dependencies

```
numpy
matplotlib
scipy          # RegularGridInterpolator for field maps
pandas         # CSV input in analysis section
ipympl         # optional — interactive 3D rotation (%matplotlib widget)
```

---

## Relation to pC_rate_sim.C

This notebook and `pC_rate_sim.C` address complementary aspects of the same
program:

| | `toy_MC_pp_pC_CNI_RHIC.ipynb` | `pC_rate_sim.C` |
|---|---|---|
| Machine | RHIC (100 GeV) | BNL Booster (0.2–1.5 GeV) |
| Reaction | pp and pC CNI | pC and ³HeC elastic |
| Detector R | 220 mm | 100 mm |
| B-field tracking | Yes (Boris push) | No |
| Rate normalization | No (raw counts) | Yes (Hz/pad) |
| Output format | Python plots | ROOT histograms + PDFs |

---

## Known limitations / open items

1. **A_N model uses Z=1** (proton target). For pC the hadronic amplitude and
   nuclear form factor need to be updated for carbon.
2. **PSTAR tables are for protons** in Si. For carbon-recoil detection use the
   SRIM ¹²C-in-Si tables from `pC_rate_sim.C`'s `CinSi` namespace instead.
3. **No rate normalization** — all histograms show raw MC counts, not physical
   rates. Use `pC_rate_sim.C` for absolute rate predictions.
4. **Field-map CSV files** are required for the B-field sections but are not
   committed to the repository; those cells will fail without them.
5. **3D voxelization and resolution smearing** (cells 78–81) are commented out.
