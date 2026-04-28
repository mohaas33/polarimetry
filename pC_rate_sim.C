// =============================================================================
//  pC_rate_sim.C  —  Booster polarimeter rate simulation  (with target smearing)
// =============================================================================
//
//  Simulates:
//      (1) p + 12C elastic at T = 500, 850, 1200 MeV       (Bonin)
//      (2) 3He + 12C elastic at T = 443 MeV                (Kamiya avg. 30 mb/sr)
//
//  Target smearing (NEW vs earlier version):
//      - For each event, sample vertex uniformly inside the fiber cross section
//        (diameter deff = 800 nm, round-fiber equivalent of Eq. 48 in the paper).
//      - Propagate carbon recoil to the fiber exit along its direction.
//      - Apply Highland multiple scattering (MS) as a Gaussian kick perpendicular
//        to the recoil direction.
//      - Apply a simple dE/dx energy loss along the path.
//      - Use the exit point as the origin for the pad ray-intersection.
//      - Forward projectile: MS through 800 nm of C is ~27 μrad -> <25 μm at
//        z=750 mm.  Negligible vs. forward pad segmentation.  Skip for speed.
//
//  Reference: Shmakova et al., Apr 2026.
//
//  Usage (ROOT):
//      .L pC_rate_sim.C+
//      run_all()                              // all 4 energies, 500k events
//      run_all(1000000, 55.0, kTRUE)          // more stats, 55um pitch, smearing on
//      run_all(500000, 55.0, kFALSE)          // disable smearing globally
//      pC_rate_sim(500., 500000)
//      He3C_rate_sim(443., 500000)
// =============================================================================

#include "TMath.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TString.h"
#include "TVector3.h"
#include "TLegend.h"
#include <vector>

// =============================================================================
//  Constants
// =============================================================================
namespace PC {
    const Double_t mp     = 938.272;
    const Double_t m3He   = 2808.391;
    const Double_t mC     = 12.0 * 931.494;
    const Double_t NA     = 6.022e23;
    const Double_t AC     = 12.0;
    const Double_t mb2cm2 = 1.0e-27;
    const Double_t c_cms  = 29.9792458;   // cm/ns
}

// =============================================================================
//  Bonin parametrization  (inclusive p+12C, proton LAB frame)
// =============================================================================
namespace Bonin {
    const Double_t Bx = -0.880;
    const Double_t Cx =  0.0567;
    const Double_t Dx = -0.00139;
    const Double_t G0 = -2.98e-23;
    const Double_t G1 =  9.32e-23;
    const Double_t G2 = -2.39e-24;

    const Double_t a[4] = {  1.626, -1.786,   1.933,  -0.285 };
    const Double_t b[4] = { 14.909, -0.594, -91.766, 207.064 };
    const Double_t c[4] = {574.719,-560.410,-610.026,1154.660};
    const Double_t d[4] = {  0.1275,  0.0550,  -0.138,  -0.639 };

    inline Double_t poly3(const Double_t k[4], Double_t x) {
        return k[0] + k[1]*x + k[2]*x*x + k[3]*x*x*x;
    }

    Double_t dSig_mb(Double_t th_lab_deg, Double_t p_GeVc) {
        Double_t G  = G0 + G1*p_GeVc + G2*p_GeVc*p_GeVc;
        Double_t ex = Bx*th_lab_deg + Cx*th_lab_deg*th_lab_deg
                      + Dx*th_lab_deg*th_lab_deg*th_lab_deg;
        Double_t xs = G * TMath::Exp(ex);
        return (xs > 0) ? xs / PC::mb2cm2 : 0.0;
    }

    Double_t Ay(Double_t th_lab_deg, Double_t p_GeVc) {
        Double_t pp = p_GeVc - 1.4;
        Double_t r  = p_GeVc * TMath::Sin(TMath::DegToRad()*th_lab_deg);
        Double_t aa = poly3(a, pp);
        Double_t bb = poly3(b, pp);
        Double_t cc = poly3(c, pp);
        Double_t dd = poly3(d, pp);
        Double_t den = 1.0 + bb*r*r + cc*r*r*r*r;
        if (TMath::Abs(den) < 1e-9) return 0.0;
        return (aa*r)/den
             + dd*p_GeVc*TMath::Sin(TMath::DegToRad()*5.0*th_lab_deg);
    }

    inline Bool_t InRange(Double_t p_GeVc) {
        return (p_GeVc >= 1.0 && p_GeVc <= 1.9);
    }
}

// =============================================================================
//  Kamiya placeholder  (paper Table 13 acceptance-averaged value)
// =============================================================================
namespace Kamiya {
    const Double_t dSig_avg_mb = 30.0;
    const Double_t Ay_avg      = 0.63;

    Double_t dSig_mb(Double_t th_proj_lab_deg) {
        if (th_proj_lab_deg >= 7.5 && th_proj_lab_deg <= 14.8)
            return dSig_avg_mb;
        return 0.0;
    }
}


// =============================================================================
//  Stopping power and range for 12C in Si  (SRIM-2013, paper Table 9)
// =============================================================================
namespace CinSi {
    // T [MeV], dE/dx [MeV/um], Range [um]
    const Int_t    N_pts = 10;
    const Double_t T_tab[N_pts]    = {  1.0,   5.0,  10.0,  20.0,  30.0,
                                       50.0,  72.0, 100.0, 150.0, 197.0 };
    const Double_t dEdx_tab[N_pts] = {0.043, 0.093, 0.125, 0.146, 0.152,
                                      0.147, 0.135, 0.119, 0.097, 0.083 };
    const Double_t Rng_tab[N_pts]  = {  8.0,  38.0,  72.0, 149.0, 228.0,
                                      415.0, 625.0, 930.0,1520.0,2200.0 };

    // Linear interpolation in log10(T); flat extrapolation outside the table
    Double_t Interp(const Double_t* y, Double_t T_MeV) {
        if (T_MeV <= T_tab[0])         return y[0];
        if (T_MeV >= T_tab[N_pts-1])   return y[N_pts-1];
        Double_t lT = TMath::Log10(T_MeV);
        for (Int_t i = 0; i < N_pts-1; ++i) {
            if (T_MeV >= T_tab[i] && T_MeV < T_tab[i+1]) {
                Double_t l1 = TMath::Log10(T_tab[i]);
                Double_t l2 = TMath::Log10(T_tab[i+1]);
                Double_t f  = (lT - l1) / (l2 - l1);
                return y[i] + f * (y[i+1] - y[i]);
            }
        }
        return y[N_pts-1];
    }
    inline Double_t dEdx_MeV_per_um(Double_t T_MeV) { return Interp(dEdx_tab, T_MeV); }
    inline Double_t Range_um       (Double_t T_MeV) { return Interp(Rng_tab,  T_MeV); }

    // Energy actually deposited in a Si sensor of given thickness:
    //   - if range < thickness: full kinetic energy (carbon stops inside)
    //   - else: dE/dx * thickness (a thin-detector approximation; OK because
    //          for these C kinematics dE/dx varies <30% across 500 um)
    Double_t Edep_MeV(Double_t T_MeV, Double_t Si_thickness_um = 500.0) {
        Double_t R = Range_um(T_MeV);
        if (R <= Si_thickness_um) return T_MeV;
        return dEdx_MeV_per_um(T_MeV) * Si_thickness_um;
    }
}


// =============================================================================
//  CM kinematics
// =============================================================================
struct CMKin {
    Double_t m1, m2;
    Double_t T;
    Double_t E1_lab;
    Double_t p1_lab;
    Double_t p1_lab_GeVc;
    Double_t sqrts;
    Double_t beta_cm, gamma_cm;
    Double_t p_cm;
    Double_t E1_cm, E2_cm;
};

CMKin MakeCMKin(Double_t T_MeV, Double_t m1, Double_t m2) {
    CMKin k;
    k.m1 = m1; k.m2 = m2; k.T = T_MeV;
    k.E1_lab = T_MeV + m1;
    k.p1_lab = TMath::Sqrt(k.E1_lab*k.E1_lab - m1*m1);
    k.p1_lab_GeVc = k.p1_lab / 1000.0;
    Double_t s = m1*m1 + m2*m2 + 2.0*m2*k.E1_lab;
    k.sqrts    = TMath::Sqrt(s);
    k.beta_cm  = k.p1_lab / (k.E1_lab + m2);
    k.gamma_cm = (k.E1_lab + m2) / k.sqrts;
    k.p_cm     = TMath::Sqrt((s - (m1+m2)*(m1+m2))
                             *(s - (m1-m2)*(m1-m2))) / (2.0*k.sqrts);
    k.E1_cm    = TMath::Sqrt(k.p_cm*k.p_cm + m1*m1);
    k.E2_cm    = TMath::Sqrt(k.p_cm*k.p_cm + m2*m2);
    return k;
}

// =============================================================================
//  Beam parameters (paper Table 1)
// =============================================================================
struct BeamTarget {
    TString  species;
    Double_t N_cyc;
    Double_t f_rev;
    Double_t n_t;
    Double_t f_ov;
    Double_t sigma_beam_mm;
    Double_t deff_nm;        // effective round-fiber diameter
};

static Double_t beta_of_T(Double_t T, Double_t m) {
    Double_t E = T + m;
    Double_t p = TMath::Sqrt(E*E - m*m);
    return p / E;
}

Double_t f_rev_proton(Double_t T) {
    const Double_t f_inj = 0.841e6, f_ext = 1.371e6;
    const Double_t T_inj = 200.0,   T_ext = 1500.0;
    Double_t b_inj = beta_of_T(T_inj, PC::mp);
    Double_t b_ext = beta_of_T(T_ext, PC::mp);
    Double_t b     = beta_of_T(T,     PC::mp);
    Double_t frac  = (b - b_inj) / (b_ext - b_inj);
    frac = TMath::Max(0.0, TMath::Min(1.0, frac));
    return f_inj + frac * (f_ext - f_inj);
}

Double_t f_rev_helion(Double_t T) {
    const Double_t f_inj = 0.097e6, f_ext = 1.260e6;
    const Double_t T_inj = 6.0,     T_ext = 2496.0;
    Double_t b_inj = beta_of_T(T_inj, PC::m3He);
    Double_t b_ext = beta_of_T(T_ext, PC::m3He);
    Double_t b     = beta_of_T(T,     PC::m3He);
    Double_t frac  = (b - b_inj) / (b_ext - b_inj);
    frac = TMath::Max(0.0, TMath::Min(1.0, frac));
    return f_inj + frac * (f_ext - f_inj);
}

BeamTarget MakeBeamTarget(TString species, Double_t T_MeV,
                           Double_t dt_ugcm2 = 4.0,
                           Double_t eN_um    = 1.0,
                           Double_t beta_f_m = 7.75,
                           Double_t deff_nm  = 800.0)
{
    BeamTarget bt;
    bt.species = species;
    bt.deff_nm = deff_nm;

    Double_t m = 0;
    if (species == "p") {
        bt.N_cyc = 1.5e11;
        bt.f_rev = f_rev_proton(T_MeV);
        m = PC::mp;
    } else if (species == "3He") {
        bt.N_cyc = 1.0e11;
        bt.f_rev = f_rev_helion(T_MeV);
        m = PC::m3He;
    } else {
        ::Error("MakeBeamTarget", "Unknown species '%s'", species.Data());
        bt.N_cyc = 0; bt.f_rev = 0; m = PC::mp;
    }

    bt.n_t = (dt_ugcm2 * 1.0e-6) * PC::NA / PC::AC;

    Double_t E  = T_MeV + m;
    Double_t p  = TMath::Sqrt(E*E - m*m);
    Double_t bg = p / m;
    Double_t eN = eN_um * 1.0e-6;
    Double_t e_geo   = eN / bg;
    Double_t sigma_m = TMath::Sqrt(e_geo * beta_f_m);
    bt.sigma_beam_mm = sigma_m * 1000.0;

    Double_t deff_cm  = deff_nm * 1.0e-7;
    Double_t sigma_cm = sigma_m * 100.0;
    bt.f_ov = deff_cm / (TMath::Sqrt(2.0*TMath::Pi()) * sigma_cm);

    return bt;
}

// =============================================================================
//  Detector geometry — physical flat pads (paper Table 6)
// =============================================================================
namespace Det {
    const Double_t recoil_r_mm      = 100.0;
    const Double_t recoil_z_mm      =  17.4;
    const Double_t recoil_side_mm   =  40.0;
    const Double_t recoil_theta_deg =  80.0;
    const Double_t recoil_thC_lo    =  68.7;
    const Double_t recoil_thC_hi    =  91.3;

    const Double_t fwd_z_mm         = 750.0;
    const Double_t fwd_r_inner_mm   = 100.0;
    const Double_t fwd_side_mm      =  98.0;
    const Double_t fwd_theta_lo     =   7.5;
    const Double_t fwd_theta_hi     =  14.8;

    const Int_t    n_stn            = 6;

    inline Double_t stn_phi(Int_t i) { return i * TMath::Pi() / 3.0; }

    inline Double_t nearest_stn_phi(Double_t phi) {
        const Double_t step = TMath::Pi() / 3.0;
        Double_t k = TMath::Nint(phi / step);
        return k * step;
    }
}

// =============================================================================
//  Target smearing
//
//  Fiber model (paper Sec. 5.2):
//    • Real ribbon t ≈ 50 nm × w ≈ 10 μm replaced by a round fiber of equal
//      cross-sectional area: deff ≈ 800 nm DIAMETER (Eq. 48).
//    • Axis chosen along +y (vertical); beam along +z.  By symmetry in y,
//      the y coordinate of the vertex doesn't affect anything, so set y_v=0.
//    • Vertex sampled uniformly in the fiber CROSS SECTION (circle in x–z).
//
//  Scattering model:
//    • Amorphous carbon:  ρ = 1.5 g/cm³, X0 = 42.7 g/cm²
//    • Highland MS for carbon recoil (z=6):
//         θ₀ = (13.6 MeV · z / βp) · √(x/X0) · (1 + 0.038 ln(x/X0))
//    • Gaussian kick in the plane perpendicular to direction.
//    • Linear dE/dx model (0.10 MeV/μm in carbon — tunable).  For 800 nm path,
//      ΔE ≈ 80 keV: negligible kinematic effect but we keep it for consistency.
// =============================================================================
struct SmearParams {
    Bool_t   enable           = kTRUE;
    Double_t rho_gcm3         = 1.50;    // amorphous C
    Double_t X0_gcm2          = 42.7;    // graphite (good enough)
    Double_t z_recoil         = 6.0;     // Z of carbon recoil ion
    Double_t dEdx_MeV_per_um  = 0.10;    // mean in C (tunable)
};

// Sample uniformly inside a disc of radius R (in mm)
inline void SampleDiscUniform(TRandom3& rng, double R_mm,
                              double& x_mm, double& z_mm) {
    double r = R_mm * std::sqrt(rng.Rndm());
    double phi = 2.0 * TMath::Pi() * rng.Rndm();
    x_mm = r * std::cos(phi);
    z_mm = r * std::sin(phi);
}

// Apply target smearing to a carbon recoil:
//   Updates (dx,dy,dz) in-place with MS-kicked unit direction.
//   Updates T,E,p,beta in-place with energy loss.
//   Writes exit point (mm) to origin_mm.
// Returns false if the trajectory can't exit the fiber (shouldn't happen in
// normal kinematics), in which case it still sets origin_mm = vertex and
// leaves direction unchanged.
Bool_t ApplyTargetSmearing(const SmearParams& sp, TRandom3& rng,
                           Double_t deff_nm,
                           Double_t& dx, Double_t& dy, Double_t& dz,
                           Double_t& T_C, Double_t& E_C,
                           Double_t& p_C, Double_t& beta_C,
                           TVector3& origin_mm)
{
    origin_mm.SetXYZ(0, 0, 0);
    if (!sp.enable) return kTRUE;

    // Fiber radius in mm (deff is DIAMETER):  R = deff/2
    const Double_t R_mm = 0.5 * deff_nm * 1.0e-6;

    // Sample vertex in fiber cross-section (x–z plane, fiber axis = y)
    Double_t x_v, z_v;
    SampleDiscUniform(rng, R_mm, x_v, z_v);
    const Double_t y_v = 0.0;

    // Unit direction
    TVector3 d(dx, dy, dz);
    if (d.Mag() < 1e-12) return kFALSE;
    d = d.Unit();

    // Intersect fiber cylinder  x^2 + z^2 = R^2  along the direction
    //   (x_v + λ dx)^2 + (z_v + λ dz)^2 = R^2
    const Double_t A = d.X()*d.X() + d.Z()*d.Z();
    const Double_t B = 2.0 * (x_v*d.X() + z_v*d.Z());
    const Double_t C = x_v*x_v + z_v*z_v - R_mm*R_mm;

    Double_t lam_mm = 0.0;
    const Double_t disc = B*B - 4.0*A*C;
    if (disc > 0.0 && A > 1e-18) {
        const Double_t sq = TMath::Sqrt(disc);
        const Double_t l1 = (-B + sq) / (2.0*A);
        const Double_t l2 = (-B - sq) / (2.0*A);
        // forward intersection
        lam_mm = (l1 > 0.0) ? l1 : ((l2 > 0.0) ? l2 : 0.0);
    }
    if (lam_mm <= 0.0) {
        // Direction parallel to fiber axis (A~0), or already outside.
        origin_mm.SetXYZ(x_v, y_v, z_v);
        return kFALSE;
    }

    const Double_t L_cm   = lam_mm * 0.1;             // mm -> cm
    const Double_t x_gcm2 = sp.rho_gcm3 * L_cm;
    const Double_t frac   = TMath::Max(x_gcm2 / sp.X0_gcm2, 1.0e-18);

    // Highland MS (projected RMS per plane)
    const Double_t bp_safe = TMath::Max(beta_C * p_C, 1.0e-6);
    const Double_t theta0  = (13.6 * sp.z_recoil / bp_safe)
                           * TMath::Sqrt(frac)
                           * (1.0 + 0.038 * TMath::Log(frac));

    // Two Gaussian kicks in orthogonal directions ⟂ d
    TVector3 n = d.Orthogonal().Unit();
    TVector3 m = d.Cross(n).Unit();
    Double_t g1 = rng.Gaus(0.0, theta0);
    Double_t g2 = rng.Gaus(0.0, theta0);
    TVector3 d_new = (d + g1*n + g2*m).Unit();

    // Energy loss (linear model)
    const Double_t L_um = L_cm * 1.0e4;
    Double_t dE = sp.dEdx_MeV_per_um * L_um;
    T_C = TMath::Max(0.0, T_C - dE);
    E_C = T_C + PC::mC;
    p_C = TMath::Sqrt(TMath::Max(0.0, E_C*E_C - PC::mC*PC::mC));
    beta_C = (E_C > 0.0) ? (p_C / E_C) : 0.0;

    // Update direction and exit origin (using ORIGINAL direction, since the
    // MS kick happens at the exit face — a thin-scatterer approximation)
    dx = d_new.X(); dy = d_new.Y(); dz = d_new.Z();
    origin_mm.SetXYZ(x_v + lam_mm*d.X(),
                     y_v + lam_mm*d.Y(),
                     z_v + lam_mm*d.Z());
    return kTRUE;
}

// =============================================================================
//  Pad intersection
//  Single definitions; default origin = (0,0,0) gives the "no smearing" case.
// =============================================================================
Bool_t HitRecoilPad(Double_t dx, Double_t dy, Double_t dz,
                    Double_t phi_stn,
                    Double_t& u_mm, Double_t& v_mm,
                    const TVector3& origin_mm = TVector3(0,0,0))
{
    const Double_t th_n = TMath::DegToRad() * Det::recoil_theta_deg;
    Double_t sn = TMath::Sin(th_n), cn = TMath::Cos(th_n);
    Double_t cs = TMath::Cos(phi_stn), ss = TMath::Sin(phi_stn);

    TVector3 P0(Det::recoil_r_mm * cs,
                Det::recoil_r_mm * ss,
                Det::recoil_z_mm);
    TVector3 n_face(sn*cs, sn*ss, cn);
    TVector3 u_loc(-ss,    cs,    0.0);
    TVector3 v_loc(cn*cs,  cn*ss, -sn);

    TVector3 d(dx, dy, dz);
    Double_t denom = d.Dot(n_face);
    if (denom <= 1e-9) return kFALSE;
    Double_t lam = (P0.Dot(n_face) - origin_mm.Dot(n_face)) / denom;
    if (lam <= 0) return kFALSE;

    TVector3 H = origin_mm + d * lam;
    TVector3 rel = H - P0;
    u_mm = rel.Dot(u_loc);
    v_mm = rel.Dot(v_loc);
    Double_t half = Det::recoil_side_mm / 2.0;
    return (TMath::Abs(u_mm) <= half && TMath::Abs(v_mm) <= half);
}

Bool_t HitFwdPad(Double_t dx, Double_t dy, Double_t dz,
                 Double_t phi_stn,
                 Double_t& u_mm, Double_t& v_mm,
                 const TVector3& origin_mm = TVector3(0,0,0))
{
    if (dz <= 1e-9) return kFALSE;
    Double_t lam = (Det::fwd_z_mm - origin_mm.Z()) / dz;
    if (lam <= 0) return kFALSE;

    Double_t hx = origin_mm.X() + lam * dx;
    Double_t hy = origin_mm.Y() + lam * dy;

    Double_t cs = TMath::Cos(phi_stn), ss = TMath::Sin(phi_stn);
    Double_t hx_loc =  cs*hx + ss*hy;
    Double_t hy_loc = -ss*hx + cs*hy;

    Double_t r_center = Det::fwd_r_inner_mm + Det::fwd_side_mm / 2.0;
    u_mm = hy_loc;
    v_mm = hx_loc - r_center;
    Double_t half = Det::fwd_side_mm / 2.0;
    return (TMath::Abs(u_mm) <= half && TMath::Abs(v_mm) <= half);
}

// =============================================================================
//  SimResult
// =============================================================================
struct SimResult {
    TString  tag;
    TString  species;
    Double_t T_MeV;
    Double_t p_lab_GeVc;
    BeamTarget bt;

    Double_t sigma_recoil_2pi_mb;
    Double_t sigma_fwd_2pi_mb;
    Double_t rate_recoil_2pi;
    Double_t rate_fwd_2pi;
    Double_t dsigdO_avg_recoil;
    Double_t dsigdO_avg_fwd;

    Double_t sigma_recoil_pad_mb;
    Double_t sigma_fwd_pad_mb;
    Double_t rate_recoil_pad;
    Double_t rate_fwd_pad;

    TH1D* h_theta_C = nullptr;
    TH1D* h_TC      = nullptr;
    TH1D* h_theta_p = nullptr;
    TH1D* h_tof     = nullptr;
    TH2D* h_TC_th   = nullptr;
    TH2D* h_TC_du   = nullptr;
    TH2D* h_TC_dv   = nullptr;

    // With smearing (if enabled)
    TH2D* h_recoil_pad = nullptr;
    TH2D* h_fwd_pad    = nullptr;

    // Without smearing (always filled; for comparison)
    TH2D* h_recoil_pad_nosm = nullptr;
    TH2D* h_fwd_pad_nosm    = nullptr;

    // Per-event difference: smeared hit - nominal hit (same pad frame)
    TH2D* h_recoil_pad_diff = nullptr;
    TH2D* h_fwd_pad_diff    = nullptr;

    TH2D* h_recoil_charge_vs_TC   = nullptr;   // cluster charge vs T_C
    TH2D* h_recoil_pad_charge   = nullptr;   // Hz x keV per pixel
    TH1D* h_clusterSize         = nullptr;   // pixels per cluster (rate-weighted)
    TH2D* h_clusterSize_vs_TC   = nullptr;   // cluster size vs T_C
    TH2D* h_ToF_vs_Edep = nullptr;   // deposited energy [keV] vs TOF [ns]
};

// =============================================================================
//  Core simulation
// =============================================================================
SimResult RunSim(TString species, Double_t T_MeV, Int_t n_events,
                 Double_t pad_pitch_um = 55.0,
                 Bool_t   smearing_on  = kTRUE,
                 Bool_t   verbose      = kTRUE)
{

    // --- Charge-sharing / readout parameters ---
    const Double_t Si_thickness_um   = 500.0;   // recoil sensor thickness
    const Double_t sigma_diff_mm     = 0.025;   // ~25 um EFFECTIVE cluster sigma.
                                            //   includes drift diffusion (~8 um),
                                            //   delta-ray range, and other broadening.
                                            //   Tune to match cluster size data:
                                            //   ~0.030 reproduces 5-9 pix proton clusters
                                            //   (Bergmann et al, doi:10.1016/j.radmeas.2024.107086)
    const Double_t Q_thr_keV         = 3.0;     // pixel firing threshold (3 keV typical Timepix3)
    const Int_t    pixel_window      = 6;       // search +/- N pixels around hit center
                                                //   set so it covers ~5*sigma_diff/pitch


    SimResult R;
    R.species = species;
    R.T_MeV   = T_MeV;
    R.tag     = TString::Format("%s_%.0fMeV_%s",
                                species.Data(), T_MeV,
                                smearing_on ? "smear" : "nosm");

    Double_t m_proj = (species == "p") ? PC::mp : PC::m3He;
    CMKin k = MakeCMKin(T_MeV, m_proj, PC::mC);
    R.p_lab_GeVc = k.p1_lab_GeVc;
    R.bt = MakeBeamTarget(species, T_MeV);

    SmearParams smear;
    smear.enable = smearing_on;

    Bool_t use_bonin  = (species == "p");
    Bool_t use_kamiya = (species == "3He");
    if (use_bonin && !Bonin::InRange(R.p_lab_GeVc)) {
        ::Warning("RunSim",
                  "p_lab = %.3f GeV/c outside Bonin validity [1.0, 1.9]",
                  R.p_lab_GeVc);
    }

    if (verbose) {
        Printf("\n================================================================");
        Printf("   %s + 12C  at  T = %.0f MeV   [smearing: %s]",
               species.Data(), T_MeV, smearing_on ? "ON" : "OFF");
        Printf("================================================================");
        Printf("   p_lab      = %.4f GeV/c", R.p_lab_GeVc);
        Printf("   sqrt(s)    = %.1f MeV", k.sqrts);
        Printf("   beta_cm    = %.5f   gamma_cm = %.5f", k.beta_cm, k.gamma_cm);
        Printf("   N_cyc      = %.2e", R.bt.N_cyc);
        Printf("   f_rev      = %.3f MHz", R.bt.f_rev / 1e6);
        Printf("   n_t        = %.3e cm^-2  (dt=4 ug/cm^2)", R.bt.n_t);
        Printf("   sigma_beam = %.2f mm", R.bt.sigma_beam_mm);
        Printf("   f_ov       = %.3e", R.bt.f_ov);
        Printf("   deff       = %.0f nm (fiber diameter)", R.bt.deff_nm);
        Printf("   dsigma/dOmega source: %s",
               use_bonin ? "Bonin (Eq. 21, Table 3)" :
               use_kamiya ? "Kamiya-average (30 mb/sr on forward cone)" :
               "NONE");
    }

    // -------- Kinematic histograms --------
    R.h_theta_C = new TH1D(TString::Format("h_theta_C_%s", R.tag.Data()),
        ";#theta_{C}^{lab} [deg];Rate [counts/s/deg]", 90, 0, 90);
    R.h_TC = new TH1D(TString::Format("h_TC_%s", R.tag.Data()),
        ";T_{C}^{lab} [MeV];Rate [counts/s/MeV]", 120, 0, 120);
    R.h_theta_p = new TH1D(TString::Format("h_theta_p_%s", R.tag.Data()),
        ";#theta_{proj}^{lab} [deg];Rate [counts/s/deg]", 180, 0, 90);
    R.h_tof = new TH1D(TString::Format("h_tof_%s", R.tag.Data()),
        ";TOF [ns];Rate [counts/s/ns]", 100, 0, 60);
    R.h_TC_th = new TH2D(TString::Format("h_TC_th_%s", R.tag.Data()),
        ";#theta_{C}^{lab} [deg];T_{C} [MeV]", 90, 60, 90, 120, 0, 120);
    R.h_TC_du = new TH2D(TString::Format("h_TC_du_%s", R.tag.Data()),
        ";T_{C} [MeV];#Delta u [mm]", 120, 0, 120, 200, -2.0, 2.0);
    R.h_TC_dv = new TH2D(TString::Format("h_TC_dv_%s", R.tag.Data()),
        ";T_{C} [MeV];#Delta v [mm]", 120, 0, 120, 200, -2.0, 2.0);

    // -------- Pad histograms --------
    const Double_t pitch_mm = pad_pitch_um * 1e-3;
    const Int_t n_pad_recoil = TMath::Nint(Det::recoil_side_mm / pitch_mm);
    const Int_t n_pad_fwd    = TMath::Nint(Det::fwd_side_mm    / pitch_mm);
    const Double_t half_recoil = Det::recoil_side_mm / 2.0;
    const Double_t half_fwd    = Det::fwd_side_mm    / 2.0;

    R.h_recoil_pad = new TH2D(
        TString::Format("h_recoil_pad_%s", R.tag.Data()),
        TString::Format("Recoil pad rate [Hz/pad/stn] (pitch %.0f #mum, %s);"
                        "u [mm];v [mm]", pad_pitch_um,
                        smearing_on ? "smeared" : "nominal"),
        n_pad_recoil, -half_recoil, half_recoil,
        n_pad_recoil, -half_recoil, half_recoil);
    R.h_fwd_pad = new TH2D(
        TString::Format("h_fwd_pad_%s", R.tag.Data()),
        TString::Format("Forward pad rate [Hz/pad/stn] (pitch %.0f #mum, %s);"
                        "u [mm];v [mm]", pad_pitch_um,
                        smearing_on ? "smeared" : "nominal"),
        n_pad_fwd, -half_fwd, half_fwd,
        n_pad_fwd, -half_fwd, half_fwd);

    R.h_recoil_pad_nosm = new TH2D(
        TString::Format("h_recoil_pad_nosm_%s", R.tag.Data()),
        TString::Format("Recoil pad rate (no smearing) [Hz/pad/stn];u [mm];v [mm]"),
        n_pad_recoil, -half_recoil, half_recoil,
        n_pad_recoil, -half_recoil, half_recoil);
    R.h_fwd_pad_nosm = new TH2D(
        TString::Format("h_fwd_pad_nosm_%s", R.tag.Data()),
        TString::Format("Forward pad rate (no smearing) [Hz/pad/stn];u [mm];v [mm]"),
        n_pad_fwd, -half_fwd, half_fwd,
        n_pad_fwd, -half_fwd, half_fwd);

    // Difference histograms: (u,v)_smeared - (u,v)_nominal per event.
    // Range chosen to cover MS-scale excursions on each ring:
    //    recoil: ±2 mm  (MS ~ mm at T_C ~ few MeV)
    //    forward: ±0.5 mm (very small, just for diagnostic)
    R.h_recoil_pad_diff = new TH2D(
        TString::Format("h_recoil_pad_diff_%s", R.tag.Data()),
        "Recoil pad displacement (smeared - nominal);#Delta u [mm];#Delta v [mm]",
        200, -2.0, 2.0, 200, -2.0, 2.0);
    R.h_fwd_pad_diff = new TH2D(
        TString::Format("h_fwd_pad_diff_%s", R.tag.Data()),
        "Forward pad displacement (smeared - nominal);#Delta u [mm];#Delta v [mm]",
        200, -0.5, 0.5, 200, -0.5, 0.5);

    R.h_recoil_pad_charge = new TH2D(
        TString::Format("h_recoil_pad_charge_%s", R.tag.Data()),
        TString::Format("Recoil pad charge density [Hz#timeskeV/pixel/stn] (pitch %.0f #mum);"
                        "u [mm];v [mm]", pad_pitch_um),
        n_pad_recoil, -half_recoil, half_recoil,
        n_pad_recoil, -half_recoil, half_recoil);

    R.h_clusterSize = new TH1D(
        TString::Format("h_clusterSize_%s", R.tag.Data()),
        "Cluster size (pixels above 5% of central charge);N_{pixels};Rate [Hz/stn]",
        20, -0.5, 19.5);

    R.h_clusterSize_vs_TC = new TH2D(
        TString::Format("h_clusterSize_vs_TC_%s", R.tag.Data()),
        "Cluster size vs T_{C};T_{C}^{lab} [MeV];N_{pixels}",
        240, 0, 120, 20, -0.5, 19.5);   
        
    R.h_recoil_charge_vs_TC = new TH2D(
        TString::Format("h_recoil_charge_vs_TC_%s", R.tag.Data()),
        "Recoil charge vs T_{C};T_{C}^{lab} [MeV];Energy deposition [MeV]",
        240, 0, 120, 1000, 0, 10000);    
    R.h_ToF_vs_Edep = new TH2D(
        TString::Format("h_ToF_vs_Edep_%s", R.tag.Data()),
        "Energy deposited vs TOF;E_{dep} [MeV];TOF [ns];Rate [Hz/bin]",
        200, 0, 100 ,  // Edep axis: 0-100 MeV in keV, 500 keV/bin      
        120, 0, 60);  // TOF axis: 0-60 ns, 0.5 ns/bin 
    // -------- MC loop --------
    TRandom3 rng(12345);
    Double_t sumw_recoil_2pi = 0, sumw_fwd_2pi = 0;
    Double_t sumw_recoil_pad = 0, sumw_fwd_pad = 0;   // with smearing
    Double_t sumw_recoil_pad_nosm = 0, sumw_fwd_pad_nosm = 0;

    for (Int_t iev = 0; iev < n_events; ++iev) {
        Double_t costh = -1.0 + 2.0 * rng.Rndm();
        Double_t phi   = 2.0 * TMath::Pi() * rng.Rndm();
        Double_t sinth = TMath::Sqrt(TMath::Max(0.0, 1.0 - costh*costh));

        Double_t px_cm = k.p_cm * sinth * TMath::Cos(phi);
        Double_t py_cm = k.p_cm * sinth * TMath::Sin(phi);
        Double_t pz_cm = k.p_cm * costh;
        Double_t pT    = TMath::Sqrt(px_cm*px_cm + py_cm*py_cm);
        Double_t pz_p_lab = k.gamma_cm * (pz_cm + k.beta_cm * k.E1_cm);
        Double_t p_p_lab  = TMath::Sqrt(pT*pT + pz_p_lab*pz_p_lab);
        Double_t th_p_lab = TMath::RadToDeg() * TMath::ATan2(pT, pz_p_lab);

        Double_t J_lab_over_cm =
            k.gamma_cm * k.p_cm * k.p_cm
            * (k.p_cm + k.beta_cm * k.E1_cm * costh)
            / (p_p_lab*p_p_lab*p_p_lab);
        if (J_lab_over_cm <= 0) continue;

        Double_t dsigdO_mb = 0;
        if (use_bonin)       dsigdO_mb = Bonin::dSig_mb(th_p_lab, R.p_lab_GeVc);
        else if (use_kamiya) dsigdO_mb = Kamiya::dSig_mb(th_p_lab);
        if (dsigdO_mb <= 0) continue;

        Double_t w = dsigdO_mb * J_lab_over_cm;

        // Recoil in CM (opposite projectile)
        Double_t px_C = -px_cm;
        Double_t py_C = -py_cm;
        Double_t pz_C_cm = -pz_cm;

        // Boost to lab
        Double_t pz_C_lab = k.gamma_cm * (pz_C_cm + k.beta_cm * k.E2_cm);
        Double_t E_C_lab  = k.gamma_cm * (k.E2_cm  + k.beta_cm * pz_C_cm);
        if (pz_C_lab <= 0) continue;
        Double_t pT_C   = TMath::Sqrt(px_C*px_C + py_C*py_C);
        Double_t p_C    = TMath::Sqrt(pT_C*pT_C + pz_C_lab*pz_C_lab);
        Double_t th_C_lab = TMath::RadToDeg() * TMath::ATan2(pT_C, pz_C_lab);
        Double_t T_C_lab  = E_C_lab - PC::mC;
        Double_t beta_C   = (p_C > 0.0) ? (p_C / E_C_lab) : 0.0;

        // TOF (nominal, pre-smearing)
        Double_t R_stn_cm = TMath::Sqrt(Det::recoil_r_mm*Det::recoil_r_mm
                                        + Det::recoil_z_mm*Det::recoil_z_mm)
                            * 0.1;
        Double_t path_cm  = R_stn_cm / TMath::Sin(TMath::DegToRad()*th_C_lab);
        Double_t tof_ns   = (beta_C > 0) ? path_cm / (beta_C * PC::c_cms) : 0;

        // Unit directions (nominal)
        Double_t dxp = px_cm/p_p_lab, dyp = py_cm/p_p_lab, dzp = pz_p_lab/p_p_lab;
        Double_t dxC_nom = px_C/p_C,   dyC_nom = py_C/p_C,   dzC_nom = pz_C_lab/p_C;

        // Kinematic hists
        R.h_theta_C->Fill(th_C_lab, w);
        R.h_TC     ->Fill(T_C_lab,  w);
        R.h_theta_p->Fill(th_p_lab, w);
        R.h_tof    ->Fill(tof_ns,   w);
        R.h_TC_th  ->Fill(th_C_lab, T_C_lab, w);

        if (th_C_lab >= Det::recoil_thC_lo && th_C_lab <= Det::recoil_thC_hi)
            sumw_recoil_2pi += w;
        if (th_p_lab >= Det::fwd_theta_lo && th_p_lab <= Det::fwd_theta_hi)
            sumw_fwd_2pi += w;

        Double_t phi_C = TMath::ATan2(dyC_nom, dxC_nom);
        Double_t phi_p = TMath::ATan2(dyp, dxp);
        Double_t phi_stn_C = Det::nearest_stn_phi(phi_C);
        Double_t phi_stn_p = Det::nearest_stn_phi(phi_p);

        // ---- Nominal (no smearing) hits ----
        Double_t u_nom, v_nom;
        Bool_t hit_r_nom = HitRecoilPad(dxC_nom, dyC_nom, dzC_nom,
                                         phi_stn_C, u_nom, v_nom);
        if (hit_r_nom) {
            sumw_recoil_pad_nosm += w;
            R.h_recoil_pad_nosm->Fill(u_nom, v_nom, w);
        }

        Double_t uf_nom, vf_nom;
        Bool_t hit_f_nom = HitFwdPad(dxp, dyp, dzp, phi_stn_p, uf_nom, vf_nom);
        if (hit_f_nom) {
            sumw_fwd_pad_nosm += w;
            R.h_fwd_pad_nosm->Fill(uf_nom, vf_nom, w);
        }

        // ---- Smeared hits (same MC event) ----
        Double_t dxC = dxC_nom, dyC = dyC_nom, dzC = dzC_nom;
        Double_t T_C_s = T_C_lab, E_C_s = E_C_lab, p_C_s = p_C, beta_C_s = beta_C;
        TVector3 origin_exit(0, 0, 0);
        ApplyTargetSmearing(smear, rng, R.bt.deff_nm,
                            dxC, dyC, dzC,
                            T_C_s, E_C_s, p_C_s, beta_C_s,
                            origin_exit);

        Double_t u_s, v_s;
        Bool_t hit_r_s = HitRecoilPad(dxC, dyC, dzC,
                                       phi_stn_C, u_s, v_s, origin_exit);
        if (hit_r_s || hit_r_nom) {
            sumw_recoil_pad += w;
            R.h_recoil_pad->Fill(u_s, v_s, w);
            if (hit_r_nom) {
                R.h_recoil_pad_diff->Fill(u_s - u_nom, v_s - v_nom, w);
                R.h_TC_du ->Fill(T_C_lab, u_s - u_nom, w);
                R.h_TC_dv ->Fill(T_C_lab, v_s - v_nom, w);
            }
        }

        // Forward ring: projectile traverses same 800 nm.  MS ~27 urad @ 500 MeV
        // at r=750 mm gives ~20 um smear. Skipping for speed — origin stays at IP.
        Double_t uf_s, vf_s;
        Bool_t hit_f_s = HitFwdPad(dxp, dyp, dzp, phi_stn_p, uf_s, vf_s);
        if (hit_f_s || hit_f_nom) {
            sumw_fwd_pad += w;
            R.h_fwd_pad->Fill(uf_s, vf_s, w);
            if (hit_f_nom) {
                R.h_fwd_pad_diff->Fill(uf_s - uf_nom, vf_s - vf_nom, w);
            }
        }


        // --- Charge deposition + Gaussian diffusion (recoil pad only) ---
        if (hit_r_s) {
            // Energy actually deposited in 500 um Si
            Double_t Edep_MeV = CinSi::Edep_MeV(T_C_s, Si_thickness_um);
            Double_t Edep_keV = Edep_MeV * 1000.0;
            
            R.h_ToF_vs_Edep->Fill(Edep_MeV, tof_ns, w);
            R.h_recoil_charge_vs_TC->Fill(T_C_s, Edep_keV, w);
            // Quick cluster-size estimate: count how many neighbour pixels
            // would receive > 5% of the central pixel's charge for this single hit.
            // Using analytic Gaussian integral over each pixel.
            Int_t pix_u0 = (Int_t) TMath::Floor((u_s + half_recoil) / pitch_mm);
            Int_t pix_v0 = (Int_t) TMath::Floor((v_s + half_recoil) / pitch_mm);
            Int_t cluster_n = 0;
            Double_t cluster_E_keV = 0.0;


            Double_t s = sigma_diff_mm * TMath::Sqrt(2.0);

            for (Int_t du = -pixel_window; du <= pixel_window; ++du) {
                for (Int_t dv = -pixel_window; dv <= pixel_window; ++dv) {
                    Int_t pu = pix_u0 + du;
                    Int_t pv = pix_v0 + dv;
                    if (pu < 0 || pv < 0)               continue;
                    if (pu >= n_pad_recoil)             continue;
                    if (pv >= n_pad_recoil)             continue;
                
                    Double_t lo_u = -half_recoil + pu     * pitch_mm;
                    Double_t hi_u = -half_recoil + (pu+1) * pitch_mm;
                    Double_t lo_v = -half_recoil + pv     * pitch_mm;
                    Double_t hi_v = -half_recoil + (pv+1) * pitch_mm;
                
                    Double_t fu = 0.5*(TMath::Erf((hi_u - u_s)/s) - TMath::Erf((lo_u - u_s)/s));
                    Double_t fv = 0.5*(TMath::Erf((hi_v - v_s)/s) - TMath::Erf((lo_v - v_s)/s));
                    Double_t frac    = fu * fv;
                    Double_t q_pixel = Edep_keV * frac;
                
                    if (q_pixel > Q_thr_keV) {
                        Double_t pix_u_center = -half_recoil + (pu + 0.5) * pitch_mm;
                        Double_t pix_v_center = -half_recoil + (pv + 0.5) * pitch_mm;
                    
                        R.h_recoil_pad_charge->Fill(pix_u_center, pix_v_center, w);// * q_pixel);
                    
                        cluster_n++;
                        cluster_E_keV += q_pixel;
                    }
                }
            }
            // Per-event cluster statistics (rate-weighted)
            if (cluster_n > 0) {
                R.h_clusterSize       ->Fill(cluster_n,            w);
                R.h_clusterSize_vs_TC ->Fill(T_C_s, cluster_n,     w);
            }
        }
    }

    // -------- Scale to rate --------
    const Double_t ps_vol = 4.0 * TMath::Pi();
    R.sigma_recoil_2pi_mb = (ps_vol / n_events) * sumw_recoil_2pi;
    R.sigma_fwd_2pi_mb    = (ps_vol / n_events) * sumw_fwd_2pi;
    R.sigma_recoil_pad_mb = (ps_vol / n_events) * sumw_recoil_pad;
    R.sigma_fwd_pad_mb    = (ps_vol / n_events) * sumw_fwd_pad;

    const Double_t flux = R.bt.N_cyc * R.bt.f_rev * R.bt.n_t * R.bt.f_ov;
    R.rate_recoil_2pi = flux * R.sigma_recoil_2pi_mb * PC::mb2cm2;
    R.rate_fwd_2pi    = flux * R.sigma_fwd_2pi_mb    * PC::mb2cm2;
    R.rate_recoil_pad = flux * R.sigma_recoil_pad_mb * PC::mb2cm2;
    R.rate_fwd_pad    = flux * R.sigma_fwd_pad_mb    * PC::mb2cm2;

    Double_t Om_recoil_2pi = 2*TMath::Pi() *
        (TMath::Cos(TMath::DegToRad()*Det::recoil_thC_lo)
         - TMath::Cos(TMath::DegToRad()*Det::recoil_thC_hi));
    Double_t Om_fwd_2pi = 2*TMath::Pi() *
        (TMath::Cos(TMath::DegToRad()*Det::fwd_theta_lo)
         - TMath::Cos(TMath::DegToRad()*Det::fwd_theta_hi));
    R.dsigdO_avg_recoil = (Om_recoil_2pi > 0)
        ? R.sigma_recoil_2pi_mb / Om_recoil_2pi : 0;
    R.dsigdO_avg_fwd    = (Om_fwd_2pi > 0)
        ? R.sigma_fwd_2pi_mb / Om_fwd_2pi : 0;

    const Double_t scale_counts_to_rate = flux * (ps_vol / n_events) * PC::mb2cm2;
    auto scale1D = [&](TH1D* h){
        h->Scale(scale_counts_to_rate / h->GetBinWidth(1));
    };
    scale1D(R.h_theta_C);
    scale1D(R.h_TC);
    scale1D(R.h_theta_p);
    scale1D(R.h_tof);
    R.h_TC_th->Scale(scale_counts_to_rate);

    R.h_recoil_pad     ->Scale(scale_counts_to_rate / Det::n_stn);
    R.h_fwd_pad        ->Scale(scale_counts_to_rate / Det::n_stn);
    R.h_recoil_pad_nosm->Scale(scale_counts_to_rate / Det::n_stn);
    R.h_fwd_pad_nosm   ->Scale(scale_counts_to_rate / Det::n_stn);
    R.h_recoil_pad_diff->Scale(scale_counts_to_rate / Det::n_stn);
    R.h_fwd_pad_diff   ->Scale(scale_counts_to_rate / Det::n_stn);
    R.h_TC_dv          ->Scale(scale_counts_to_rate / Det::n_stn);
    R.h_TC_du          ->Scale(scale_counts_to_rate / Det::n_stn);

    R.h_recoil_pad_charge->Scale(scale_counts_to_rate / Det::n_stn);
    R.h_clusterSize      ->Scale(scale_counts_to_rate / Det::n_stn);
    R.h_clusterSize_vs_TC->Scale(scale_counts_to_rate / Det::n_stn);  
    R.h_recoil_charge_vs_TC->Scale(scale_counts_to_rate / Det::n_stn);
    R.h_ToF_vs_Edep->Scale(scale_counts_to_rate);
     
    if (verbose) {
        Printf("\n   --- 2π-az θ-only ring (ideal) ---");
        Printf("     sigma(recoil,2π) = %.3f mb    (Ω = %.3f sr)",
               R.sigma_recoil_2pi_mb, Om_recoil_2pi);
        Printf("     <dσ/dΩ>_recoil   = %.1f mb/sr", R.dsigdO_avg_recoil);
        Printf("     sigma(fwd,   2π) = %.3f mb    (Ω = %.3f sr)",
               R.sigma_fwd_2pi_mb, Om_fwd_2pi);
        Printf("     <dσ/dΩ>_fwd      = %.1f mb/sr", R.dsigdO_avg_fwd);

        Printf("\n   --- Physical 6-pad acceptance (%s) ---",
               smearing_on ? "with smearing" : "no smearing");
        Printf("     sigma_recoil_pads = %.3f mb  (6 stations)",
               R.sigma_recoil_pad_mb);
        Printf("     sigma_fwd_pads    = %.3f mb  (6 stations)",
               R.sigma_fwd_pad_mb);

        Printf("\n   --- Rates ---");
        Printf("     R_recoil (2π az)  = %.3e Hz    ÷6 = %.3e Hz",
               R.rate_recoil_2pi, R.rate_recoil_2pi/Det::n_stn);
        Printf("     R_recoil (6 pads) = %.3e Hz    ÷6 = %.3e Hz   <<  physical",
               R.rate_recoil_pad, R.rate_recoil_pad/Det::n_stn);
        Printf("     R_fwd    (2π az)  = %.3e Hz    ÷6 = %.3e Hz",
               R.rate_fwd_2pi, R.rate_fwd_2pi/Det::n_stn);
        Printf("     R_fwd    (6 pads) = %.3e Hz    ÷6 = %.3e Hz   <<  physical",
               R.rate_fwd_pad, R.rate_fwd_pad/Det::n_stn);

        Double_t drms_u = R.h_recoil_pad_diff->GetRMS(1);
        Double_t drms_v = R.h_recoil_pad_diff->GetRMS(2);
        Printf("\n   --- Smearing diagnostics (recoil pad) ---");
        Printf("     RMS(Δu) = %.3f mm,  RMS(Δv) = %.3f mm",
               drms_u, drms_v);

        Printf("\n   --- Charge / cluster summary ---");
        Printf("     Mean cluster size      = %.2f pixels (rate-weighted)",
               R.h_clusterSize->GetMean());
        Printf("     Total deposited charge = %.3e keV/s/stn",
               R.h_recoil_pad_charge->Integral());
    }

    return R;
}

// =============================================================================
//  Plots
// =============================================================================
static void StyleH1(TH1* h, Int_t col) {
    h->SetLineColor(col);
    h->SetLineWidth(2);
    h->SetFillColorAlpha(col, 0.15);
}

void MakeKinPlot(const SimResult& R) {
    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas(
        TString::Format("c_kin_%s", R.tag.Data()),
        TString::Format("Kinematics %s", R.tag.Data()),
        1400, 900);
    c->Divide(2, 2, 0.004, 0.004);

    TString spc_name = (R.species == "p") ? "pC" : "^{3}HeC";

    c->cd(1); gPad->SetGrid();
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    StyleH1(R.h_theta_C, kRed+1);
    R.h_theta_C->Draw("HIST");
    TBox* b1 = new TBox(Det::recoil_thC_lo, 0,
                        Det::recoil_thC_hi, R.h_theta_C->GetMaximum());
    b1->SetFillColorAlpha(kOrange, 0.15);
    b1->SetLineColor(kOrange+2);
    b1->Draw();
    R.h_theta_C->Draw("HIST SAME");
    R.h_theta_C->SetTitle(TString::Format(
        "%s, T=%.0f MeV: rate vs #theta_{C}^{lab};"
        "#theta_{C}^{lab} [deg];Rate [counts/s/deg]",
        spc_name.Data(), R.T_MeV));

    c->cd(2); gPad->SetGrid();
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    StyleH1(R.h_TC, kBlue+1);
    R.h_TC->Draw("HIST");
    R.h_TC->SetTitle(TString::Format(
        "%s, T=%.0f MeV: rate vs T_{C}^{lab};"
        "T_{C}^{lab} [MeV];Rate [counts/s/MeV]",
        spc_name.Data(), R.T_MeV));

    c->cd(3); gPad->SetGrid();
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    StyleH1(R.h_theta_p, kGreen+2);
    R.h_theta_p->Draw("HIST");
    TBox* b3 = new TBox(Det::fwd_theta_lo, 0,
                        Det::fwd_theta_hi, R.h_theta_p->GetMaximum());
    b3->SetFillColorAlpha(kOrange, 0.15);
    b3->SetLineColor(kOrange+2);
    b3->Draw();
    R.h_theta_p->Draw("HIST SAME");
    R.h_theta_p->SetTitle(TString::Format(
        "%s, T=%.0f MeV: rate vs #theta_{proj}^{lab};"
        "#theta_{proj}^{lab} [deg];Rate [counts/s/deg]",
        spc_name.Data(), R.T_MeV));

    c->cd(4); gPad->SetGrid();
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.14);
    R.h_TC_th->Draw("COLZ");
    R.h_TC_th->SetTitle(TString::Format(
        "%s, T=%.0f MeV: elastic locus;"
        "#theta_{C}^{lab} [deg];T_{C}^{lab} [MeV]",
        spc_name.Data(), R.T_MeV));

    c->cd(1);
    TPaveText* pt = new TPaveText(0.16, 0.66, 0.62, 0.92, "NDC");
    pt->SetFillStyle(1001); pt->SetFillColor(kWhite);
    pt->SetBorderSize(1); pt->SetTextSize(0.030); pt->SetTextAlign(12);
    pt->AddText(TString::Format("R_{recoil, 2#pi}   = %.2e /s",
                                 R.rate_recoil_2pi));
    pt->AddText(TString::Format("R_{recoil, pads}  = %.2e /s   (%.2e /stn)",
                                 R.rate_recoil_pad,
                                 R.rate_recoil_pad/Det::n_stn));
    pt->AddText(TString::Format("R_{fwd,    2#pi}  = %.2e /s",
                                 R.rate_fwd_2pi));
    pt->AddText(TString::Format("R_{fwd,    pads} = %.2e /s   (%.2e /stn)",
                                 R.rate_fwd_pad,
                                 R.rate_fwd_pad/Det::n_stn));
    pt->AddText(TString::Format("<d#sigma/d#Omega>_{recoil} = %.1f mb/sr",
                                 R.dsigdO_avg_recoil));
    pt->AddText(TString::Format("<d#sigma/d#Omega>_{fwd}    = %.1f mb/sr",
                                 R.dsigdO_avg_fwd));
    pt->Draw();

    TString out = TString::Format("kin_%s", R.tag.Data());
    c->SaveAs(out + ".png");
    c->SaveAs(out + ".pdf");
}

void MakePadPlot(const SimResult& R) {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(64);

    TCanvas* c = new TCanvas(
        TString::Format("c_pads_%s", R.tag.Data()),
        TString::Format("Pad maps %s", R.tag.Data()),
        1400, 1300);
    c->Divide(2, 2, 0.004, 0.004);

    TString spc_name = (R.species == "p") ? "pC" : "^{3}HeC";

    // Row 1: smeared pad maps
    c->cd(1);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.16);
    gPad->SetBottomMargin(0.13);
    R.h_recoil_pad->Draw("COLZ");
    R.h_recoil_pad->SetTitle(TString::Format(
        "%s, T=%.0f MeV: recoil pad (smeared) [Hz/pad/stn];"
        "u [mm];v [mm]", spc_name.Data(), R.T_MeV));
    TLatex lt;
    lt.SetTextSize(0.030); lt.SetNDC();
    lt.DrawLatex(0.16, 0.85,
        TString::Format("Total/stn: %.2e Hz   Peak: %.2e Hz/pad",
                        R.h_recoil_pad->Integral(),
                        R.h_recoil_pad->GetMaximum()));

    c->cd(2);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.16);
    gPad->SetBottomMargin(0.13);
    R.h_fwd_pad->Draw("COLZ");
    R.h_fwd_pad->SetTitle(TString::Format(
        "%s, T=%.0f MeV: forward pad (smeared) [Hz/pad/stn];"
        "u [mm];v [mm]", spc_name.Data(), R.T_MeV));
    lt.DrawLatex(0.16, 0.85,
        TString::Format("Total/stn: %.2e Hz   Peak: %.2e Hz/pad",
                        R.h_fwd_pad->Integral(),
                        R.h_fwd_pad->GetMaximum()));

    // Row 2: displacement distributions (smeared - nominal)
    c->cd(3);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.16);
    gPad->SetBottomMargin(0.13);
    R.h_recoil_pad_diff->Draw("COLZ");
    lt.DrawLatex(0.16, 0.85,
        TString::Format("RMS(#Deltau)=%.3f mm, RMS(#Deltav)=%.3f mm",
                        R.h_recoil_pad_diff->GetRMS(1),
                        R.h_recoil_pad_diff->GetRMS(2)));

    c->cd(4);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.16);
    gPad->SetBottomMargin(0.13);
    R.h_fwd_pad_diff->Draw("COLZ");
    lt.DrawLatex(0.16, 0.85,
        TString::Format("RMS(#Deltau)=%.3f mm, RMS(#Deltav)=%.3f mm",
                        R.h_fwd_pad_diff->GetRMS(1),
                        R.h_fwd_pad_diff->GetRMS(2)));

    TString out = TString::Format("pads_%s", R.tag.Data());
    c->SaveAs(out + ".png");
    c->SaveAs(out + ".pdf");
}

void MakeDiffPlot(const SimResult& R) {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(64);

    TCanvas* c = new TCanvas(
        TString::Format("c_diff_%s", R.tag.Data()),
        TString::Format("Displacement distributions %s", R.tag.Data()),
        1400, 1300);
    c->Divide(2, 2, 0.004, 0.004);
    TString spc_name = (R.species == "p") ? "pC" : "^{3}HeC";
    c->cd(1); gPad->SetGrid();
    TH1D* h_recoil_pad_x = (TH1D*)R.h_recoil_pad->ProjectionX();
    TH1D* h_recoil_pad_y = (TH1D*)R.h_recoil_pad->ProjectionY();
    TH1D* h_recoil_pad_nosm_x = (TH1D*)R.h_recoil_pad_nosm->ProjectionX();
    TH1D* h_recoil_pad_nosm_y = (TH1D*)R.h_recoil_pad_nosm->ProjectionY();

    h_recoil_pad_x->Rebin(4);  // rebin for smoother comparison
    h_recoil_pad_y->Rebin(4);  // rebin for smoother comparison
    h_recoil_pad_nosm_x->Rebin(4);  // rebin for smoother comparison
    h_recoil_pad_nosm_y->Rebin(4);  // rebin for smoother comparison


    h_recoil_pad_x->Scale(1./4);  // compensate for rebinning);
    h_recoil_pad_y->Scale(1./4);  // compensate for rebinning);
    h_recoil_pad_nosm_x->Scale(1./4);  // compensate for rebinning);
    h_recoil_pad_nosm_y->Scale(1./4);  // compensate for rebinning);

    h_recoil_pad_nosm_x->SetLineColor(kRed+1);
    h_recoil_pad_nosm_y->SetLineColor(kRed+1);

    h_recoil_pad_nosm_x->GetYaxis()->SetRangeUser(0, TMath::Max(h_recoil_pad_nosm_x->GetMaximum(), h_recoil_pad_x->GetMaximum())*1.5);

    h_recoil_pad_nosm_x->Draw("HIST");
    h_recoil_pad_x->Draw("HISTsame");

    TLegend* leg1 = new TLegend(0.16, 0.66, 0.62, 0.92);
    leg1->AddEntry(R.h_recoil_pad, "smeared", "l");
    leg1->AddEntry(h_recoil_pad_nosm_x, "nominal (no smear)", "l");
    leg1->SetTextSize(0.048);
    leg1-> SetTextFont(42);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->Draw();

    c->cd(2); gPad->SetGrid();

    h_recoil_pad_nosm_y->Draw("HIST");
    h_recoil_pad_y->Draw("HISTsame");
    TLegend* leg2 = new TLegend(0.16, 0.66, 0.62, 0.92);
    leg2->AddEntry(R.h_recoil_pad, "smeared", "l");
    leg2->AddEntry(h_recoil_pad_nosm_y, "nominal (no smear)", "l");
    leg2->SetTextSize(0.048);
    leg2-> SetTextFont(42);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->Draw();

    c->cd(3);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.16);
    gPad->SetBottomMargin(0.13);
    R.h_TC_dv->Draw("COLZ");
    c->cd(4);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.16);
    gPad->SetBottomMargin(0.13);
    R.h_TC_du->Draw("COLZ");

    TString out = TString::Format("diff_%s", R.tag.Data());
    c->SaveAs(out + ".png");
    c->SaveAs(out + ".pdf");
}
void MakeChargePlot(const SimResult& R) {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(64);

    TCanvas* c = new TCanvas(
        TString::Format("c_charge_%s", R.tag.Data()),
        TString::Format("Charge distributions %s", R.tag.Data()),
        1400, 2000);
    c->Divide(2, 3, 0.004, 0.004);
    TString spc_name = (R.species == "p") ? "pC" : "^{3}HeC";
    c->cd(1); gPad->SetGrid();
    TH1D* h_recoil_pad_x = (TH1D*)R.h_recoil_pad->ProjectionX();
    TH1D* h_recoil_pad_y = (TH1D*)R.h_recoil_pad->ProjectionY();
    TH1D* h_recoil_pad_charge_x = (TH1D*)R.h_recoil_pad_charge->ProjectionX();
    TH1D* h_recoil_pad_charge_y = (TH1D*)R.h_recoil_pad_charge->ProjectionY();

    h_recoil_pad_x->Rebin(4);  // rebin for smoother comparison
    h_recoil_pad_y->Rebin(4);  // rebin for smoother comparison
    h_recoil_pad_charge_x->Rebin(4);  // rebin for smoother comparison
    h_recoil_pad_charge_y->Rebin(4);  // rebin for smoother comparison


    h_recoil_pad_x->Scale(1./4);  // compensate for rebinning);
    h_recoil_pad_y->Scale(1./4);  // compensate for rebinning);
    h_recoil_pad_charge_x->Scale(1./4);  // compensate for rebinning);
    h_recoil_pad_charge_y->Scale(1./4);  // compensate for rebinning);

    h_recoil_pad_charge_x->SetLineColor(kRed+1);
    h_recoil_pad_charge_y->SetLineColor(kRed+1);

    h_recoil_pad_charge_x->GetYaxis()->SetRangeUser(0, TMath::Max(h_recoil_pad_charge_x->GetMaximum(), h_recoil_pad_x->GetMaximum())*1.5);

    h_recoil_pad_charge_x->Draw("HIST");
    h_recoil_pad_x->Draw("HISTsame");

    TLegend* leg1 = new TLegend(0.16, 0.66, 0.62, 0.92);
    leg1->SetBorderSize(1); leg1->SetFillColor(kWhite);
    leg1->SetTextSize(0.030); leg1->SetTextAlign(12);
    leg1->AddEntry(R.h_recoil_pad, "smeared", "l");
    leg1->AddEntry(h_recoil_pad_charge_x, "charge", "l");
    leg1->SetTextSize(0.048);
    leg1-> SetTextFont(42);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->Draw();

    c->cd(2); gPad->SetGrid();

    h_recoil_pad_charge_y->Draw("HIST");
    h_recoil_pad_y->Draw("HISTsame");
    TLegend* leg2 = new TLegend(0.16, 0.66, 0.62, 0.92);
    leg2->SetBorderSize(1); leg2->SetFillColor(kWhite);
    leg2->SetTextSize(0.030); leg2->SetTextAlign(12);
    leg2->AddEntry(R.h_recoil_pad, "smeared", "l");
    leg2->AddEntry(h_recoil_pad_charge_y, "charge", "l");
    leg2->SetTextSize(0.048);
    leg2-> SetTextFont(42);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->Draw();

    c->cd(3);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.16);
    gPad->SetBottomMargin(0.13);
    R.h_ToF_vs_Edep->Draw("COLZ");
    c->cd(4);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.16);
    gPad->SetBottomMargin(0.13);
    R.h_recoil_pad_charge->Draw("COLZ");
        TLatex lt;
    lt.SetTextSize(0.030); lt.SetNDC();
    lt.DrawLatex(0.16, 0.85,
    TString::Format("Total/stn: %.2e Hz   Peak: %.2e Hz/pad",
                    R.h_recoil_pad_charge->Integral(),
                    R.h_recoil_pad_charge->GetMaximum()));
    
    c->cd(5);
    R.h_clusterSize       ->Draw("HIST");
    
    c->cd(6);            
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.16);
    gPad->SetBottomMargin(0.13);
    R.h_clusterSize_vs_TC ->Draw("COLZ");
    
    TString out = TString::Format("charge_%s", R.tag.Data());
    c->SaveAs(out + ".png");
    c->SaveAs(out + ".pdf");
}
// =============================================================================
//  Wrappers
// =============================================================================
SimResult pC_rate_sim(Double_t T_MeV = 500.0, Int_t n_events = 500000,
                      Double_t pad_pitch_um = 55.0,
                      Bool_t smearing_on = kTRUE) {
    SimResult R = RunSim("p", T_MeV, n_events, pad_pitch_um, smearing_on);
    MakeKinPlot(R);
    MakePadPlot(R);
    MakeDiffPlot(R);
    MakeChargePlot(R);
    return R;
}

SimResult He3C_rate_sim(Double_t T_MeV = 443.0, Int_t n_events = 500000,
                         Double_t pad_pitch_um = 55.0,
                         Bool_t smearing_on = kTRUE) {
    SimResult R = RunSim("3He", T_MeV, n_events, pad_pitch_um, smearing_on);
    MakeKinPlot(R);
    MakePadPlot(R);
    MakeDiffPlot(R);

    return R;
}

void run_all(Int_t n_events = 5000000, Double_t pad_pitch_um = 55.0,
             Bool_t smearing_on = kTRUE) {
    std::vector<SimResult> Rs;
    Rs.push_back(pC_rate_sim( 500.0, n_events, pad_pitch_um, smearing_on));
    Rs.push_back(pC_rate_sim( 850.0, n_events, pad_pitch_um, smearing_on));
    Rs.push_back(pC_rate_sim(1200.0, n_events, pad_pitch_um, smearing_on));
    //Rs.push_back(He3C_rate_sim(443.0, n_events, pad_pitch_um, smearing_on));

    Printf("\n\n===========================================================================================");
    Printf("  SUMMARY vs paper Table 13     [smearing: %s]",
           smearing_on ? "ON" : "OFF");
    Printf("===========================================================================================");
    Printf("  Species   T    <dσ/dΩ>      R_2π       R/stn_2π    R_pads     R/stn_pad");
    Printf("          [MeV]  [mb/sr]      [Hz]        [Hz]        [Hz]       [Hz]");
    Printf("-------------------------------------------------------------------------------------------");
    for (const auto& r : Rs) {
        TString spc = (r.species == "p") ? "pC" : "3HeC";
        Bool_t is_p = (r.species == "p");
        Double_t avg  = is_p ? r.dsigdO_avg_recoil : r.dsigdO_avg_fwd;
        Double_t r2p  = is_p ? r.rate_recoil_2pi   : r.rate_fwd_2pi;
        Double_t rpad = is_p ? r.rate_recoil_pad   : r.rate_fwd_pad;
        Printf("  %-6s %5.0f  %9.1f   %.3e   %.3e   %.3e   %.3e",
               spc.Data(), r.T_MeV, avg, r2p, r2p/Det::n_stn,
               rpad, rpad/Det::n_stn);
    }
    Printf("-------------------------------------------------------------------------------------------");
    Printf("  Paper Table 13 (per station):");
    Printf("    pC   500 MeV: <dσ/dΩ>=2470 mb/sr   R/stn = 1.6 MHz");
    Printf("    pC   850 MeV: <dσ/dΩ>=3820 mb/sr   R/stn = 3.3 MHz");
    Printf("    pC  1200 MeV: <dσ/dΩ>=5036 mb/sr   R/stn = 5.1 MHz");
    Printf("    3HeC 443 MeV: <dσ/dΩ>=  30 mb/sr   R/stn = 635 Hz   (forward ring, coincidence)");
    Printf("===========================================================================================");
}