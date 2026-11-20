"""Revision studies (addressing referee-style comments). Each study writes data_v3/<name>.csv.
Not meant to be called directly: use run_revision.py, which runs one DEVSIM device per process.

  uq_gated          Monte Carlo of the gated on state at the literature alignment (S>=0.05, dEV sigma 0.3 eV)
  thermal_lit       self-heating at the measured / calculated pinning (Pt), incl. contact sink and larger R_TBR
  mesh_convergence  mesh-factor 1,2,3 for the tunnelling-dominated gated cases and the Table 3 operating points
  thickness_lit     layer-number study at the literature alignment, and valence-band-share sensitivity
  inse_electron     electron injection into InSe (measured S=0.05) versus T and V_G, for comparison with Ref. [11]
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import studies as s  # noqa: E402
import tcad2d as m  # noqa: E402

NM = s.NM
GATE = s.GATE
PT = s.METALS["Pt"]


def _cnl_cases(t=5 * NM):
    """Pinning cases defined relative to the band edges of a film of thickness t (same rule as studies.pin_cases)."""
    P0 = dict(m.DEFAULTS); P0.update(s.LIT)
    Eg, Ec, _ = m.band(P0, t)
    Ev = Ec + Eg
    return {"direct (measured)": dict(S=0.05, CNL=Ec + 0.30),
            "ideal MIGS (DFT)": dict(S=0.18, CNL=Ev - 0.58),
            "interlayer (measured)": dict(S=0.50, CNL=Ec + 0.60),
            "unpinned": dict(S=1.0, CNL=Ec + 0.30)}


# ----------------------------------------------------------------------------- worker
def run_dev(spec):
    kind = spec["kind"]
    tag = spec["tag"]
    try:
        gate = spec.get("gate_dict", GATE) if spec.get("gate") else None
        kw = dict(Lch=spec["Lch"]) if "Lch" in spec else {}
        dv = m.Device(kind, P=spec["P"], Phi_M=spec["Phi"], Lc=spec["Lc"], t=spec.get("t", 5 * NM),
                      gate=gate, mf=spec.get("mf", 1.0), **kw)
        rows = []
        for st in spec["steps"]:
            if "VG" in st:
                if "vg_step" in spec:  # thick oxide: larger gate steps (the ramp still halves on failure)
                    dv._ramp("gate_bias", "vg", st["VG"], spec["vg_step"], max_step=2 * spec["vg_step"])
                else:
                    dv.set_gate(st["VG"])
            if "V" in st:
                dv.set_bias(st["V"])
            if "heat_step" in st:
                arr = dv.heat_step(st["heat_step"])
                rows += [dict(tag, geom=s.lab(kind, spec["Lc"]), Phi_M=spec["Phi"], t=a, Tmax=b, P_el=c_, P_out=d_)
                         for a, b, c_, d_ in arr]
                continue
            if st.get("save", True):  # a step with save=False only moves the operating point
                rows.append(s.summarize(dv, dict(tag, **st.get("tag", {}))))
        return spec["out"], rows
    except Exception as e:  # keep going; failures are recorded, not hidden
        import traceback
        print("RES FAILED", spec["out"], kind, tag, str(e)[:200], flush=True)
        print("TRACEBACK", traceback.format_exc()[-2000:], flush=True)
        return spec["out"], [dict(tag, geom=s.lab(kind, spec["Lc"]), Phi_M=spec["Phi"], failed=1, err=str(e)[:120])]


def _spec(out, kind, Phi, P, tag, steps, gate=True, t=5 * NM, mf=1.0):
    return dict(out=out, kind=kind, Phi=Phi, Lc=50 * NM, P=P, tag=tag, steps=steps, gate=gate, t=t, mf=mf)


GEOMS = ["top", "edge"]


# ----------------------------------------------------------------------------- tasks
def tasks_uq_gated():
    N = int(os.environ.get("UQ_N", "240"))
    rng = np.random.default_rng(7)
    T = []
    for k in range(N):
        lg = lambda a, b: 10 ** rng.uniform(math.log10(a), math.log10(b))
        smp = dict(sample=k, dEV=rng.normal(0, 0.3), S=lg(0.05, 0.5), S_edge_indep=lg(0.05, 0.5),
                   cnl_off=rng.uniform(0.1, 0.6), d_vdW=rng.uniform(0.2, 0.4), Phi_vdW=rng.uniform(1, 3),
                   mu_p=lg(20, 100), rmu=lg(0.002, 0.1), NA=lg(1e15, 1e17), eps_r=rng.uniform(5, 8),
                   mh=rng.uniform(0.1, 0.3))
        P = dict(s.LIT, chi_bulk=s.LIT["chi_bulk"] - smp["dEV"], d_vdW=smp["d_vdW"] * NM,
                 Phi_vdW=smp["Phi_vdW"], mu_p=smp["mu_p"], mu_perp_ratio=smp["rmu"], NA=smp["NA"],
                 eps_r=smp["eps_r"], mh_dos=smp["mh"], mh_tun=smp["mh"], thermal=False)
        P0 = dict(m.DEFAULTS); P0.update(P)
        Ec = m.band(P0, 5 * NM)[1]
        P["Phi_CNL"] = Ec + smp["cnl_off"]
        steps = [{"V": 0.05}, {"VG": -3.0}]
        for kind in GEOMS:  # shared pinning
            T.append(_spec("uq_gated", kind, PT, dict(P, S_top=smp["S"], S_edge=smp["S"]),
                           dict(smp, variant="shared"), steps))
        T.append(_spec("uq_gated", "edge", PT, dict(P, S_top=smp["S"], S_edge=smp["S_edge_indep"]),
                       dict(smp, variant="indep_edge"), steps))
    return T


def tasks_thermal_lit():
    T = []
    variants = [("full", {}), ("contact_sink", dict(h_contact=1e4)), ("R_TBR_3e-8", dict(R_TBR=3e-8 * 1e4))]
    for case, c in s.pin_cases().items():
        for kind in GEOMS:
            for name, extra in variants:
                P = dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"], **extra)
                steps = [{"VG": -3.0, "save": False}, {"V": 0.5}, {"V": 1.0}]
                T.append(_spec("thermal_lit", kind, PT, P, dict(case=case, metal="Pt", variant=name), steps))
    return T


def tasks_mesh_convergence():
    T = []
    for mf in (1.0, 2.0, 3.0):
        for kind in GEOMS:
            # (a) the Table 3 operating point (reference alignment, thermal on)
            T.append(_spec("mesh_convergence", kind, 5.0, {}, dict(case="table3_reference", mf=mf),
                           [{"VG": -2.0, "save": False}, {"V": 0.5}, {"V": 1.0}], mf=mf))
            # (b) tunnelling-dominated measured-pinning on state (direct Pt, interlayer Pt)
            for case in ("direct (measured)", "interlayer (measured)"):
                c = s.pin_cases()[case]
                P = dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"], thermal=False)
                T.append(_spec("mesh_convergence", kind, PT, P, dict(case=case, mf=mf),
                               [{"V": 0.05}, {"VG": -3.0}], mf=mf))
    return T


def tasks_mesh_ungated():
    """Mesh-factor 0.5-3 for the ungated device at the reference alignment (Rtot at 0.05 V, 4.8 and 5.0 eV)."""
    T = []
    for mf in (0.5, 1.0, 2.0, 3.0):
        for kind in ("top", "edge"):
            for Phi in (4.8, 5.0):
                T.append(_spec("mesh_ungated", kind, Phi, dict(thermal=False), dict(case="ungated_reference", mf=mf),
                               [{"V": 0.05}], gate=False, mf=mf))
    return T


def tasks_top_retry():
    """Gated on-state (V_G = -3 V, V_D = 0.05 V, literature alignment) for the cases of metals_pinning that did not
    converge, with two other continuation paths: (A) gate ramped at V_D = 0, then the bias applied; (B) bias first,
    then a very small gate step (0.02 V)."""
    pc = s.pin_cases()
    cases = [("direct (measured)", m_, "top") for m_ in ("Ti", "Cr", "Au", "Pd")]
    cases += [("interlayer (measured)", m_, "top") for m_ in ("Cr", "Ti")]
    cases += [("unpinned", "Cr", "top"), ("unpinned", "Ti", "edge")]
    T = []
    for case, metal, kind in cases:
        c = pc[case]
        P = dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"], thermal=False)
        tag = dict(case=case, metal=metal, S=c["S"], CNL=c["CNL"], study="gated")
        a = _spec("top_retry", kind, s.METALS[metal], P, dict(tag, variant="A_gate_first"),
                  [{"VG": -3.0, "save": False}, {"V": 0.05}])
        a["vg_step"] = 0.1
        b = _spec("top_retry", kind, s.METALS[metal], P, dict(tag, variant="B_small_gate_step"),
                  [{"V": 0.05, "save": False}, {"VG": -3.0}])
        b["vg_step"] = 0.02
        T += [a, b]
    return T


def tasks_thickness_lit():
    T = []
    for N in (3, 4, 5, 6, 8, 12):
        t = N * 0.83 * NM
        for case in ("direct (measured)", "interlayer (measured)"):
            c = _cnl_cases(t)[case]
            for kind in GEOMS:
                P = dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"], thermal=False)
                T.append(_spec("thickness_lit", kind, PT, P, dict(layers=N, case=case, study="literature"),
                               [{"V": 0.05}, {"VG": -3.0}], t=t))
        for frac in (0.2, 0.32, 0.5):  # valence-band share of the confinement shift, reference alignment
            for kind in GEOMS:
                P = dict(dEv_frac=frac, thermal=False)
                T.append(_spec("thickness_lit", kind, 5.0, P, dict(layers=N, dEv_frac=frac, study="reference_S1"),
                               [{"V": 0.05}, {"VG": -3.0}], t=t))
    return T


def tasks_inse_electron():
    T = []
    metals = dict(In=4.12, Ti=4.33, Pd=5.12)
    cases = {"direct (measured)": (0.05, 0.30), "interlayer (measured)": (0.50, None)}
    for T0 in (200.0, 250.0, 300.0, 350.0, 400.0):
        for metal, Phi in metals.items():
            for kind in GEOMS:
                P0 = dict(m.DEFAULTS); P0.update(s.LIT)
                Ec = m.band(P0, 5 * NM)[1]
                P = dict(s.LIT, NA=0.0, ND=1e16, T0=T0, thermal=False, S_top=0.05, S_edge=0.05, Phi_CNL=Ec + 0.30)
                steps = [{"V": 0.05}] + [{"VG": v} for v in (0.0, 1.0, 2.0, 3.0, 4.0)]
                T.append(_spec("inse_electron", kind, Phi, P, dict(metal=metal, S=0.05, T0=T0, carrier="electron"),
                               steps))
    return T


def tasks_chen2019():
    """Reproduce the measurement protocol of Chen et al., npj 2D Mater. Appl. 3, 49 (2019) [Ref. 11]:
    n-type InSe, 10-20 nm thick, top contacts (In, Ti, Cr, Pd), 300 nm SiO2 back gate, Arrhenius barriers
    from 300-400 K, with the reported pinning of control (S=0.05, Phi_CNL=4.9 eV) and oxidized-monolayer
    (S=0.5, Phi_CNL=5.2 eV) contacts. Part B: In-contact on-current at V_G=80 V, V_D=1 V (reported 103 uA/um
    for the oxidized-monolayer device, 3-6x the control)."""
    T = []
    t = 18 * 0.83 * NM                      # 15 nm
    g300 = dict(t_ox=300e-7)
    metals = dict(In=4.12, Ti=4.33, Cr=4.50, Pd=5.12)
    cases = {"control": (0.05, 4.9), "oxidized monolayer": (0.5, 5.2)}
    for case, (S, cnl) in cases.items():
        for metal, Phi in metals.items():
            for T0 in (300.0, 350.0, 400.0):
                P = dict(s.LIT, NA=0.0, ND=1e16, T0=T0, thermal=False, S_top=S, S_edge=S, Phi_CNL=cnl)
                steps = [{"V": 0.05, "save": False}] + [{"VG": v} for v in (-20.0, 0.0, 20.0)]
                sp = _spec("chen2019", "top", Phi, P, dict(part="A_arrhenius", case=case, metal=metal, S=S, CNL=cnl, T0=T0),
                           steps, t=t)
                sp.update(gate_dict=g300, Lc=200 * NM, Lch=500 * NM, vg_step=1.0)
                T.append(sp)
        for Lch in (1000,):
            P = dict(s.LIT, NA=0.0, ND=1e16, thermal=False, S_top=S, S_edge=S, Phi_CNL=cnl)
            steps = [{"V": 1.0, "save": False}] + [{"VG": v} for v in (40.0, 80.0)]
            sp = _spec("chen2019", "top", metals["In"], P, dict(part="B_oncurrent", case=case, metal="In", S=S, CNL=cnl,
                                                               Lch_nm=Lch), steps, t=t)
            sp.update(gate_dict=g300, Lc=500 * NM, Lch=Lch * NM, vg_step=1.0)
            T.append(sp)
    return T


def tasks_chen2019_300K():
    """Re-run of the 300 K points of chen2019 part A without the unreachable negative gate (deep depletion at
    V_G < -5 V in the n-type film): V_G = 0 and +20 V, shared with the 350 K and 400 K runs."""
    T = []
    t = 18 * 0.83 * NM
    metals = dict(In=4.12, Ti=4.33, Cr=4.50, Pd=5.12)
    cases = {"control": (0.05, 4.9), "oxidized monolayer": (0.5, 5.2)}
    for case, (S, cnl) in cases.items():
        for metal, Phi in metals.items():
            P = dict(s.LIT, NA=0.0, ND=1e16, T0=300.0, thermal=False, S_top=S, S_edge=S, Phi_CNL=cnl)
            steps = [{"V": 0.05, "save": False}] + [{"VG": v} for v in (0.0, 20.0)]
            sp = _spec("chen2019_300K", "top", Phi, P, dict(part="A_arrhenius", case=case, metal=metal, S=S, CNL=cnl, T0=300.0),
                       steps, t=t)
            sp.update(gate_dict=dict(t_ox=300e-7), Lc=200 * NM, Lch=500 * NM, vg_step=1.0)
            T.append(sp)
    return T


def tasks_chen2019_mu():
    """Sensitivity of the chen2019 part-B on-currents to the assumed electron mobility (Chen et al. used InSe on
    SiO2, where mobilities are typically far below the 1000 cm2/Vs baseline) and to the unreported channel length."""
    T = []
    t = 18 * 0.83 * NM
    cases = {"control": (0.05, 4.9), "oxidized monolayer": (0.5, 5.2)}
    grid = [(mu, 1000) for mu in (100.0, 200.0, 500.0)] + [(1000.0, 3000)]
    for case, (S, cnl) in cases.items():
        for mu, Lch in grid:
            P = dict(s.LIT, NA=0.0, ND=1e16, thermal=False, S_top=S, S_edge=S, Phi_CNL=cnl, mu_n=mu)
            steps = [{"V": 1.0, "save": False}] + [{"VG": v} for v in (40.0, 80.0)]
            sp = _spec("chen2019_mu", "top", 4.12, P, dict(part="B_oncurrent", case=case, metal="In", S=S, CNL=cnl,
                                                           Lch_nm=Lch, mu_n=mu), steps, t=t)
            sp.update(gate_dict=dict(t_ox=300e-7), Lc=500 * NM, Lch=Lch * NM, vg_step=1.0)
            T.append(sp)
    return T


def tasks_mh_sweep():
    """Sensitivity of the gated on state (V_G = -3 V, V_D = 0.05 V, literature alignment, Pt) to the hole mass.
    The baseline 0.17 m0 is light compared with the flat, ring-shaped ('Mexican-hat') valence band of few-layer
    InSe. Variant 'tun': only the tunnelling mass changes (DOS mass 0.17); variant 'both': DOS and tunnelling."""
    pc = s.pin_cases()
    T = []
    for case in ("direct (measured)", "interlayer (measured)", "unpinned"):
        c = pc[case]
        for kind in GEOMS:
            for mh in (0.17, 0.5, 1.0, 2.0):
                for var in (("tun", "both") if mh != 0.17 else ("tun",)):
                    P = dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"], thermal=False, mh_tun=mh,
                             mh_dos=(mh if var == "both" else 0.17))
                    sp = _spec("mh_sweep", kind, PT, P, dict(case=case, metal="Pt", S=c["S"], mh=mh, variant=var,
                                                              study="gated"),
                               [{"VG": -3.0, "save": False}, {"V": 0.05}])
                    sp["vg_step"] = 0.1
                    T.append(sp)
            if kind == "top":  # direction-dependent: light c-axis tunnelling mass, heavier in-plane DOS mass
                for mdos in (0.5, 1.0, 2.0):
                    P = dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"], thermal=False, mh_tun=0.17,
                             mh_dos=mdos)
                    sp = _spec("mh_sweep", kind, PT, P, dict(case=case, metal="Pt", S=c["S"], mh=mdos, variant="dir",
                                                              study="gated"),
                               [{"VG": -3.0, "save": False}, {"V": 0.05}])
                    sp["vg_step"] = 0.1
                    T.append(sp)
    return T


def tasks_oblique():
    """Gated on state (V_G=-3 V, V_D=0.05 V, literature alignment) with tunnelling along the contact normal and
    along the minimum-action straight path within 75 degrees of the normal (both geometries)."""
    Ev = s.lit_Ev()
    pc = s.pin_cases()
    cases = [(f"phiB={b}", Ev - b, dict(s.LIT)) for b in (0.2, 0.6, 1.0)]
    for case in ("direct (measured)", "interlayer (measured)"):
        c = pc[case]
        cases.append((case, PT, dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"])))
    T = []
    for case, Phi, P0 in cases:
        for kind in GEOMS:
            for obl in (False, True):
                P = dict(P0, thermal=False, oblique=obl)
                sp = _spec("oblique", kind, Phi, P, dict(case=case, oblique=int(obl), study="gated"),
                           [{"VG": -3.0, "save": False}, {"V": 0.05}])
                sp["vg_step"] = 0.1
                T.append(sp)
    return T


def tasks_uq_cnl():
    """Edge contacts of the gated Monte Carlo with a pinning level drawn independently of the top contact.
    Reproduces the samples of tasks_uq_gated (same seed and draw order); top-contact results are reused."""
    N = int(os.environ.get("UQ_N", "240"))
    rng = np.random.default_rng(7)
    rng2 = np.random.default_rng(11)
    T = []
    for k in range(N):
        lg = lambda a, b: 10 ** rng.uniform(math.log10(a), math.log10(b))
        smp = dict(sample=k, dEV=rng.normal(0, 0.3), S=lg(0.05, 0.5), S_edge_indep=lg(0.05, 0.5),
                   cnl_off=rng.uniform(0.1, 0.6), d_vdW=rng.uniform(0.2, 0.4), Phi_vdW=rng.uniform(1, 3),
                   mu_p=lg(20, 100), rmu=lg(0.002, 0.1), NA=lg(1e15, 1e17), eps_r=rng.uniform(5, 8),
                   mh=rng.uniform(0.1, 0.3))
        smp["cnl_edge"] = rng2.uniform(0.1, 0.6)
        P = dict(s.LIT, chi_bulk=s.LIT["chi_bulk"] - smp["dEV"], d_vdW=smp["d_vdW"] * NM,
                 Phi_vdW=smp["Phi_vdW"], mu_p=smp["mu_p"], mu_perp_ratio=smp["rmu"], NA=smp["NA"],
                 eps_r=smp["eps_r"], mh_dos=smp["mh"], mh_tun=smp["mh"], thermal=False)
        P0 = dict(m.DEFAULTS); P0.update(P)
        Ec = m.band(P0, 5 * NM)[1]
        P["Phi_CNL"] = Ec + smp["cnl_edge"]
        steps = [{"V": 0.05}, {"VG": -3.0}]
        T.append(_spec("uq_cnl", "edge", PT, dict(P, S_top=smp["S"], S_edge=smp["S"]),
                       dict(smp, variant="indep_cnl"), steps))
        T.append(_spec("uq_cnl", "edge", PT, dict(P, S_top=smp["S"], S_edge=smp["S_edge_indep"]),
                       dict(smp, variant="indep_cnl_S"), steps))
    return T


def tasks_thermo_const():
    """Table-3 operating point with the contact Peltier constant of the SG flux (1.5) instead of 2."""
    T = []
    for pcst in (2.0, 1.5):
        for kind in GEOMS:
            T.append(_spec("thermo_const", kind, 5.0, dict(pelt_const=pcst), dict(pelt_const=pcst),
                           [{"VG": -2.0, "save": False}, {"V": 0.5}, {"V": 1.0}]))
    return T


def tasks_heat_step():
    """Thermal response to switching on the heat sources at fixed bias, with and without tunnelling."""
    T = []
    for kind in GEOMS:
        for tun in ("nonlocal", "off"):
            T.append(_spec("heat_step", kind, 5.0, dict(tunneling=tun), dict(case="reference", tunneling=tun),
                           [{"VG": -2.0, "save": False}, {"V": 1.0, "save": False}, {"heat_step": 1e-6}]))
    c = s.pin_cases()["unpinned"]
    for kind in GEOMS:
        T.append(_spec("heat_step", kind, PT, dict(s.LIT, S_top=1.0, S_edge=1.0, Phi_CNL=c["CNL"]),
                       dict(case="unpinned Pt", tunneling="nonlocal"),
                       [{"VG": -3.0, "save": False}, {"V": 1.0, "save": False}, {"heat_step": 1e-6}]))
    return T


def _dir_mass(P, kind, m_in, m_c=0.17):
    """Direction-dependent hole mass: top contacts tunnel along c (light), edge contacts in-plane (heavy);
    the density-of-states mass is the in-plane value for both."""
    return dict(P, mh_dos=m_in, mh_tun=(m_c if kind == "top" else m_in))


def tasks_mh_dir_map():
    """Top/edge on-state ratio with direction-dependent hole masses versus in-plane mass and barrier."""
    Ev = s.lit_Ev()
    pc = s.pin_cases()
    cases = [(f"phiB={b}", Ev - b, dict(s.LIT)) for b in (0.2, 0.6, 1.0)]
    for case in ("direct (measured)", "interlayer (measured)", "unpinned"):
        c = pc[case]
        cases.append((case, PT, dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"])))
    T = []
    for case, Phi, P0 in cases:
        for m_in in (0.3, 0.5, 0.7, 0.9, 1.2):
            for kind in GEOMS:
                P = _dir_mass(dict(P0, thermal=False), kind, m_in)
                sp = _spec("mh_dir_map", kind, Phi, P, dict(case=case, m_in=m_in, study="gated"),
                           [{"VG": -3.0, "save": False}, {"V": 0.05}])
                sp["vg_step"] = 0.1
                T.append(sp)
    return T


def tasks_uq_dir():
    """Gated Monte Carlo (same samples as tasks_uq_gated) with direction-dependent hole masses:
    c-axis tunnelling mass 0.17 m0 for top contacts, in-plane mass log-uniform 0.3-1.5 m0 for edge tunnelling
    and for the density of states of both."""
    N = int(os.environ.get("UQ_N", "240"))
    rng = np.random.default_rng(7)
    rng2 = np.random.default_rng(11)
    rng3 = np.random.default_rng(13)
    T = []
    for k in range(N):
        lg = lambda a, b: 10 ** rng.uniform(math.log10(a), math.log10(b))
        smp = dict(sample=k, dEV=rng.normal(0, 0.3), S=lg(0.05, 0.5), S_edge_indep=lg(0.05, 0.5),
                   cnl_off=rng.uniform(0.1, 0.6), d_vdW=rng.uniform(0.2, 0.4), Phi_vdW=rng.uniform(1, 3),
                   mu_p=lg(20, 100), rmu=lg(0.002, 0.1), NA=lg(1e15, 1e17), eps_r=rng.uniform(5, 8),
                   mh=rng.uniform(0.1, 0.3))
        smp["m_in"] = 10 ** rng3.uniform(math.log10(0.3), math.log10(1.5))
        smp["cnl_edge"] = rng2.uniform(0.1, 0.6)
        P = dict(s.LIT, chi_bulk=s.LIT["chi_bulk"] - smp["dEV"], d_vdW=smp["d_vdW"] * NM,
                 Phi_vdW=smp["Phi_vdW"], mu_p=smp["mu_p"], mu_perp_ratio=smp["rmu"], NA=smp["NA"],
                 eps_r=smp["eps_r"], thermal=False)
        P0 = dict(m.DEFAULTS); P0.update(P)
        Ec = m.band(P0, 5 * NM)[1]
        P["Phi_CNL"] = Ec + smp["cnl_off"]
        steps = [{"V": 0.05, "save": False}, {"VG": -3.0}]
        for kind in GEOMS:
            Pk = _dir_mass(dict(P, S_top=smp["S"], S_edge=smp["S"]), kind, smp["m_in"])
            T.append(_spec("uq_dir", kind, PT, Pk, dict(smp, variant="shared"), steps))
        Pk = _dir_mass(dict(P, S_top=smp["S"], S_edge=smp["S_edge_indep"]), "edge", smp["m_in"])
        T.append(_spec("uq_dir", "edge", PT, Pk, dict(smp, variant="indep_edge"), steps))
        Pk = _dir_mass(dict(P, S_top=smp["S"], S_edge=smp["S"], Phi_CNL=Ec + smp["cnl_edge"]), "edge", smp["m_in"])
        T.append(_spec("uq_dir", "edge", PT, Pk, dict(smp, variant="indep_cnl"), steps))
    return T


def tasks_thickness_dir():
    """Thickness dependence at the literature alignment (Pt, V_G=-3 V) with direction-dependent hole mass
    (top: 0.17 m0 tunnelling; edge and DOS: 0.9 m0) and, for reference, the light isotropic mass."""
    T = []
    for N in (3, 4, 6, 8, 12):
        t = N * 0.83 * NM
        for case in ("direct (measured)", "interlayer (measured)"):
            c = _cnl_cases(t)[case]
            for kind in GEOMS:
                P0 = dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"], thermal=False)
                sp = _spec("thickness_dir", kind, PT, _dir_mass(P0, kind, 0.9), dict(layers=N, case=case, masses="dir0.9"),
                           [{"VG": -3.0, "save": False}, {"V": 0.05}], t=t)
                sp["vg_step"] = 0.1
                T.append(sp)
    return T


def _tensor(P, m_par, m_c=0.17, disp="kane", obl=False, m_dos=None):
    """Anisotropic hole mass: tunnelling mass along each path from the inverse-mass tensor (c axis m_c,
    in-plane m_par), density-of-states mass (m_par^2 m_c)^(1/3), imaginary dispersion 'parabolic' or 'kane'."""
    md = m_dos if m_dos is not None else (m_par * m_par * m_c) ** (1.0 / 3.0)
    return dict(P, mh_tun_c=m_c, mh_tun_par=m_par, mh_dos=md, tun_disp=disp, oblique=obl)


def _mass_cases():
    Ev = s.lit_Ev()
    pc = s.pin_cases()
    cases = [(f"phiB={b}", Ev - b, dict(s.LIT)) for b in (0.2, 0.6, 1.0)]
    for case in ("direct (measured)", "interlayer (measured)", "unpinned"):
        c = pc[case]
        cases.append((case, PT, dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"])))
    return cases


def _gs(out, kind, Phi, P, tag, vg_step=0.1, t=5 * NM):
    sp = _spec(out, kind, Phi, P, tag, [{"VG": -3.0, "save": False}, {"V": 0.05}], t=t)
    sp["vg_step"] = vg_step
    return sp


def tasks_mass_controls():
    """Separate the hole-mass anisotropy from the direction of the gate field (m=0.9, Kane, normal paths):
    'iso': 0.9 m0 in every direction; 'swapped': heavy along c (0.9), light in-plane (0.17); DOS mass fixed."""
    md = (0.9 * 0.9 * 0.17) ** (1.0 / 3.0)
    T = []
    for case, Phi, P0 in _mass_cases():
        for var, mc, mp in (("iso", 0.9, 0.9), ("swapped", 0.9, 0.17)):
            for kind in GEOMS:
                P = _tensor(dict(P0, thermal=False), mp, m_c=mc, m_dos=md)
                T.append(_gs("mass_controls", kind, Phi, P, dict(case=case, variant=var)))
    return T


def tasks_edge_T():
    """Edge-contact interface transmission needed to reverse the low-barrier ranking (Kane, tensor 0.9)."""
    Ev = s.lit_Ev()
    pc = s.pin_cases()
    cases = [(f"phiB={b}", Ev - b, dict(s.LIT)) for b in (0.2, 0.6)]
    for case in ("interlayer (measured)", "unpinned"):
        c = pc[case]
        cases.append((case, PT, dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"])))
    T = []
    for case, Phi, P0 in cases:
        for Te in (0.3, 0.1, 0.03, 0.01, 0.003, 0.001):
            P = _tensor(dict(P0, thermal=False, T_edge=Te), 0.9)
            T.append(_gs("edge_T", "edge", Phi, P, dict(case=case, T_edge=Te)))
    return T


def tasks_kane_map():
    """Top/edge ratio vs in-plane mass with the mass tensor: parabolic and Kane dispersion, and Kane with
    minimum-action oblique paths."""
    T = []
    for case, Phi, P0 in _mass_cases():
        for m_par in (0.3, 0.5, 0.7, 0.9, 1.2):
            for var, disp, obl in (("parabolic", "parabolic", False), ("kane", "kane", False), ("kane_oblique", "kane", True)):
                for kind in GEOMS:
                    P = _tensor(dict(P0, thermal=False), m_par, disp=disp, obl=obl)
                    T.append(_gs("kane_map", kind, Phi, P, dict(case=case, m_par=m_par, variant=var)))
    return T


def tasks_thickness_kane():
    T = []
    for N in (3, 4, 6, 8, 12):
        t = N * 0.83 * NM
        for case in ("direct (measured)", "interlayer (measured)"):
            c = _cnl_cases(t)[case]
            for kind in GEOMS:
                P = _tensor(dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"], thermal=False), 0.9)
                T.append(_gs("thickness_kane", kind, PT, P, dict(layers=N, case=case), vg_step=0.02, t=t))
    return T


def tasks_uq_kane():
    """Gated Monte Carlo (same samples as uq_gated/uq_dir) with the mass tensor and Kane dispersion."""
    N = int(os.environ.get("UQ_N", "240"))
    rng = np.random.default_rng(7)
    rng2 = np.random.default_rng(11)
    rng3 = np.random.default_rng(13)
    T = []
    for k in range(N):
        lg = lambda a, b: 10 ** rng.uniform(math.log10(a), math.log10(b))
        smp = dict(sample=k, dEV=rng.normal(0, 0.3), S=lg(0.05, 0.5), S_edge_indep=lg(0.05, 0.5),
                   cnl_off=rng.uniform(0.1, 0.6), d_vdW=rng.uniform(0.2, 0.4), Phi_vdW=rng.uniform(1, 3),
                   mu_p=lg(20, 100), rmu=lg(0.002, 0.1), NA=lg(1e15, 1e17), eps_r=rng.uniform(5, 8),
                   mh=rng.uniform(0.1, 0.3))
        smp["m_in"] = 10 ** rng3.uniform(math.log10(0.3), math.log10(1.5))
        smp["cnl_edge"] = rng2.uniform(0.1, 0.6)
        P = dict(s.LIT, chi_bulk=s.LIT["chi_bulk"] - smp["dEV"], d_vdW=smp["d_vdW"] * NM,
                 Phi_vdW=smp["Phi_vdW"], mu_p=smp["mu_p"], mu_perp_ratio=smp["rmu"], NA=smp["NA"],
                 eps_r=smp["eps_r"], thermal=False)
        P0 = dict(m.DEFAULTS); P0.update(P)
        Ec = m.band(P0, 5 * NM)[1]
        P["Phi_CNL"] = Ec + smp["cnl_off"]
        steps = [{"V": 0.05, "save": False}, {"VG": -3.0}]
        for kind in GEOMS:
            T.append(_spec("uq_kane", kind, PT, _tensor(dict(P, S_top=smp["S"], S_edge=smp["S"]), smp["m_in"]),
                           dict(smp, variant="shared"), steps))
        T.append(_spec("uq_kane", "edge", PT, _tensor(dict(P, S_top=smp["S"], S_edge=smp["S_edge_indep"]), smp["m_in"]),
                       dict(smp, variant="indep_edge"), steps))
        T.append(_spec("uq_kane", "edge", PT, _tensor(dict(P, S_top=smp["S"], S_edge=smp["S"], Phi_CNL=Ec + smp["cnl_edge"]),
                                                      smp["m_in"]), dict(smp, variant="indep_cnl"), steps))
    return T


def tasks_mass_phase():
    """Geometry phase diagram in (m_c, m_par) at four barriers: two-band dispersion, most favourable straight path.
    m_c spans the GW c-axis value (0.078), the bulk value (0.17) and heavy effective values up to the out-of-plane
    transport mass reported for thick flakes (2.5)."""
    Ev = s.lit_Ev()
    pc = s.pin_cases()
    cases = [("phiB=0.8", Ev - 0.8, dict(s.LIT))]
    for case in ("direct (measured)", "interlayer (measured)", "unpinned"):
        c = pc[case]
        cases.append((case, PT, dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"])))
    T = []
    for case, Phi, P0 in cases:
        for m_c in (0.078, 0.17, 0.4, 0.9, 2.5):
            for m_par in (0.3, 0.6, 0.9, 1.9):
                for kind in GEOMS:
                    P = _tensor(dict(P0, thermal=False), m_par, m_c=m_c, disp="kane", obl=True)
                    T.append(_gs("mass_phase", kind, Phi, P, dict(case=case, m_c=m_c, m_par=m_par)))
    return T


STUDIES = {"uq_gated": tasks_uq_gated, "thermal_lit": tasks_thermal_lit, "mesh_convergence": tasks_mesh_convergence,
           "thickness_lit": tasks_thickness_lit, "inse_electron": tasks_inse_electron,
           "chen2019": tasks_chen2019, "chen2019_300K": tasks_chen2019_300K,
           "chen2019_mu": tasks_chen2019_mu,
           "mesh_ungated": tasks_mesh_ungated,
           "top_retry": tasks_top_retry, "mh_sweep": tasks_mh_sweep, "oblique": tasks_oblique,
           "uq_cnl": tasks_uq_cnl, "thermo_const": tasks_thermo_const, "heat_step": tasks_heat_step,
           "mh_dir_map": tasks_mh_dir_map, "uq_dir": tasks_uq_dir, "thickness_dir": tasks_thickness_dir,
           "mass_controls": tasks_mass_controls, "edge_T": tasks_edge_T, "kane_map": tasks_kane_map,
           "thickness_kane": tasks_thickness_kane, "uq_kane": tasks_uq_kane, "mass_phase": tasks_mass_phase}
