"""v3 studies. Usage: python studies.py <study> [...]  -> ../data_v3/<study>.csv"""
import csv, math, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import tcad2d as m  # noqa

OUT = os.path.join(os.path.dirname(__file__), "..", "data_v3")
os.makedirs(OUT, exist_ok=True)
NM = 1e-7
GATE = dict(t_ox=10e-7)


def write(name, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(OUT, name + ".csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print("RES wrote", name, len(rows), flush=True)


def lab(kind, Lc):
    return "edge" if kind == "edge" else f"top{Lc / NM:.0f}"


def summarize(dv, extra=None):
    I, Ip, In, It = dv.current()
    dec = dv.decompose()
    pb = dv.power_balance()
    r = dict(geom=lab(dv.kind, dv.geo["Lc"]), Lc_nm=dv.geo["Lc"] / NM, Lch_nm=dv.geo["Lch"] / NM,
             t_nm=dv.t / NM, Eg=dv.Eg, chi=dv.chi, Phi_M=dv.Phi_M, phiBp=dv.phiBp, V=dv.bias, VG=dv.vg,
             I=I, I_tun=It, tun_share=It / I if I else np.nan, R_tot=dec["R_tot"], R_c=dec["R_c"],
             R_ch=dec["R_ch"], frac_contact=dec["frac_contact"], Tmax=pb.get("Tmax", np.nan),
             P_el=pb["P_el"], P_out=pb.get("P_out", np.nan), T_source=pb.get("T_source", np.nan),
             T_drain=pb.get("T_drain", np.nan), P_junction=pb.get("P_junction", np.nan),
             P_tunnel=pb.get("P_tunnel", np.nan), P_bulk=pb.get("P_bulk", np.nan))
    if dv.gate:
        ps, ns = dv.sheet_density()
        r.update(p_sheet=ps, n_sheet=ns)
    if extra:
        r.update(extra)
    print("RES", {k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    return r


def point(kind, Phi, V=0.05, P=None, Lc=50 * NM, Lch=50 * NM, t=5 * NM, gate=None, VG=None, extra=None):
    try:
        dv = m.Device(kind, P=P, Phi_M=Phi, Lc=Lc, Lch=Lch, t=t, gate=gate)
        if VG is not None:
            dv.set_gate(VG)
        dv.set_bias(V)
        return summarize(dv, extra), dv
    except Exception as e:  # record and continue
        print("RES FAILED", kind, Phi, P, extra, str(e)[:200], flush=True)
        r = dict(geom=lab(kind, Lc), Lc_nm=Lc / NM, Phi_M=Phi, V=V, failed=1)
        r.update(extra or {})
        return r, None


GEOMS = [("top", 10 * NM), ("top", 50 * NM), ("edge", 50 * NM)]
WF = [4.33, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.4, 5.65]


def s_wf():
    rows = []
    for kind, Lc in GEOMS:
        for Phi in WF:
            rows.append(point(kind, Phi, Lc=Lc)[0])
        for Phi in (5.0, 5.2, 5.4):  # Boltzmann (clipped) comparison
            rows.append(point(kind, Phi, Lc=Lc, P=dict(statistics="boltzmann"), extra=dict(stat="boltzmann"))[0])
        for Phi in (4.6, 4.8, 5.0):  # tunnelling model comparison (ungated)
            for tun in ("off", "local"):
                rows.append(point(kind, Phi, Lc=Lc, P=dict(tunneling=tun), extra=dict(tunneling=tun))[0])
    write("wf_sweep", rows)


def s_aniso():
    rows = []
    for rmu in (1.0, 0.1, 0.01, 0.002):
        for kind, Lc in GEOMS:
            for Phi in (4.8, 5.0):
                rows.append(point(kind, Phi, Lc=Lc, P=dict(mu_perp_ratio=rmu), extra=dict(rmu=rmu))[0])
    for rmu in (1.0, 0.01):
        for Lc in (5, 10, 20, 50, 100):
            for Phi in (4.8, 5.0):
                rows.append(point("top", Phi, Lc=Lc * NM, P=dict(mu_perp_ratio=rmu), extra=dict(rmu=rmu, study="Lc"))[0])
    for kp in (0.085, 0.0076):  # thermal anisotropy at high bias (ungated)
        for kind, Lc in GEOMS:
            rows.append(point(kind, 5.0, V=1.0, Lc=Lc, P=dict(kappa_perp=kp), extra=dict(kperp=kp, study="thermal"))[0])
    write("anisotropy", rows)


def s_pinning():
    rows = []
    for Se in (1.0, 0.7, 0.5, 0.3):
        for Phi in (4.6, 4.8, 5.0, 5.2, 5.4, 5.65):
            rows.append(point("edge", Phi, P=dict(S_edge=Se), extra=dict(S_edge=Se))[0])
    for St in (1.0, 0.5):
        for Phi in (4.6, 4.8, 5.0, 5.2, 5.4, 5.65):
            rows.append(point("top", Phi, P=dict(S_top=St), extra=dict(S_top=St))[0])
    write("pinning_top_edge", rows)


def transfer(kind, Phi, VD, P=None, Lc=50 * NM, VGs=None, tag=None):
    rows = []
    try:
        dv = m.Device(kind, P=P, Phi_M=Phi, Lc=Lc, gate=GATE)
        dv.set_bias(VD)
        for VG in VGs:
            dv.set_gate(VG)
            rows.append(summarize(dv, dict(tag or {})))
    except Exception as e:
        print("RES FAILED transfer", kind, Phi, P, str(e)[:200], flush=True)
    return rows


def s_gated():
    rows = []
    VGs = list(np.round(np.arange(1.0, -4.01, -0.5), 3))
    for Phi in (4.8, 5.0, 5.2):
        for kind, Lc in [("top", 50 * NM), ("edge", 50 * NM)]:
            rows += transfer(kind, Phi, 0.05, Lc=Lc, VGs=VGs, tag=dict(study="transfer"))
    for kind, Lc in [("top", 50 * NM), ("edge", 50 * NM)]:
        for tun in ("off", "local"):
            rows += transfer(kind, 4.8, 0.05, P=dict(tunneling=tun), Lc=Lc, VGs=VGs,
                             tag=dict(study="transfer", tunneling=tun))
    for kind, Lc in [("top", 50 * NM), ("edge", 50 * NM)]:  # output characteristics, on state
        for Phi in (4.8, 5.0):
            try:
                dv = m.Device(kind, Phi_M=Phi, Lc=Lc, gate=GATE)
                dv.set_gate(-3.0)
                for V in (0.05, 0.25, 0.5, 1.0, 1.5, 2.0):
                    dv.set_bias(V)
                    rows.append(summarize(dv, dict(study="output")))
            except Exception as e:
                print("RES FAILED output", kind, Phi, str(e)[:200], flush=True)
        write("gated", rows)
    write("gated", rows)


def s_thickness():
    rows = []
    for N in (3, 4, 6, 8, 12):
        t = N * 0.83 * NM
        for kind, Lc in [("top", 50 * NM), ("edge", 50 * NM)]:
            for Phi in (4.8, 5.0, 5.2):
                rows.append(point(kind, Phi, Lc=Lc, t=t, extra=dict(layers=N, study="ungated"))[0])
            rows.append(point(kind, 5.0, Lc=Lc, t=t, gate=GATE, VG=-3.0, extra=dict(layers=N, study="gated_on"))[0])
    write("thickness", rows)


def s_compact():
    rows = []
    for Phi in (4.75, 4.8, 4.85, 4.9, 4.95):
        for d in (0.2, 0.3, 0.4):
            for Lc in (10, 20, 50, 100):
                for rmu in (1.0, 0.01):
                    P = dict(d_vdW=d * NM, mu_perp_ratio=rmu, tunneling="off", thermal=False)
                    r, dv = point("top", Phi, Lc=Lc * NM, P=P, extra=dict(d_vdW_nm=d, rmu=rmu, T_vdW=m.wkb_vdw(dict(m.DEFAULTS, **P))))
                    # under-contact hole density (centre of contact, mid-film) for the compact model
                    cl = dv.centerline()
                    r["p_under"] = float(np.interp(0.5 * Lc * NM, cl["x"], cl["p"]))
                    r["p_channel"] = float(np.interp(dv.geo["x_ch0"] + 0.5 * dv.geo["Lch"], cl["x"], cl["p"]))
                    rows.append(r)
    write("compact_model_data", rows)


def s_uq():
    rng = np.random.default_rng(2026)
    rows = []
    N = int(os.environ.get("UQ_N", "120"))
    for k in range(N):
        smp = dict(dEV=rng.normal(0, 0.1), S=rng.uniform(0.3, 1.0), d_vdW=rng.uniform(0.2, 0.4),
                   Phi_vdW=rng.uniform(1.0, 3.0), mu_p=10 ** rng.uniform(math.log10(20), 2),
                   rmu=10 ** rng.uniform(math.log10(0.002), -1), NA=10 ** rng.uniform(15, 17),
                   eps_r=rng.uniform(5, 8), mh=rng.uniform(0.1, 0.3))
        for kind, Lc in [("top", 50 * NM), ("edge", 50 * NM)]:
            for Phi in (4.8, 5.0):
                P = dict(chi_bulk=m.INSE["chi_bulk"] - smp["dEV"], S_top=smp["S"], S_edge=smp["S"],
                         d_vdW=smp["d_vdW"] * NM, Phi_vdW=smp["Phi_vdW"], mu_p=smp["mu_p"],
                         mu_perp_ratio=smp["rmu"], NA=smp["NA"], eps_r=smp["eps_r"], mh_dos=smp["mh"],
                         mh_tun=smp["mh"], thermal=False)
                try:
                    r, _ = point(kind, Phi, Lc=Lc, P=P, extra=dict(sample=k, **smp))
                except Exception as e:  # noqa
                    r = dict(geom=lab(kind, Lc), Phi_M=Phi, sample=k, failed=1, **smp)
                rows.append(r)
        if k % 10 == 9:
            write("uq_samples", rows)
    write("uq_samples", rows)


def s_mos2():
    """Calibrate the vdW transmission of Au/MoS2 top contacts to R_c = 740 Ohm um at n = 1e13 cm^-2,
    300 K (English et al. 2016), then predict L_T(300 K) and rho_c, L_T at 400 K."""
    rows = []
    base = dict(m.MOS2, NA=0.0, ND=1e16, thermal=False, S_top=1.0, statistics="fd", tunneling="nonlocal")
    gate = dict(t_ox=90e-7)
    t = 4.5 * NM
    phiBn = 0.15
    Phi = base["chi"] + phiBn  # S = 1 -> phi_Bn = Phi - chi

    def run(TvdW, Lc, T0=300.0, n_target=1e13, Lch=200 * NM, VG0=None):
        P = dict(base, T_vdW_top=TvdW, T0=T0)
        dv = m.Device("top", P=P, Phi_M=Phi, Lc=Lc, Lch=Lch, t=t, gate=gate)
        dv.set_bias(0.02)
        VG = VG0 if VG0 is not None else 42.0
        for _ in range(5):  # secant on V_G for the target channel density
            dv._ramp("gate_bias", "vg", VG, 1.0, max_step=4.0)
            ns = dv.sheet_density()[1]
            if abs(ns / n_target - 1) < 0.01:
                break
            VG = dv.vg + (n_target - ns) * m.q / (3.9 * m.eps0 / 90e-7)
        dec = dv.decompose()
        Rsh = dec["R_ch"] / dv.geo["Lch"]
        return dict(R_c=dec["R_c"], R_c_ohm_um=dec["R_c"] * 1e4, Rsh=Rsh, n=dv.sheet_density()[1], VG=dv.vg,
                    I=dv.current()[0], tun_share=dv.current()[3] / dv.current()[0])

    target = 740e-4  # Ohm cm
    logT, hist, mode = -2.0, [], "T_vdW"
    for it in range(6):  # log-secant, R_c ~ T_vdW^-1/2 when L_c >> L_T
        r = run(10 ** logT, 200 * NM, Lch=100 * NM)
        hist.append((logT, r["R_c_ohm_um"]))
        rows.append(dict(step="calibration_iter", T_vdW=10 ** logT, **r))
        print("RES calib", logT, r, flush=True)
        if abs(r["R_c"] / target - 1) < 0.03:
            break
        if len(hist) >= 2 and hist[-1][0] != hist[-2][0]:
            sl = (math.log10(hist[-1][1]) - math.log10(hist[-2][1])) / (hist[-1][0] - hist[-2][0])
            sl = min(sl, -0.05)
        else:
            sl = -0.5
        logT = min(0.0, logT + (math.log10(target * 1e4) - math.log10(r["R_c_ohm_um"])) / sl)
        if logT >= 0.0 and len(hist) > 1 and hist[-1][0] >= 0.0:
            mode = "T_vdW saturated at 1"
            break
    Tcal = 10 ** logT
    rows.append(dict(step="calibration", T_vdW=Tcal, mode=mode, **r))
    for T0 in (300.0, 400.0):
        for Lc in (20, 35, 50, 100, 200):
            rr = run(Tcal, Lc * NM, T0=T0, Lch=100 * NM)
            rows.append(dict(step="Lc_sweep", T0=T0, Lc_nm=Lc, T_vdW=Tcal, **rr))
            print("RES Lc", T0, Lc, rr, flush=True)
    write("mos2_calibration", rows)


def s_thermal():
    rows = []
    for kind, Lc in [("top", 50 * NM), ("edge", 50 * NM)]:
        for P, tag in [({}, "full"), (dict(peltier=False), "no_peltier"), (dict(interface_heat=False, peltier=False), "bulk_joule_only"),
                       (dict(kappa_perp=0.085), "isotropic_kappa"), (dict(h_contact=1e4), "contact_sink"),
                       (dict(R_TBR=1e-9 * 1e4), "R_TBR_1e-9"), (dict(nonisothermal=False), "isothermal_stats")]:
            try:
                dv = m.Device(kind, P=P, Phi_M=5.0, Lc=Lc, gate=GATE)
                dv.set_gate(-2.0)
                for V in (0.5, 1.0):
                    dv.set_bias(V)
                    rows.append(summarize(dv, dict(case=tag)))
            except Exception as e:
                print("RES FAILED thermal", kind, tag, str(e)[:200], flush=True)
    write("thermal_gated", rows)
    tr = []
    for kind, Lc in [("top", 50 * NM), ("edge", 50 * NM)]:
        dv = m.Device(kind, Phi_M=5.0, Lc=Lc, gate=GATE, P=dict(tunneling="off"))
        dv.set_gate(-2.0); dv.set_bias(0.0)
        res = dv.transient_step(1.0, 1e-6, nsteps=50)
        for t_, Tm, Pel, Pout in res:
            tr.append(dict(geom=lab(kind, Lc), t=t_, Tmax=Tm, P_el=Pel, P_out=Pout))
    write("thermal_transient", tr)


if __name__ == "__main__":
    for s in sys.argv[1:]:
        t0 = time.time()
        globals()["s_" + s]()
        print("RES study", s, "done in", round(time.time() - t0), "s", flush=True)
