"""
Parameter studies for the revised manuscript. Usage:  python studies.py <study> [...]
Each study writes ../data_v2/<study>.csv (and profiles as .npz where relevant).

Studies
  wf        work-function sweep, top (Lc=10,20,50 nm) and edge, V=0.05 V, with decomposition
  wf_iso    same as wf for Lc=50/edge with the thermal equation switched off (consistency)
  pinning   work-function sweep for pinning factors S = 1, 0.5, 0.3
  iv        electro-thermal I-V to 2 V (Phi_M = 4.8, 5.0; top Lc=10/50, edge)
  lch       channel-length sweep (TLM) + I-V for SCLC analysis
  profiles  centre-line and 2D profiles for SCLC / thermal figures
  tunnel    tunnelling on/off, constant-Gamma test (v1 comparison), vdW-gap sweep
  doping    N_A sweep
  material  mu_p, E_V (chi), eps_r sensitivity
  thermalbc kappa, h_sub, contact heat-sinking at 2 V
  mesh      mesh-convergence sweep + Richardson extrapolation
"""
import csv
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import inse2d as m  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "data_v2")
os.makedirs(OUT, exist_ok=True)
NM = 1e-7
V_LOW = 0.05


def write(name, rows):
    if not rows:
        return
    keys = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(OUT, name + ".csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print("RES wrote", name, len(rows))


def label(kind, Lc):
    return "edge" if kind == "edge" else f"top{Lc / NM:.0f}"


def point(kind, Phi, V=V_LOW, P=None, Lc=50 * NM, Lch=50 * NM, mf=1.0, extra=None):
    t0 = time.time()
    dv = m.Device(kind, P=P, Phi_M=Phi, Lc=Lc, Lch=Lch, mf=mf)
    dv.set_bias(V)
    I = dv.current()[0]
    dec = dv.decompose()
    pb = dv.power_balance()
    r = dict(geom=label(kind, Lc), Lc_nm=Lc / NM if kind == "top" else 0, Lch_nm=Lch / NM,
             Phi_M=Phi, Phi_eff=dv.bh["Phi_used"], phiBp=dv.bh["phiBp"],
             clipped=int(dv.bh["clipped"]), V=V, I_A_per_cm=I, R_tot=dec["R_tot"],
             R_ch=dec["R_ch"], R_c=dec["R_c"], frac_contact=dec["frac_contact"],
             R_c_source=dec["R_c_source"], R_c_drain=dec["R_c_drain"], R_if=dec["R_if"],
             R_ch_slope=dec["R_ch_slope"],
             Tmax=pb.get("Tmax", 300.0), P_el=pb["P_el"], P_bulk=pb.get("P_bulk", np.nan),
             P_interface=pb.get("P_interface", np.nan), nodes=dv.nn,
             Gp_max=float(dv.nodev("Gp").max()), T_vdW=dv.T_vdW, secs=time.time() - t0)
    if extra:
        r.update(extra)
    print("RES", {k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()})
    return r, dv


WFS = [4.33, 4.5, 4.6, 4.7, 4.8, 4.9, 4.925, 4.95, 4.975, 5.0, 5.1, 5.2]
GEOMS = [("top", 10 * NM), ("top", 20 * NM), ("top", 50 * NM), ("edge", 0.0)]


def s_wf():
    rows = []
    for kind, Lc in GEOMS:
        for Phi in WFS:
            rows.append(point(kind, Phi, Lc=Lc or 50 * NM)[0])
    write("wf_sweep", rows)


def s_wf_iso():
    rows = []
    for kind, Lc in [("top", 50 * NM), ("edge", 0.0)]:
        for Phi in [4.33, 4.8, 5.0]:
            for th in (True, False):
                rows.append(point(kind, Phi, P=dict(thermal=th), Lc=Lc or 50 * NM,
                                  extra=dict(thermal=int(th)))[0])
    write("wf_thermal_vs_isothermal", rows)


def s_pinning():
    rows = []
    for S in (1.0, 0.5, 0.3):
        for kind, Lc in [("top", 50 * NM), ("edge", 0.0)]:
            for Phi in [4.33, 4.6, 4.8, 5.0, 5.2, 5.4]:
                rows.append(point(kind, Phi, P=dict(S_pin=S), Lc=Lc or 50 * NM,
                                  extra=dict(S_pin=S))[0])
    write("pinning_sweep", rows)


def iv_curve(kind, Phi, Lc=50 * NM, Lch=50 * NM, P=None, Vmax=2.0, dV=0.1, tag=None):
    dv = m.Device(kind, P=P, Phi_M=Phi, Lc=Lc, Lch=Lch)
    rows = []
    Vs = [0.02, 0.05] + list(np.round(np.arange(dV, Vmax + 1e-9, dV), 4))
    for V in Vs:
        dv.set_bias(float(V))
        I = dv.current()[0]
        pb = dv.power_balance()
        dec = dv.decompose()
        r = dict(geom=label(kind, Lc), Lch_nm=Lch / NM, Phi_M=Phi, V=float(V), I_A_per_cm=I,
                 R_tot=V / I, R_ch=dec["R_ch"], R_c=dec["R_c"], R_if=dec["R_if"], Tmax=pb.get("Tmax", 300.0),
                 P_el=pb["P_el"], P_bulk=pb.get("P_bulk", np.nan),
                 P_interface=pb.get("P_interface", np.nan), P_out=pb.get("P_out", np.nan))
        if tag:
            r.update(tag)
        rows.append(r)
        print("RES", kind, Lc, Phi, V, I, pb.get("Tmax"))
    # differential resistance
    Vv = np.array([r["V"] for r in rows]); Iv = np.array([r["I_A_per_cm"] for r in rows])
    Rd = np.gradient(Vv, Iv)
    for r, x in zip(rows, Rd):
        r["R_diff"] = float(x)
    return rows, dv


def s_iv():
    rows = []
    for Phi in (4.8, 5.0):
        for kind, Lc in [("top", 10 * NM), ("top", 50 * NM), ("edge", 50 * NM)]:
            rows += iv_curve(kind, Phi, Lc=Lc)[0]
    write("iv_electrothermal", rows)


def s_lch():
    rows, iv = [], []
    for Phi in (4.8, 5.0):
        for kind in ("top", "edge"):
            for Lch in (25, 50, 100, 200):
                rows.append(point(kind, Phi, Lch=Lch * NM)[0])
                if Phi == 5.0:
                    iv += iv_curve(kind, Phi, Lch=Lch * NM, Vmax=2.0, dV=0.2,
                                   P=dict(thermal=False))[0]
    write("lch_sweep", rows)
    write("lch_iv_isothermal", iv)


def s_profiles():
    prof = {}
    for kind in ("top", "edge"):
        for Phi in (4.8, 5.0):
            for Lch in (50, 100):
                dv = m.Device(kind, Phi_M=Phi, Lch=Lch * NM)
                for V in (0.05, 0.5, 1.0, 2.0):
                    dv.set_bias(V)
                    cl = dv.centerline()
                    key = f"{kind}_Phi{Phi}_L{Lch}_V{V}"
                    E = -np.gradient(cl["Potential"], cl["x"])
                    phip = cl["phip"]
                    I = dv.current()[0]
                    prof[key + "_x"] = cl["x"]
                    prof[key + "_p"] = cl["Holes"]
                    prof[key + "_psi"] = cl["Potential"]
                    prof[key + "_E"] = E
                    prof[key + "_phip"] = phip
                    prof[key + "_T"] = cl.get("T", np.full_like(cl["x"], 300.0))
                    prof[key + "_I"] = np.array([I])
                    if V == 2.0 and Lch == 50:
                        prof[key + "_2d_x"] = dv.x
                        prof[key + "_2d_y"] = dv.y
                        prof[key + "_2d_p"] = dv.nodev("Holes")
                        prof[key + "_2d_T"] = dv.nodev("T")
                    print("RES profile", key, I)
    np.savez_compressed(os.path.join(OUT, "profiles.npz"), **prof)
    print("RES wrote profiles")


def s_tunnel():
    rows = []
    for kind, Lc in [("top", 50 * NM), ("edge", 0.0)]:
        for Phi in (4.33, 4.6, 4.8, 4.9, 5.0):
            for tun in (True, False):
                rows.append(point(kind, Phi, P=dict(tunneling=tun), Lc=Lc or 50 * NM,
                                  extra=dict(study="tunneling", tunneling=int(tun)))[0])
    # v1-style constant multiplier (1+Gamma) with Gamma = 1, 10, 100, applied as a
    # constant interface-velocity enhancement in the corrected Robin BC
    for kind, Lc in [("top", 50 * NM), ("edge", 0.0)]:
        for Phi in (4.6, 4.8, 5.0):
            for G in (1.0, 10.0, 100.0):
                dvP = dict(tunneling=False)
                r, dv = point(kind, Phi, P=dvP, Lc=Lc or 50 * NM)
                # rebuild with constant Gamma
                dv = m.Device(kind, P=dvP, Phi_M=Phi, Lc=Lc or 50 * NM)
                import devsim as ds
                ds.set_node_values(device=dv.name, region=m.REGION, name="Gp", values=[G] * dv.nn)
                ds.set_node_values(device=dv.name, region=m.REGION, name="Gn", values=[G] * dv.nn)
                dv._solve()
                dv.set_bias(V_LOW)
                I = dv.current()[0]
                rows.append(dict(geom=label(kind, Lc or 50 * NM), Phi_M=Phi, V=V_LOW,
                                 I_A_per_cm=I, R_tot=V_LOW / I, study="const_gamma",
                                 Gamma_const=G))
                print("RES constG", kind, Phi, G, V_LOW / I)
    for d in (0.0, 0.2, 0.3, 0.4, 0.5):
        for Phi in (4.8, 5.0):
            rows.append(point("top", Phi, P=dict(d_vdW=d * NM), Lc=50 * NM,
                              extra=dict(study="vdw", d_vdW_nm=d))[0])
    write("tunneling_vdw", rows)


def s_doping():
    rows = []
    for NA in (1e15, 1e16, 1e17):
        for kind, Lc in [("top", 50 * NM), ("edge", 0.0)]:
            for Phi in (4.6, 4.8, 5.0):
                rows.append(point(kind, Phi, P=dict(NA=NA), Lc=Lc or 50 * NM,
                                  extra=dict(NA=NA))[0])
    write("doping_sweep", rows)


def s_material():
    rows = []
    cases = [("mu_p", v, dict(mu_p=v)) for v in (20.0, 50.0, 100.0)]
    cases += [("E_V", 4.9, dict(chi=3.6)), ("E_V", 5.0, dict(chi=3.7)), ("E_V", 5.1, dict(chi=3.8))]
    cases += [("E_g", 1.26, dict(Eg=1.26)), ("E_g", 1.45, dict(Eg=1.45))]
    cases += [("eps_r", v, dict(eps_r=v)) for v in (5.0, 6.4, 8.0)]
    cases += [("mh_dos", v, dict(mh_dos=v, mh_tun=v)) for v in (0.1, 0.17, 0.3)]
    for par, val, P in cases:
        for kind, Lc in [("top", 50 * NM), ("edge", 0.0)]:
            for Phi in (4.8, 5.0):
                rows.append(point(kind, Phi, P=P, Lc=Lc or 50 * NM,
                                  extra=dict(param=par, value=val))[0])
    write("material_sensitivity", rows)


def s_thermalbc():
    rows = []
    cases = [("baseline", {})]
    cases += [("kappa", dict(kappa=v)) for v in (0.02, 0.2)]
    cases += [("h_sub", dict(h_sub=v)) for v in (1.35e2, 1.35e4)]
    cases += [("h_contact", dict(h_contact=1.0e4))]
    for tag, P in cases:
        for kind, Lc in [("top", 50 * NM), ("edge", 50 * NM)]:
            dv = m.Device(kind, P=P, Phi_M=5.0, Lc=Lc)
            dv.set_bias(2.0)
            pb = dv.power_balance()
            I = dv.current()[0]
            r = dict(case=tag, **{k: v for k, v in P.items()}, geom=label(kind, Lc), V=2.0,
                     I_A_per_cm=I, R_tot=2.0 / I, Tmax=pb["Tmax"], dT=pb["Tmax"] - 300,
                     P_el=pb["P_el"], P_out=pb["P_out"])
            rows.append(r)
            print("RES", r)
    write("thermal_bc_sensitivity", rows)


def s_mesh():
    rows = []
    for kind, Lc in [("top", 50 * NM), ("edge", 0.0)]:
        for Phi in (4.8, 5.0):
            for mf in (0.5, 1.0, 2.0, 3.0):
                rows.append(point(kind, Phi, Lc=Lc or 50 * NM, mf=mf, extra=dict(mf=mf))[0])
    write("mesh_convergence", rows)


if __name__ == "__main__":
    for s in sys.argv[1:]:
        globals()["s_" + s]()
