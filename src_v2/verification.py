"""
Verification of the v2 model against closed-form results (addresses R2 'verification',
R1.1 / R1.2 implementation concerns). Writes ../data_v2/verification_*.csv.

 V1  Ohmic p-type resistor:   R = L / (q N_A mu_p)                       (drift)
 V2  Robin contact resistor:  R = L/(q N_A mu_p) + 2 V_t/(q v N_A)       (finite interfacial injection)
 V3  Mott-Gurney SCLC:        J -> (9/8) eps mu V^2 / L^3                (Poisson-carrier coupling)
 V4  Schottky diode:          J_0 -> A* T^2 exp(-phi_B/kT), ideality 1   (thermionic emission)
 V5  Energy balance (2D):     I V = P_bulk + P_interface = heat out      (electro-thermal coupling)
"""
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import inse2d as m  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "data_v2")
os.makedirs(OUT, exist_ok=True)
rows = []


def rec(test, quantity, simulated, analytic, note=""):
    err = (simulated - analytic) / analytic
    rows.append(dict(test=test, quantity=quantity, simulated=simulated, analytic=analytic,
                     rel_error=err, note=note))
    print(f"RES {test:4s} {quantity:28s} sim={simulated:.5e} ana={analytic:.5e} err={err:+.2e} {note}")


base = dict(thermal=False, tunneling=False)

# V1 ohmic resistor (1D, per cm^2)
L = 1e-4
dv = m.Device("1d", P=base, L1d=L, contact_bc=("ohmic", "ohmic"), ps1d=1e-6)
dv.set_bias(0.01)
I = dv.current()[0]
rec("V1", "R ohmic (Ohm cm^2)", 0.01 / I, L / (m.q * 1e16 * 50.0))

# V2 Robin resistor: zero-barrier contact (Phi_M = E_V would be degenerate); use the
# ohmic potential BC but a finite interface velocity, implemented by a Phi_M giving
# p_eq = N_A exactly.
P2 = dict(base)
Dd = m.derived(m.DEFAULTS)
Phi_flat = m.DEFAULTS["chi"] + m.DEFAULTS["Eg"] - Dd["Vt"] * math.log(Dd["Nv"] / 1e16)
for v in (1e9, 1e5, 1e4):
    dv = m.Device("1d", P=P2, Phi_M=Phi_flat, L1d=L, ps1d=1e-6, v_override=v)
    dv.set_bias(0.005)
    I = dv.current()[0]
    ana = L / (m.q * 1e16 * 50.0) + 2 * Dd["Vt"] / (m.q * v * 1e16)
    rec("V2", f"R Robin v={v:.0e} (Ohm cm^2)", 0.005 / I, ana)

# V3 Mott-Gurney: nearly intrinsic film, barrier-free injecting contacts (p_eq = N_V).
# Note: with p_eq = 1e16 the same test becomes contact (injection) limited once
# J_MG > q v_R p_eq, which is the physically expected behaviour of a Robin contact.
P3 = dict(base, NA=1e8)
L3 = 1e-4
Phi0 = m.DEFAULTS["chi"] + m.DEFAULTS["Eg"]
dv = m.Device("1d", P=P3, Phi_M=Phi0, L1d=L3, ps1d=2e-7)
mg = []
for V in (1, 2, 5, 10, 20, 40, 80, 160):
    dv.set_bias(V, step=0.2)
    J = dv.current()[0]
    Jmg = 9 / 8 * m.DEFAULTS["eps_r"] * m.eps0 * 50.0 * V ** 2 / L3 ** 3
    mg.append((V, J, Jmg))
    rec("V3", f"J(V={V}) / J_MG", J / Jmg, 1.0, "diffusion raises J at low V")
V = np.array([r[0] for r in mg]); J = np.array([r[1] for r in mg])
slope = np.polyfit(np.log(V[-3:]), np.log(J[-3:]), 1)[0]
rec("V3", "dlnJ/dlnV (40-160 V)", slope, 2.0)
with open(os.path.join(OUT, "verification_mott_gurney.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["V", "J_sim_A_cm2", "J_MottGurney_A_cm2"]); w.writerows(mg)

# V4 Schottky diode (thermionic emission limit: high mobility so v_D >> v_R)
phiB = 0.40
P4 = dict(base, mu_p=1e5, mu_n=1e5)
PhiM4 = m.DEFAULTS["chi"] + m.DEFAULTS["Eg"] - phiB
dv = m.Device("1d", P=P4, Phi_M=PhiM4, L1d=1e-4, contact_bc=("schottky", "ohmic"), ps1d=5e-7)
iv = []
for Vb in np.arange(0.0, 0.301, 0.025):
    dv.set_bias(float(Vb), step=0.025)
    iv.append((float(Vb), dv.current()[0]))
dv.set_bias(0.0, step=0.05)
for Vb in np.arange(-0.025, -0.301, -0.025):
    dv.set_bias(float(Vb), step=0.025)
    iv.append((float(Vb), dv.current()[0]))
iv = np.array(sorted(iv))
# forward branch = the polarity on which |J| grows exponentially
pos = iv[iv[:, 0] > 0]; neg = iv[iv[:, 0] < 0]
fwd = pos if abs(pos[-1, 1]) > abs(neg[0, 1]) else neg[::-1]
rev = neg if fwd is pos else pos
sel = (np.abs(fwd[:, 0]) >= 0.075) & (np.abs(fwd[:, 0]) <= 0.2)
Vf = np.abs(fwd[sel, 0]); Jf = np.abs(fwd[sel, 1])
p = np.polyfit(Vf, np.log(Jf), 1)
n_ideal = 1 / (p[0] * Dd["Vt"])
J0 = math.exp(p[1])
J0_te = Dd["A_p"] * 300 ** 2 * math.exp(-phiB / Dd["Vt"])
rec("V4", "ideality factor", n_ideal, 1.0)
rec("V4", "J0 (A/cm^2)", J0, J0_te, "A* T^2 exp(-phiB/kT)")
rec("V4", "|J_rev(0.3 V)| (A/cm^2)", abs(rev[-1, 1]) if fwd is neg else abs(rev[0, 1]), J0_te, "saturation")
with open(os.path.join(OUT, "verification_schottky_diode.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["V_drain", "J_A_cm2"]); w.writerows(iv.tolist())

# V5 energy balance in 2D electro-thermal devices
for kind in ("top", "edge"):
    for PhiM, Vb in ((4.8, 0.5), (5.0, 1.0)):
        dv = m.Device(kind, Phi_M=PhiM)
        dv.set_bias(Vb)
        pb = dv.power_balance()
        rec("V5", f"{kind} {PhiM}eV {Vb}V (P_bulk+P_if)/IV", (pb["P_bulk"] + pb["P_interface"]), pb["P_el"])
        rec("V5", f"{kind} {PhiM}eV {Vb}V P_out/P_bulk", pb["P_out"], pb["P_bulk"])
        rec("V5", f"{kind} {PhiM}eV {Vb}V I_src/I_drn", -dv.current_source(), dv.current()[0])

with open(os.path.join(OUT, "verification_summary.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
