"""
v3 verification suite -> ../data_v3/verification_*.csv
 V1  ohmic resistor                         R = L/(q N_A mu)
 V2  finite-velocity (Robin) contact        R = L/(q N_A mu) + 2V_t/(q v N_A)
 V3  Mott-Gurney SCLC (Boltzmann and FD)    J -> 9/8 eps mu V^2/L^3
 V4  Schottky diode TE at 300 and 350 K     J_0 = A* T^2 exp(-phi_B/kT), n = 1
 V5  2D energy balance (all heat terms), top/edge, gated/ungated
 V6  nonlocal tunnelling vs independent semi-analytic Tsu-Esaki/WKB (depletion approx.)
 V7  Fermi-Dirac degenerate contact at equilibrium: I = 0, p_s = N_V F_1/2(-phi_B/kT)
 V8  transient heating: energy conservation and approach to the DC state
 V9  Gauss law with a gate: sum of contact charges + semiconductor charge = 0
"""
import csv, math, os, sys
import numpy as np
from scipy.integrate import quad

sys.path.insert(0, os.path.dirname(__file__))
import tcad2d as m  # noqa
import devsim as ds  # noqa

OUT = os.path.join(os.path.dirname(__file__), "..", "data_v3")
os.makedirs(OUT, exist_ok=True)
rows = []


def rec(test, quantity, sim, ana, note=""):
    err = (sim - ana) / ana if ana != 0 else sim
    rows.append(dict(test=test, quantity=quantity, simulated=sim, analytic=ana, rel_error=err, note=note))
    print(f"RES {test:4s} {quantity:40s} sim={sim:.5e} ref={ana:.5e} err={err:+.2e} {note}", flush=True)


BZ = dict(thermal=False, tunneling="off", statistics="boltzmann", band_from_thickness=False,
          Eg=1.30, chi=3.70)
Dd = m.derived(dict(m.DEFAULTS, **BZ), 1.30)
Ev = 3.70 + 1.30

# V1
L = 1e-4
dv = m.Device("1d", P=BZ, L1d=L, contact_bc=("ohmic", "ohmic"), ps1d=1e-6)
dv.set_bias(0.01)
rec("V1", "R ohmic (Ohm cm^2)", 0.01 / dv.current()[0], L / (m.q * 1e16 * 50.0))

# V2
Phi_flat = Ev - Dd["Vt"] * math.log(Dd["Nv"] / 1e16)
for v in (1e9, 1e5, 1e4):
    dv = m.Device("1d", P=BZ, Phi_M=Phi_flat, L1d=L, ps1d=1e-6, v_override=v)
    dv.set_bias(0.005)
    rec("V2", f"R Robin v={v:.0e} (Ohm cm^2)", 0.005 / dv.current()[0],
        L / (m.q * 1e16 * 50.0) + 2 * Dd["Vt"] / (m.q * v * 1e16))

# V3
mg = []
for stat in ("boltzmann", "fd"):
    P3 = dict(BZ, NA=1e8, statistics=stat)
    dv = m.Device("1d", P=P3, Phi_M=Ev, L1d=1e-4, ps1d=2e-7)
    for V in (5, 20, 80, 160):
        dv.set_bias(V, step=0.2)
        J = dv.current()[0]
        Jmg = 9 / 8 * 6.4 * m.eps0 * 50.0 * V ** 2 / 1e-12
        mg.append((stat, V, J, Jmg))
        rec("V3", f"{stat} J/J_MG at {V} V", J / Jmg, 1.0)
with open(os.path.join(OUT, "verification_mott_gurney.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["statistics", "V", "J_sim", "J_MG"]); w.writerows(mg)

# V4 (300 K and 350 K)
phiB = 0.40
for T0 in (300.0, 350.0):
    P4 = dict(BZ, mu_p=1e5, mu_n=1e5, T0=T0)
    D4 = m.derived(dict(m.DEFAULTS, **P4), 1.30, T=T0)
    dv = m.Device("1d", P=P4, Phi_M=Ev - phiB, L1d=1e-4, contact_bc=("schottky", "ohmic"), ps1d=5e-7)
    iv = []
    for Vb in np.arange(0.025, 0.301, 0.025):
        dv.set_bias(float(Vb), step=0.025); iv.append((float(Vb), dv.current()[0]))
    dv.set_bias(0.0)
    for Vb in np.arange(-0.025, -0.301, -0.025):
        dv.set_bias(float(Vb), step=0.025); iv.append((float(Vb), dv.current()[0]))
    iv = np.array(sorted(iv))
    pos, neg = iv[iv[:, 0] > 0], iv[iv[:, 0] < 0]
    fwd, rev = (pos, neg) if abs(pos[-1, 1]) > abs(neg[0, 1]) else (neg[::-1], pos)
    sel = (np.abs(fwd[:, 0]) >= 0.075) & (np.abs(fwd[:, 0]) <= 0.2)
    pfit = np.polyfit(np.abs(fwd[sel, 0]), np.log(np.abs(fwd[sel, 1])), 1)
    J0te = D4["A_p"] * T0 ** 2 * math.exp(-phiB / D4["Vt"])
    rec("V4", f"ideality ({T0:.0f} K)", 1 / (pfit[0] * D4["Vt"]), 1.0)
    jr = abs(rev[0, 1]) if rev is neg else abs(rev[-1, 1])
    rec("V4", f"|J_rev(0.3V)| ({T0:.0f} K)", jr, J0te, "A*T^2 exp(-phiB/kT)")
    with open(os.path.join(OUT, f"verification_schottky_{T0:.0f}K.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["V", "J"]); w.writerows(iv.tolist())

# V5 energy balance, full v3 physics
for kind in ("top", "edge"):
    for gate in (None, dict(t_ox=10e-7)):
        for Phi, Vb in ((4.8, 0.05), (5.0, 0.5)):
            dv = m.Device(kind, Phi_M=Phi, gate=gate)
            if gate:
                dv.set_gate(-2.0)
            dv.set_bias(Vb)
            pb = dv.power_balance()
            tag = f"{kind}{'-gated' if gate else ''} {Phi}eV {Vb}V"
            rec("V5", f"{tag} P_out/IV", pb["P_out"], pb["P_el"])
            rec("V5", f"{tag} -I_src/I_drn", -dv.current_source(), dv.current()[0])

# V6 nonlocal tunnelling vs independent semi-analytic calculation (1D, heavy doping)
NA6, phiB6, Vr = 1e18, 0.40, -0.3
P6 = dict(BZ, NA=NA6, tunneling="nonlocal", mu_p=1e4, mu_n=1e4, tun_maxlen=30e-7)
dv = m.Device("1d", P=P6, Phi_M=Ev - phiB6, L1d=0.2e-4, contact_bc=("schottky", "ohmic"), ps1d=0.1e-7)
dv.set_bias(Vr, step=0.02)
I_sim_tun = dv.tunnel_current()["source"]
# independent: band profile from the simulation's own potential (tests the integral and
# bookkeeping), and from the depletion approximation (tests the coupled solution)
x = dv.x; psi = dv.nodev("Potential"); order = np.argsort(x); x, psi = x[order], psi[order]
Vt = Dd["Vt"]
D6 = m.derived(dict(m.DEFAULTS, **P6), 1.30)
Bh = psi - dv.Cv  # hole band edge (eps units)
phip = dv._phi()["phip"][order]


def ts_integral(xx, B, eFm, eFs_fun, mt):
    kap0 = math.sqrt(2 * mt * m.m0 * m.q) / m.hbar * 1e-2
    B0 = B[0]; Bmin = np.minimum.accumulate(B)

    def integrand(eps):
        over = np.clip(B - eps, 0, None)
        i_t = np.argmax(B <= eps) if np.any(B <= eps) else len(B) - 1
        integ = np.trapezoid(np.sqrt(over[: i_t + 1]), xx[: i_t + 1])
        Tr = math.exp(-2 * kap0 * integ)
        return Tr * (np.logaddexp(0, (eFm - eps) / Vt) - np.logaddexp(0, (eFs_fun(i_t) - eps) / Vt))
    lo = Bmin[np.searchsorted(xx, 30e-7) - 1]
    val, _ = quad(integrand, lo, B0, limit=400)
    return D6["A_p"] * 300.0 / m.kq * val


J_semi_sim = ts_integral(x, Bh, 0.0 - 0.0, lambda i: phip[i], 0.17)
rec("V6", "J_tun: model vs quad on same profile", I_sim_tun, J_semi_sim, "integration/bookkeeping")
# depletion approximation: parabolic band bending with the simulated total bending
eps_s = 6.4 * m.eps0
dB = Bh[0] - Bh[-1]
W = math.sqrt(2 * eps_s * abs(dB) / (m.q * NA6))
xx = np.linspace(0, x[-1], 4000)
Bdep = Bh[-1] + np.sign(dB) * np.where(xx < W, (m.q * NA6 / (2 * eps_s)) * (W - xx) ** 2, 0.0)
phip_i = np.interp(xx, x, phip)
J_dep = ts_integral(xx, Bdep, 0.0, lambda i: phip_i[i], 0.17)
rec("V6", "J_tun: model vs depletion-approx. quad", I_sim_tun, J_dep, "full Poisson vs depletion")
rec("V6", "J_tun / J_total", I_sim_tun / dv.current_source(), 1.0, "tunnelling share (info)")

# V7 FD degenerate contact equilibrium
phiB7 = -0.30
P7 = dict(BZ, statistics="fd")
dv = m.Device("1d", P=P7, Phi_M=Ev - phiB7, L1d=1e-5, ps1d=0.5e-7)
I0 = dv.current()[0]
ps = dv.nodev("Holes")[np.argmin(dv.x)]
f12 = 2 / math.sqrt(math.pi) * quad(lambda e: math.sqrt(e) / (1 + math.exp(e + phiB7 / Dd["Vt"])), 0, 80)[0]
rec("V7", "p_s / (N_V F1/2(-phiB/kT)) [exact F]", ps, Dd["Nv"] * f12)
rows.append(dict(test="V7", quantity="|I(V=0)| (A/cm^2)", simulated=abs(I0), analytic=0.0, rel_error=abs(I0), note=""))
print("RES V7 |I(V=0)|", abs(I0))
phi_eq = dv._phi()["phip"]
rows.append(dict(test="V7", quantity="max|phi_p| at equilibrium (V)", simulated=float(np.abs(phi_eq).max()), analytic=0.0,
                 rel_error=float(np.abs(phi_eq).max()), note="uniform quasi-Fermi level"))
print("RES V7 max|phi_p|", float(np.abs(phi_eq).max()))

# V8 transient energy conservation (gated edge, 5.0 eV, step 0 -> 1 V at V_G=-2 V)
dv = m.Device("edge", Phi_M=5.0, gate=dict(t_ox=10e-7), P=dict(tunneling="off"))
dv.set_gate(-2.0)
dv.set_bias(0.0)
T_before = {r: dv.nodev("T", r).copy() for r in (m.SEMI, m.OX)}
res = dv.transient_step(1.0, 1e-6, nsteps=60)
t, Tm, Pel, Pout = res.T
E_in = np.trapezoid(np.concatenate([[0], Pel]), np.concatenate([[0], t]))
E_out = np.trapezoid(np.concatenate([[0], Pout]), np.concatenate([[0], t]))
E_st = 0.0
for r, cvn in ((m.SEMI, "cv"), (m.OX, "cv_ox")):
    E_st += float(np.sum(dv.P[cvn] * (dv.nodev("T", r) - T_before[r]) * dv.nodev("NodeVolume", r)))
rec("V8", "(E_out + E_stored)/E_in (0-1 us)", E_out + E_st, E_in, "bdf1, 60 log steps")
Tdc = dv.power_balance()["Tmax"]
dv2 = m.Device("edge", Phi_M=5.0, gate=dict(t_ox=10e-7), P=dict(tunneling="off"))
dv2.set_gate(-2.0); dv2.set_bias(1.0)
rec("V8", "T_max(1 us) - 300 vs DC", Tm[-1] - 300, dv2.power_balance()["Tmax"] - 300)
with open(os.path.join(OUT, "verification_transient.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["t_s", "Tmax_K", "P_el", "P_out"]); w.writerows(res.tolist())

# V9 Gauss law
dv = m.Device("edge", Phi_M=4.8, gate=dict(t_ox=10e-7), P=dict(thermal=False))
dv.set_gate(-2.0); dv.set_bias(0.2)
Qc = sum(ds.get_contact_charge(device=dv.name, contact=c, equation="PE") for c in ("source", "drain", "gate"))
rho = m.q * (dv.nodev("Holes") - dv.nodev("Electrons") + (dv.P["ND"] - dv.P["NA"]))
Qs = float(np.sum(rho * dv.nodev("NodeVolume")))
Qg = ds.get_contact_charge(device=dv.name, contact="gate", equation="PE")
rec("V9", "(sum contact charge + semi charge)/|Q_gate|", (Qc + Qs) / abs(Qg) + 1.0, 1.0, "should be 1")

with open(os.path.join(OUT, "verification_summary.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("RES done")
