"""
tcad2d.py -- v3 DEVSIM model for metal / few-layer 2D-semiconductor contacts (InSe, MoS2).

Extends the v2 model (commit 53bbd08, src_v2/inse2d.py, removed from the tree) with:
  * anisotropic (in-plane / out-of-plane) mobility and thermal conductivity
  * optional back gate (oxide region + gate electrode, flat band at V_G = 0)
  * Fermi-Dirac statistics (Nilsson inverse / Bednarczyk F_1/2), degenerate contacts allowed
  * temperature-dependent statistics and thermionic injection (non-isothermal option)
  * interface (contact) Joule heat and Peltier heat (bulk -T J.grad(P), junction J T P)
  * nonlocal WKB / Tsu-Esaki tunnelling along contact normals (Ieong-type), or local Gamma(F)
  * separate pinning factors for top and edge contacts
  * transient lattice heating (heat capacity)
  * thickness-dependent band gap / affinity for InSe (layer number)
  * material presets: InSe, MoS2
Units: cm, V, A, s, K. 2D results per unit width (A/cm, Ohm*cm).
"""
import math
import os

os.environ.setdefault(
    "DEVSIM_MATH_LIBS",
    "/usr/lib/x86_64-linux-gnu/liblapack.so.3:/usr/lib/x86_64-linux-gnu/libblas.so.3",
)
import numpy as np
import devsim as ds

ds.set_parameter(name="debug_level", value="error")

q = 1.602176634e-19
kB = 1.380649e-23
kq = kB / q
h_planck = 6.62607015e-34
hbar = h_planck / (2 * math.pi)
m0 = 9.1093837015e-31
eps0 = 8.8541878128e-14  # F/cm
SEMI, OX = "semi", "ox"
SQPI = math.sqrt(math.pi)

# ---------------------------------------------------------------------------
# Materials
# ---------------------------------------------------------------------------
INSE = dict(
    material="InSe",
    band_from_thickness=True,
    Eg_bulk=1.26, chi_bulk=3.70,        # bulk gamma-InSe; E_V,bulk = 4.96 eV
    Eg_A=2.78, Eg_alpha=1.70,           # dEg = A N^-alpha (fit to PL: 1.69 eV @3L, ~1.44 eV @5L)
    dEv_frac=0.32,                      # share of dEg taken by the valence band (~ m_e/(m_e+m_h))
    layer_nm=0.83,
    Eg=1.39, chi=3.61,                  # used if band_from_thickness=False (5 nm values)
    me_dos=0.07, mh_dos=0.17, me_tun=0.07, mh_tun=0.17,
    eps_r=6.4, mu_n=1000.0, mu_p=50.0, mu_T_exp=1.5,
    mu_perp_ratio=0.01,                 # out-of-plane / in-plane mobility (bulk InSe ~0.01; VFET 0.002)
    kappa_par=0.085, kappa_perp=0.0076, # W/cmK, TDTR (Rai et al. 2021)
    cv=1.43,                            # J/cm^3K (Dulong-Petit, 5.55 g/cm^3)
)
MOS2 = dict(
    material="MoS2",
    band_from_thickness=False,
    Eg=1.30, chi=4.00, layer_nm=0.65,
    me_dos=0.55, mh_dos=0.65, me_tun=0.50, mh_tun=0.60,
    eps_r=7.0, mu_n=35.0, mu_p=35.0, mu_T_exp=1.6,
    mu_perp_ratio=0.01, kappa_par=0.5, kappa_perp=0.02, cv=1.9,
)

DEFAULTS = dict(
    INSE,
    T0=300.0, tau=1e-9, NA=1e16, ND=0.0,
    h_sub=1.35e3,            # W/cm^2K, ungated substrate Robin coefficient
    R_TBR=1e-8 * 1e4,        # K cm^2/W  (1e-8 m^2K/W), used at the gate heat sink when gated
    h_contact=0.0,
    S_top=1.0, S_edge=1.0, S_1d=1.0, Phi_CNL=4.35,
    tunneling="nonlocal",    # "off" | "local" | "nonlocal"
    tun_maxlen=40e-7,        # cm
    T_vdW_top=None, d_vdW=0.30e-7, Phi_vdW=2.0,
    thermal=True, statistics="fd", nonisothermal=True, peltier=True, interface_heat=True,
    eps_ox=3.9, kappa_ox=0.014, cv_ox=1.6,
)


def band(P, t_cm):
    """Band gap / affinity, optionally from the layer number of the film."""
    if P.get("band_from_thickness") and t_cm:
        N = t_cm / (P["layer_nm"] * 1e-7)
        dEg = P["Eg_A"] * N ** (-P["Eg_alpha"])
        Eg = P["Eg_bulk"] + dEg
        chi = P["chi_bulk"] - (1 - P["dEv_frac"]) * dEg
        return Eg, chi, N
    return P["Eg"], P["chi"], (t_cm / (P["layer_nm"] * 1e-7) if t_cm else None)


def derived(P, Eg, T=None):
    T = T or P["T0"]
    Vt = kq * T
    fac = 2 * (2 * math.pi * kB * T / h_planck ** 2) ** 1.5 * 1e-6
    Nc = fac * (P["me_dos"] * m0) ** 1.5
    Nv = fac * (P["mh_dos"] * m0) ** 1.5
    ni = math.sqrt(Nc * Nv) * math.exp(-Eg / (2 * Vt))
    A_p, A_n = 120.17 * P["mh_dos"], 120.17 * P["me_dos"]
    return dict(Vt=Vt, Nc=Nc, Nv=Nv, ni=ni, A_p=A_p, A_n=A_n,
                vR_p=A_p * T ** 2 / (q * Nv), vR_n=A_n * T ** 2 / (q * Nc))


def F12(eta):
    """Normalised Fermi-Dirac integral F_1/2 (Bednarczyk & Bednarczyk 1978, <0.4 %)."""
    xi = 0.75 * SQPI * (eta ** 4 + 50 + 33.6 * eta * (1 - 0.68 * math.exp(-0.17 * (eta + 1) ** 2))) ** (-0.375)
    return 1.0 / (math.exp(-eta) + xi)


def F12_expr(e):
    return (f"(1/(exp(-({e}))+0.75*{SQPI}*pow(pow({e},4)+50+33.6*({e})*(1-0.68*exp(-0.17*pow(({e})+1,2))),-0.375)))")


def eta_inv(u):
    """Nilsson (1973) inverse of F_1/2 (python, for checks)."""
    v = (0.75 * SQPI * u) ** (2 / 3)
    a = -0.5 + 0.5 * (u - 1) if abs(u - 1) < 1e-4 else math.log(u) / (1 - u * u)
    return a + v / (1 + (0.24 + 1.08 * v) ** -2)


def delta_expr(u):
    """Degeneracy correction Delta(u) = eta(u) - ln(u), u = n/N_C (or p/N_V)."""
    v = f"pow(0.75*{SQPI}*({u}),0.6666666666666666)"
    a = f"ifelse(abs(({u})-1)<1e-4,-0.5+0.5*(({u})-1),log({u})/(1-({u})*({u})))"
    return f"(({a})+({v})/(1+pow(0.24+1.08*({v}),-2))-log({u}))"


def phi_eff(P, Phi_M, S):
    return P["Phi_CNL"] + S * (Phi_M - P["Phi_CNL"])


def wkb_vdw(P):
    if P["T_vdW_top"] is not None:
        return P["T_vdW_top"]
    if P["d_vdW"] <= 0:
        return 1.0
    k = math.sqrt(2 * m0 * P["Phi_vdW"] * q) / hbar
    return math.exp(-2 * k * P["d_vdW"] * 1e-2)


def gamma_tfe(F_Vcm, phiB_eV, m_tun, T=300.0):
    """Local TFE enhancement (triangular-barrier WKB), kept for comparison."""
    if phiB_eV <= 0 or F_Vcm <= 0:
        return 0.0
    kT = kB * T
    F = F_Vcm * 100.0
    b = (4.0 / 3.0) * math.sqrt(2 * m_tun * m0) / (q * hbar * F)
    E = np.linspace(0.0, phiB_eV * q, 4001)
    return float(np.trapezoid(np.exp(E / kT - b * E ** 1.5), E) / kT)


# ---------------------------------------------------------------------------
# model-building helpers
# ---------------------------------------------------------------------------
def _nm(d, r, name, expr, deps=()):
    ds.node_model(device=d, region=r, name=name, equation=expr)
    for v in deps:
        ds.node_model(device=d, region=r, name=f"{name}:{v}", equation=f"diff({expr},{v})")


def _em(d, r, name, expr, deps=()):
    ds.edge_model(device=d, region=r, name=name, equation=expr)
    for v in deps:
        for k in ("n0", "n1"):
            ds.edge_model(device=d, region=r, name=f"{name}:{v}@{k}", equation=f"diff({expr},{v}@{k})")


def _cnm(d, c, name, expr, deps=()):
    ds.contact_node_model(device=d, contact=c, name=name, equation=expr)
    for v in deps:
        ds.contact_node_model(device=d, contact=c, name=f"{name}:{v}", equation=f"diff({expr},{v})")


def _im(d, i, name, expr, deps=()):
    ds.interface_model(device=d, interface=i, name=name, equation=expr)
    for v in deps:
        for r in ("r0", "r1"):
            ds.interface_model(device=d, interface=i, name=f"{name}:{v}@{r}", equation=f"diff({expr},{v}@{r})")


class Device:
    def __init__(self, kind, P=None, Phi_M=5.0, name="dev", mf=1.0, Lc=50e-7, Lch=50e-7,
                 t=5e-7, L1d=1e-4, contact_bc=("schottky", "schottky"), ps1d=1e-7,
                 v_override=None, gate=None):
        ds.reset_devsim()
        self.P = dict(DEFAULTS)
        if P:
            self.P.update(P)
        P = self.P
        self.kind, self.name, self.Phi_M = kind, name, Phi_M
        self.gate = gate  # dict(t_ox=cm, eps_ox=..., h_gate=...) or None
        self.t = t if kind != "1d" else None
        self.Eg, self.chi, self.Nlayers = band(P, self.t)
        self.D = derived(P, self.Eg)
        S = {"top": P["S_top"], "edge": P["S_edge"], "1d": P["S_1d"]}[kind]
        Pe = phi_eff(P, Phi_M, S)
        Ev = self.chi + self.Eg
        self.fd = P["statistics"] == "fd"
        self.phiBp_raw = Ev - Pe
        self.phiBp = self.phiBp_raw if self.fd else max(self.phiBp_raw, 0.0)
        self.Phi_used = Ev - self.phiBp
        self.phiBn = self.Eg - self.phiBp
        self.contact_bc = dict(zip(("source", "drain"), contact_bc))
        self.v_override = v_override
        self.T_vdW = wkb_vdw(P) if kind == "top" else 1.0
        self.thermal = P["thermal"] and kind != "1d"
        self.noniso = P["nonisothermal"] and self.thermal
        self.bias = 0.0
        self.vg = 0.0
        if kind == "top":
            self._mesh_top(Lc, Lch, t, mf)
        elif kind == "edge":
            self._mesh_edge(Lch, t, mf)
        else:
            self._mesh_1d(L1d, ps1d)
        self._build()

    # ------------------------------------------------------------------ meshes
    def _gate_lines(self, mf):
        g = self.gate
        tox = g["t_ox"]
        return [(self.t + 0.5 * tox, max(tox / 8, 0.5e-7)), (self.t + tox, max(tox / 8, 0.5e-7)),
                (self.t + tox + 1e-7, 0.5e-7)]

    def _mesh_top(self, Lc, Lch, t, mf):
        n = self.name
        ds.create_2d_mesh(mesh=n)
        L = 2 * Lc + Lch
        po = min(max(2e-7, Lc / 25), 10e-7) if Lc > 8e-7 else min(2e-7, Lc / 4)
        xs = [(0.0, po), (Lc, 0.25e-7), (Lc + 0.5 * Lch, 2e-7), (Lc + Lch, 0.25e-7), (L, po)]
        for x, ps in xs:
            ds.add_2d_mesh_line(mesh=n, dir="x", pos=x, ps=ps / mf)
        ys = [(-1e-7, 0.5e-7), (0.0, 0.1e-7), (0.5 * t, 0.5e-7), (t, 0.2e-7)]
        if self.gate:
            ys += self._gate_lines(mf)
        else:
            ys += [(t + 1e-7, 0.5e-7)]
        for y, ps in ys:
            ds.add_2d_mesh_line(mesh=n, dir="y", pos=y, ps=ps / mf)
        ds.add_2d_region(mesh=n, material="semi", region=SEMI, xl=0, xh=L, yl=0, yh=t)
        ds.add_2d_region(mesh=n, material="metal", region="m_src", xl=0, xh=Lc, yl=-1e-7, yh=0)
        ds.add_2d_region(mesh=n, material="metal", region="m_drn", xl=Lc + Lch, xh=L, yl=-1e-7, yh=0)
        b = 1e-10
        ds.add_2d_contact(mesh=n, name="source", material="metal", region=SEMI, xl=0, xh=Lc, yl=0, yh=0, bloat=b)
        ds.add_2d_contact(mesh=n, name="drain", material="metal", region=SEMI, xl=Lc + Lch, xh=L, yl=0, yh=0, bloat=b)
        self._bottom(n, 0, L, t)
        ds.finalize_mesh(mesh=n)
        ds.create_device(mesh=n, device=n)
        self.geo = dict(x_ch0=Lc, x_ch1=Lc + Lch, Lc=Lc, Lch=Lch, t=t, geom="top", L=L)

    def _mesh_edge(self, Lch, t, mf):
        n = self.name
        ds.create_2d_mesh(mesh=n)
        for x, ps in [(-1e-7, 0.5e-7), (0.0, 0.25e-7), (0.5 * Lch, 2e-7), (Lch, 0.25e-7), (Lch + 1e-7, 0.5e-7)]:
            ds.add_2d_mesh_line(mesh=n, dir="x", pos=x, ps=ps / mf)
        ys = [(0.0, 0.2e-7), (0.5 * t, 0.5e-7), (t, 0.2e-7)]
        ys += self._gate_lines(mf) if self.gate else [(t + 1e-7, 0.5e-7)]
        for y, ps in ys:
            ds.add_2d_mesh_line(mesh=n, dir="y", pos=y, ps=ps / mf)
        ds.add_2d_region(mesh=n, material="semi", region=SEMI, xl=0, xh=Lch, yl=0, yh=t)
        ds.add_2d_region(mesh=n, material="metal", region="m_src", xl=-1e-7, xh=0, yl=0, yh=t)
        ds.add_2d_region(mesh=n, material="metal", region="m_drn", xl=Lch, xh=Lch + 1e-7, yl=0, yh=t)
        b = 1e-10
        ds.add_2d_contact(mesh=n, name="source", material="metal", region=SEMI, xl=0, xh=0, yl=0, yh=t, bloat=b)
        ds.add_2d_contact(mesh=n, name="drain", material="metal", region=SEMI, xl=Lch, xh=Lch, yl=0, yh=t, bloat=b)
        self._bottom(n, 0, Lch, t)
        ds.finalize_mesh(mesh=n)
        ds.create_device(mesh=n, device=n)
        self.geo = dict(x_ch0=0.0, x_ch1=Lch, Lc=0.0, Lch=Lch, t=t, geom="edge", L=Lch)

    def _bottom(self, n, x0, x1, t):
        b = 1e-10
        if self.gate:
            tox = self.gate["t_ox"]
            ds.add_2d_region(mesh=n, material="oxide", region=OX, xl=x0, xh=x1, yl=t, yh=t + tox)
            ds.add_2d_region(mesh=n, material="metal", region="m_gate", xl=x0, xh=x1, yl=t + tox, yh=t + tox + 1e-7)
            ds.add_2d_interface(mesh=n, name="semi_ox", region0=SEMI, region1=OX)
            ds.add_2d_contact(mesh=n, name="gate", material="metal", region=OX, xl=x0, xh=x1,
                              yl=t + tox, yh=t + tox, bloat=b)
        else:
            ds.add_2d_region(mesh=n, material="metal", region="m_sub", xl=x0, xh=x1, yl=t, yh=t + 1e-7)
            ds.add_2d_contact(mesh=n, name="substrate", material="metal", region=SEMI, xl=x0, xh=x1,
                              yl=t, yh=t, bloat=b)

    def _mesh_1d(self, L, ps):
        n = self.name
        ds.create_1d_mesh(mesh=n)
        ds.add_1d_mesh_line(mesh=n, pos=0.0, ps=ps, tag="l")
        ds.add_1d_mesh_line(mesh=n, pos=0.5 * L, ps=ps)
        ds.add_1d_mesh_line(mesh=n, pos=L, ps=ps, tag="h")
        ds.add_1d_contact(mesh=n, name="source", tag="l", material="metal")
        ds.add_1d_contact(mesh=n, name="drain", tag="h", material="metal")
        ds.add_1d_region(mesh=n, material="semi", region=SEMI, tag1="l", tag2="h")
        ds.finalize_mesh(mesh=n)
        ds.create_device(mesh=n, device=n)
        self.geo = dict(x_ch0=0.0, x_ch1=L, Lc=0.0, Lch=L, t=None, geom="1d", L=L)

    # ------------------------------------------------------------------ expressions
    def Tn(self):
        return "T" if self.noniso else "T0"

    def Te(self):
        return "(0.5*(T@n0+T@n1))" if self.noniso else "T0"

    def _node_stat(self, carrier, at=""):
        """(N_eff expression, Delta expression) at a node; at='' (node) or '@n0'/'@n1'."""
        T = ("T" + at) if self.noniso else "T0"
        N = f"(Nv300*pow({T}/300,1.5))" if carrier == "p" else f"(Nc300*pow({T}/300,1.5))"
        var = ("Holes" if carrier == "p" else "Electrons") + at
        dl = delta_expr(f"{var}/{N}") if self.fd else "0"
        return N, dl, var, T

    def eta(self, carrier, at=""):
        N, dl, var, T = self._node_stat(carrier, at)
        return f"(log({var}/{N})+{dl})"

    # ------------------------------------------------------------------ build
    def _build(self):
        d, P, D = self.name, self.P, self.D
        two_d = self.kind != "1d"
        par = dict(q=q, kq=kq, Vt300=kq * P["T0"], T0=P["T0"], eps=P["eps_r"] * eps0,
                   ni=D["ni"], Nc300=D["Nc"], Nv300=D["Nv"], Eg=self.Eg,
                   NetDoping=P["ND"] - P["NA"], mu_n0=P["mu_n"], mu_p0=P["mu_p"], aT=P["mu_T_exp"],
                   tau=P["tau"], kpar=P["kappa_par"], kperp=P["kappa_perp"], rmu=P["mu_perp_ratio"],
                   h_sub=P["h_sub"], h_ct=P["h_contact"], cv=P["cv"], A_p=D["A_p"], A_n=D["A_n"],
                   phiBp=self.phiBp, phiBn=self.phiBn)
        for k, v in par.items():
            ds.set_parameter(device=d, name=k, value=v)
        regions = [SEMI] + ([OX] if self.gate else [])
        for s in ("Potential", "Electrons", "Holes", "Gp", "Gn", "Gtp_a", "Gtp_b", "Gtn_a", "Gtn_b", "Htun"):
            ds.node_solution(device=d, region=SEMI, name=s)
            ds.edge_from_node_model(device=d, region=SEMI, node_model=s)
        nn = len(ds.get_node_model_values(device=d, region=SEMI, name="x"))
        self.nn = nn
        Vt = D["Vt"]
        self.psi_bi = self.chi - self.Phi_used + self.Eg / 2 + 0.5 * Vt * math.log(D["Nc"] / D["Nv"])
        self.Cv = Vt * math.log(D["ni"] / D["Nv"])  # E_V = -psi + Cv  (eV)
        Nnet = P["NA"] - P["ND"]
        psi_bulk = -Vt * math.log(Nnet / D["ni"]) if Nnet > 0 else Vt * math.log(-Nnet / D["ni"])
        self.psi_bulk = psi_bulk
        ds.set_node_values(device=d, region=SEMI, name="Potential", values=[psi_bulk] * nn)
        for s in ("Gp", "Gn", "Gtp_a", "Gtp_b", "Gtn_a", "Gtn_b", "Htun"):
            ds.set_node_values(device=d, region=SEMI, name=s, values=[0.0] * nn)
        # ---- Poisson (equilibrium first)
        _em(d, SEMI, "EField", "(Potential@n0-Potential@n1)*EdgeInverseLength", ["Potential"])
        _em(d, SEMI, "DField", "eps*EField", ["Potential"])
        ds.set_parameter(device=d, name="Cv", value=self.Cv)
        if self.fd:  # equilibrium (E_F = 0) Fermi-Dirac densities, consistent with the contacts
            _nm(d, SEMI, "p_i", f"Nv300*{F12_expr('(Cv-Potential)/Vt300')}", ["Potential"])
            _nm(d, SEMI, "n_i", f"Nc300*{F12_expr('(Potential-Cv-Eg)/Vt300')}", ["Potential"])
        else:
            _nm(d, SEMI, "p_i", "ni*exp(-Potential/Vt300)", ["Potential"])
            _nm(d, SEMI, "n_i", "ni*exp(Potential/Vt300)", ["Potential"])
        _nm(d, SEMI, "RhoEq", "-q*(p_i-n_i+NetDoping)", ["Potential"])
        ds.equation(device=d, region=SEMI, name="PE", variable_name="Potential",
                    node_model="RhoEq", edge_model="DField", variable_update="log_damp")
        if self.gate:
            g = self.gate
            ds.set_parameter(device=d, region=OX, name="eps", value=g.get("eps_ox", P["eps_ox"]) * eps0)
            ds.node_solution(device=d, region=OX, name="Potential")
            ds.edge_from_node_model(device=d, region=OX, node_model="Potential")
            nox = len(ds.get_node_model_values(device=d, region=OX, name="x"))
            ds.set_node_values(device=d, region=OX, name="Potential", values=[psi_bulk] * nox)
            _em(d, OX, "EField", "(Potential@n0-Potential@n1)*EdgeInverseLength", ["Potential"])
            _em(d, OX, "DField", "eps*EField", ["Potential"])
            ds.equation(device=d, region=OX, name="PE", variable_name="Potential",
                        edge_model="DField", variable_update="default")
            _im(d, "semi_ox", "contPot", "Potential@r0-Potential@r1", ["Potential"])
            ds.interface_equation(device=d, interface="semi_ox", name="PE",
                                  interface_model="contPot", type="continuous")
            ds.set_parameter(device=d, name="gate_bias", value=0.0)
            ds.set_parameter(device=d, name="gate_psi", value=psi_bulk + g.get("V_FB_shift", 0.0))
            ds.contact_node_model(device=d, contact="gate", name="gate_pbc",
                                  equation="Potential-gate_bias-gate_psi")
            ds.contact_node_model(device=d, contact="gate", name="gate_pbc:Potential", equation="1")
            ds.contact_equation(device=d, contact="gate", name="PE", node_model="gate_pbc",
                                edge_charge_model="DField")
        for c in ("source", "drain"):
            ds.set_parameter(device=d, name=f"{c}_bias", value=0.0)
            psi_c = psi_bulk if self.contact_bc[c] == "ohmic" else self.psi_bi
            ds.set_parameter(device=d, name=f"{c}_psi", value=psi_c)
            _cnm(d, c, f"{c}_pbc", f"Potential-{c}_bias-{c}_psi", ["Potential"])
            ds.contact_equation(device=d, contact=c, name="PE", node_model=f"{c}_pbc", edge_charge_model="DField")
        self._solve(rel_err=1e-10, it=100)
        # ---- carriers
        ds.set_node_values(device=d, region=SEMI, name="Holes", init_from="p_i")
        ds.set_node_values(device=d, region=SEMI, name="Electrons", init_from="n_i")
        tdeps = []
        if self.thermal:
            for r in regions:
                ds.node_solution(device=d, region=r, name="T")
                ds.edge_from_node_model(device=d, region=r, node_model="T")
                nr = len(ds.get_node_model_values(device=d, region=r, name="x"))
                ds.set_node_values(device=d, region=r, name="T", values=[P["T0"]] * nr)
            tdeps = ["T"]
        ani = "(unitx*unitx+rmu*unity*unity)" if two_d else "1"
        Tf = "pow((0.5*(T@n0+T@n1))/300,-aT)" if self.thermal else "pow(T0/300,-aT)"
        mup = f"(mu_p0*{ani}*{Tf})"
        mun = f"(mu_n0*{ani}*{Tf})"
        Vte = f"(kq*{self.Te()})"
        _, dp0, _, _ = self._node_stat("p", "@n0"); _, dp1, _, _ = self._node_stat("p", "@n1")
        _, dn0, _, _ = self._node_stat("n", "@n0"); _, dn1, _, _ = self._node_stat("n", "@n1")
        vdp = f"((Potential@n0-Potential@n1)/{Vte}+{dp0}-{dp1})"
        vdn = f"((Potential@n0-Potential@n1)/{Vte}-{dn0}+{dn1})"
        self.Jp_expr = f"(-q*{mup}*EdgeInverseLength*{Vte}*(Holes@n1*B({vdp})-Holes@n0*B(-{vdp})))"
        self.Jn_expr = f"(q*{mun}*EdgeInverseLength*{Vte}*(Electrons@n1*B(-{vdn})-Electrons@n0*B({vdn})))"
        _em(d, SEMI, "Jp", self.Jp_expr, ["Holes", "Potential"] + tdeps)
        _em(d, SEMI, "Jn", self.Jn_expr, ["Electrons", "Potential"] + tdeps)
        Np, dpn, _, Tn_ = self._node_stat("p", "")
        Nn, dnn, _, _ = self._node_stat("n", "")
        ni2 = f"({Nn}*{Np}*exp(-Eg/(kq*{Tn_}))*exp(-({dnn})-({dpn})))"
        U = f"((Electrons*Holes-{ni2})/(tau*(Electrons+ni)+tau*(Holes+ni)))"
        ndeps = ["Electrons", "Holes"] + (["T"] if self.noniso else [])
        _nm(d, SEMI, "U", U, ndeps)
        # tunnelling generation, semi-implicit: G = a - b*carrier (A/cm^3)
        _nm(d, SEMI, "Gen_p", f"q*{U}-(Gtp_a-Gtp_b*Holes)", ndeps)
        _nm(d, SEMI, "Gen_n", f"-q*{U}+(Gtn_a-Gtn_b*Electrons)", ndeps)
        ds.equation(device=d, region=SEMI, name="HCE", variable_name="Holes", edge_model="Jp",
                    node_model="Gen_p", variable_update="positive")
        ds.equation(device=d, region=SEMI, name="ECE", variable_name="Electrons", edge_model="Jn",
                    node_model="Gen_n", variable_update="positive")
        _nm(d, SEMI, "Rho", "-q*(Holes-Electrons+NetDoping)", ["Holes", "Electrons"])
        ds.equation(device=d, region=SEMI, name="PE", variable_name="Potential", node_model="Rho",
                    edge_model="DField", variable_update="log_damp")
        # ---- carrier contacts
        self.cexpr = {}
        for c in ("source", "drain"):
            if self.contact_bc[c] == "ohmic":
                _cnm(d, c, f"{c}_hbc", f"Holes-ni*exp(-(Potential-{c}_bias)/Vt300)", ["Holes", "Potential"])
                _cnm(d, c, f"{c}_ebc", f"Electrons-ni*exp((Potential-{c}_bias)/Vt300)", ["Electrons", "Potential"])
                ds.contact_equation(device=d, contact=c, name="HCE", node_model=f"{c}_hbc", edge_current_model="Jp")
                ds.contact_equation(device=d, contact=c, name="ECE", node_model=f"{c}_ebc", edge_current_model="Jn")
                continue
            T = self.Tn()
            ds.set_parameter(device=d, name=f"{c}_tv", value=self.T_vdW)
            if self.v_override:
                vp = f"({self.v_override})"; vn = vp
            else:
                vp = f"(A_p*{T}*{T}/(q*Nv300*pow({T}/300,1.5)))"
                vn = f"(A_n*{T}*{T}/(q*Nc300*pow({T}/300,1.5)))"
            if self.fd:
                peq = f"(Nv300*pow({T}/300,1.5)*{F12_expr(f'-phiBp/(kq*{T})')})"
                neq = f"(Nc300*pow({T}/300,1.5)*{F12_expr(f'-phiBn/(kq*{T})')})"
            else:
                peq = f"(Nv300*pow({T}/300,1.5)*exp(-phiBp/(kq*{T})))"
                neq = f"(Nc300*pow({T}/300,1.5)*exp(-phiBn/(kq*{T})))"
            Jps = f"(q*{vp}*{c}_tv*(1+Gp)*(Holes-{peq}))"           # hole current density leaving semi
            Jns = f"(-q*{vn}*{c}_tv*(1+Gn)*(Electrons-{neq}))"      # conventional current leaving semi
            self.cexpr[c] = dict(Jps=Jps, Jns=Jns, peq=peq, neq=neq)
            cdeps_p = ["Holes"] + (["T"] if self.noniso else [])
            cdeps_n = ["Electrons"] + (["T"] if self.noniso else [])
            _cnm(d, c, f"{c}_hbc", f"{Jps}*ContactSurfaceArea/NodeVolume", cdeps_p)
            _cnm(d, c, f"{c}_ebc", f"{Jns}*ContactSurfaceArea/NodeVolume", cdeps_n)
            ds.contact_equation(device=d, contact=c, name="HCE", node_model=f"{c}_hbc", edge_model="Jp",
                                edge_current_model="Jp")
            ds.contact_equation(device=d, contact=c, name="ECE", node_model=f"{c}_ebc", edge_model="Jn",
                                edge_current_model="Jn")
        self._paths()
        self._solve()            # carriers with T held at T0 (no heat equation yet)
        if self.thermal:
            self._build_thermal(two_d)
            self._solve()
        self.update_tunnel()

    def _build_thermal(self, two_d):
        d, P = self.name, self.P
        kap = "(kpar*unitx*unitx+kperp*unity*unity)" if two_d else "kpar"
        _em(d, SEMI, "Qflux", f"{kap}*(T@n0-T@n1)*EdgeInverseLength", ["T"])
        Vt0 = "(kq*T@n0)" if self.noniso else "Vt300"
        Vt1 = "(kq*T@n1)" if self.noniso else "Vt300"
        ep0, ep1 = self.eta("p", "@n0"), self.eta("p", "@n1")
        en0, en1 = self.eta("n", "@n0"), self.eta("n", "@n1")
        phip0 = f"(Potential@n0+{Vt0}*{ep0})"; phip1 = f"(Potential@n1+{Vt1}*{ep1})"
        phin0 = f"(Potential@n0-{Vt0}*{en0})"; phin1 = f"(Potential@n1-{Vt1}*{en1})"
        H = (f"({self.Jp_expr})*({phip0}-{phip1})*EdgeInverseLength"
             f"+({self.Jn_expr})*({phin0}-{phin1})*EdgeInverseLength")
        if P["peltier"]:
            Te = "(0.5*(T@n0+T@n1))"
            H += (f"+{Te}*kq*({self.Jp_expr})*({ep1}-{ep0})*EdgeInverseLength"
                  f"-{Te}*kq*({self.Jn_expr})*({en1}-{en0})*EdgeInverseLength")
        self.H_expr = H
        ds.edge_model(device=d, region=SEMI, name="Heat", equation=H)
        _em(d, SEMI, "NegHeat", f"-({H})", ["Potential", "Holes", "Electrons", "T"])
        _em(d, SEMI, "EdgeNodeVol", "0.5*EdgeCouple*EdgeLength")
        ds.set_parameter(name="edge_node0_volume_model", value="EdgeNodeVol")
        ds.set_parameter(name="edge_node1_volume_model", value="EdgeNodeVol")
        _nm(d, SEMI, "Cth", "cv*T", ["T"])
        _nm(d, SEMI, "NegHtun", "-Htun")
        ds.equation(device=d, region=SEMI, name="TE", variable_name="T", edge_model="Qflux",
                    edge_volume_model="NegHeat", node_model="NegHtun", time_node_model="Cth",
                    variable_update="default")
        # metal-contact thermal BC with interface Joule + Peltier junction heat
        for c in ("source", "drain"):
            Qs = "0"
            if c in self.cexpr:
                ce = self.cexpr[c]
                T = "T"
                es = self.eta("p", ""); en = self.eta("n", "")
                eqp = f"(-phiBp/(kq*{T}))"; eqn = f"(-phiBn/(kq*{T}))"
                if P["interface_heat"]:
                    Qs += (f"+{ce['Jps']}*kq*{T}*({es}-{eqp})"
                           f"+(-{ce['Jns']})*kq*{T}*({en}-{eqn})")
                if P["peltier"]:
                    Qs += (f"+{ce['Jps']}*{T}*kq*(2-{es})"
                           f"+{ce['Jns']}*{T}*(-kq)*(2-{en})")
            self.Qs_expr = getattr(self, "Qs_expr", {})
            self.Qs_expr[c] = Qs
            expr = f"(h_ct*(T-T0)-({Qs}))*ContactSurfaceArea/NodeVolume-Htun"
            _cnm(d, c, f"{c}_tbc", expr, ["T", "Holes", "Electrons"])
            ds.contact_equation(device=d, contact=c, name="TE", node_model=f"{c}_tbc", edge_model="Qflux",
                                edge_volume_model="NegHeat", edge_current_model="Qflux")
        if self.gate:
            ds.set_parameter(device=d, region=OX, name="kox", value=P["kappa_ox"])
            ds.set_parameter(device=d, region=OX, name="cvox", value=P["cv_ox"])
            _em(d, OX, "Qflux", "kox*(T@n0-T@n1)*EdgeInverseLength", ["T"])
            _nm(d, OX, "Cth", "cvox*T", ["T"])
            ds.equation(device=d, region=OX, name="TE", variable_name="T", edge_model="Qflux",
                        time_node_model="Cth", variable_update="default")
            _im(d, "semi_ox", "contT", "T@r0-T@r1", ["T"])
            ds.interface_equation(device=d, interface="semi_ox", name="TE", interface_model="contT",
                                  type="continuous")
            hg = 1.0 / P["R_TBR"] if P["R_TBR"] > 0 else 1e8
            ds.set_parameter(device=d, name="h_gate", value=hg)
            ds.contact_node_model(device=d, contact="gate", name="gate_tbc",
                                  equation="h_gate*(T-T0)*ContactSurfaceArea/NodeVolume")
            ds.contact_node_model(device=d, contact="gate", name="gate_tbc:T",
                                  equation="h_gate*ContactSurfaceArea/NodeVolume")
            ds.contact_equation(device=d, contact="gate", name="TE", node_model="gate_tbc",
                                edge_model="Qflux", edge_current_model="Qflux")
        else:
            ds.contact_node_model(device=d, contact="substrate", name="sub_tbc",
                                  equation="h_sub*(T-T0)*ContactSurfaceArea/NodeVolume-Htun")
            ds.contact_node_model(device=d, contact="substrate", name="sub_tbc:T",
                                  equation="h_sub*ContactSurfaceArea/NodeVolume")
            ds.contact_equation(device=d, contact="substrate", name="TE", node_model="sub_tbc",
                                edge_model="Qflux", edge_volume_model="NegHeat", edge_current_model="Qflux")

    # ------------------------------------------------------------------ geometry
    def nodev(self, name, region=SEMI):
        return np.array(ds.get_node_model_values(device=self.name, region=region, name=name))

    def _paths(self):
        """Contact nodes and straight inward paths (structured grid)."""
        x = self.nodev("x")
        y = self.nodev("y") if self.kind != "1d" else np.zeros_like(x)
        self.x, self.y = x, y
        tol = 1e-12
        self.cnodes, self.cpaths = {}, {}
        if self.kind == "top":
            top = np.where(np.abs(y) < tol)[0]
            sets = {"source": top[x[top] <= self.geo["Lc"] + tol],
                    "drain": top[x[top] >= self.geo["Lc"] + self.geo["Lch"] - tol]}
            for c, ids in sets.items():
                paths = []
                for i in ids:
                    col = np.where(np.abs(x - x[i]) < tol)[0]
                    col = col[np.argsort(y[col])]
                    paths.append(col)
                self.cnodes[c], self.cpaths[c] = ids, paths
        else:
            xs = np.unique(x)
            for c, xc, sgn in (("source", xs[0], 1), ("drain", xs[-1], -1)):
                ids = np.where(np.abs(x - xc) < tol)[0]
                paths = []
                for i in ids:
                    row = np.where(np.abs(y - y[i]) < tol)[0]
                    row = row[np.argsort(sgn * x[row])]
                    paths.append(row)
                self.cnodes[c], self.cpaths[c] = ids, paths

    # ------------------------------------------------------------------ tunnelling
    def _tsu_esaki(self, path, carrier, c, psi, eF):
        """Bins of the nonlocal Tsu-Esaki tunnelling integral along one contact path.
        Returns list of (node index, eps, transmission*weight [eV]) for energies below the
        barrier top at the contact; eps is the carrier energy (holes: -E)."""
        P = self.P
        dist = np.hypot(self.x[path] - self.x[path[0]], self.y[path] - self.y[path[0]])
        keep = dist <= P["tun_maxlen"]
        path, dist = path[keep], dist[keep]
        Ev = -psi[path] + self.Cv
        B = -Ev if carrier == "p" else Ev + self.Eg
        m = P["mh_tun"] if carrier == "p" else P["me_tun"]
        out = []
        if len(path) < 3:
            return path, out, B
        run = np.minimum.accumulate(B)
        kap0 = math.sqrt(2 * m * m0 * q) / hbar * 1e-2
        for i in range(1, len(path)):
            lo, hi = run[i], run[i - 1]
            if hi - lo <= 1e-7:
                continue
            for w, f in ((0.5, 0.211324865), (0.5, 0.788675134)):
                eps = lo + f * (hi - lo)
                over = np.clip(B[: i + 1] - eps, 0.0, None)
                integ = np.trapezoid(np.sqrt(over), dist[: i + 1])
                out.append((i, eps, math.exp(-2 * kap0 * integ) * w * (hi - lo)))
        return path, out, B

    def update_tunnel(self, iters=40, tol=1e-3, alpha=0.5):
        mode = self.P["tunneling"]
        if mode == "off" or self.kind == "1d" and mode == "local":
            return
        d = self.name
        for it in range(iters):
            psi = self.nodev("Potential")
            if mode == "local":
                Gp = np.zeros(self.nn); Gn = np.zeros(self.nn)
                for c in ("source", "drain"):
                    if self.contact_bc[c] != "schottky":
                        continue
                    for path in self.cpaths[c]:
                        i, j = path[0], path[1]
                        dist = np.hypot(self.x[i] - self.x[j], self.y[i] - self.y[j])
                        F = abs(psi[j] - psi[i]) / dist
                        Gp[i] = gamma_tfe(F, self.phiBp, self.P["mh_tun"])
                        Gn[i] = gamma_tfe(F, self.phiBn, self.P["me_tun"])
                old = self.nodev("Gp")
                ds.set_node_values(device=d, region=SEMI, name="Gp", values=list(Gp))
                ds.set_node_values(device=d, region=SEMI, name="Gn", values=list(Gn))
                self._solve()
                if np.max(np.abs(Gp - old) / (1 + np.abs(Gp))) < tol:
                    return
                continue
            vol = self.nodev("NodeVolume"); csa = self.nodev("ContactSurfaceArea")
            ph = self._phi()
            p, n = self.nodev("Holes"), self.nodev("Electrons")
            Vt = self.D["Vt"]
            A = {"p": self.D["A_p"], "n": self.D["A_n"]}
            # per-contact, per-carrier supply a and back-flow S (A/cm^3), and heat H (W/cm^3)
            fresh = {c: {k: np.zeros(self.nn) for k in ("pa", "pS", "na", "nS")} for c in ("source", "drain")}
            H = np.zeros(self.nn)
            for c in ("source", "drain"):
                if self.contact_bc[c] != "schottky":
                    continue
                Vc = ds.get_parameter(device=d, name=f"{c}_bias")
                for path in self.cpaths[c]:
                    acsa = csa[path[0]]
                    for car in ("p", "n"):
                        eFs = ph["phip"] if car == "p" else -ph["phin"]
                        eFm = Vc if car == "p" else -Vc
                        pp, bins, B = self._tsu_esaki(path, car, c, psi, eFs)
                        pref = A[car] * self.P["T0"] / kq * self.T_vdW
                        for i, eps, tw in bins:
                            node = pp[i]
                            sup = pref * tw * np.logaddexp(0, (eFm - eps) / Vt)
                            back = pref * tw * np.logaddexp(0, (eFs[node] - eps) / Vt)
                            fresh[c][car + "a"][node] += sup * acsa / vol[node]
                            fresh[c][car + "S"][node] += back * acsa / vol[node]
                            j = sup - back
                            # energy bookkeeping: metal side (eFm - eps); semiconductor side
                            # (eps - mean energy) with Peltier, (eps - local Fermi level) without
                            emean = B[i] + 2 * Vt if self.P["peltier"] else eFs[node]
                            H[node] += j * (eps - emean) * acsa / vol[node]
                            H[path[0]] += j * (eFm - eps) * acsa / vol[path[0]]
            first = not hasattr(self, "tun_src")
            a_ = 1.0 if first else alpha
            if first:
                self.tun_src = fresh
                self.tun_H = H
            else:
                for c in fresh:
                    for k in fresh[c]:
                        self.tun_src[c][k] += a_ * (fresh[c][k] - self.tun_src[c][k])
                self.tun_H += a_ * (H - self.tun_H)
            tot = {k: sum(self.tun_src[c][k] for c in self.tun_src) for k in ("pa", "pS", "na", "nS")}
            self.tun_lag = {"p": p.copy(), "n": n.copy()}
            self.tun_vol = vol
            ds.set_node_values(device=d, region=SEMI, name="Gtp_a", values=list(tot["pa"]))
            ds.set_node_values(device=d, region=SEMI, name="Gtn_a", values=list(tot["na"]))
            ds.set_node_values(device=d, region=SEMI, name="Gtp_b", values=list(np.where(p > 0, tot["pS"] / p, 0.0)))
            ds.set_node_values(device=d, region=SEMI, name="Gtn_b", values=list(np.where(n > 0, tot["nS"] / n, 0.0)))
            ds.set_node_values(device=d, region=SEMI, name="Htun", values=list(self.tun_H))
            self._solve()
            # convergence: carrier densities at nodes receiving tunnelling sources vs. lag
            mk = (tot["pa"] + tot["na"]) > 0
            if not mk.any():
                break
            dev = max(np.abs(self.nodev("Holes")[mk] / p[mk] - 1).max(),
                      np.abs(self.nodev("Electrons")[mk] / n[mk] - 1).max())
            self.tun_dev = dev
            if it > 0 and dev < tol:
                break

    def tunnel_current(self):
        """Terminal tunnelling current per contact (A/cm, + = conventional current into the
        device): G_c = a_c - (S_c/carrier_lag)*carrier_now, integrated over node volumes."""
        out = {"source": 0.0, "drain": 0.0}
        if self.P["tunneling"] != "nonlocal" or not hasattr(self, "tun_src"):
            return out
        now = {"p": self.nodev("Holes"), "n": self.nodev("Electrons")}
        for c, s_ in self.tun_src.items():
            for car, sgn in (("p", 1.0), ("n", -1.0)):
                ratio = np.where(self.tun_lag[car] > 0, now[car] / self.tun_lag[car], 0.0)
                G = s_[car + "a"] - s_[car + "S"] * ratio
                out[c] += sgn * float(np.sum(G * self.tun_vol))
        return out

    def tunnel_heat(self):
        return float(np.sum(self.nodev("Htun") * self.nodev("NodeVolume"))) if hasattr(self, "tun_H") else 0.0

    def _phi(self):
        """Absolute quasi-Fermi potentials (V) at nodes: phi_p = psi - Cv + Vt eta_p etc."""
        psi = self.nodev("Potential"); p = self.nodev("Holes"); n = self.nodev("Electrons")
        T = self.nodev("T") if self.noniso else np.full(self.nn, self.P["T0"])
        Vt = kq * T
        Nv = self.D["Nv"] * (T / 300) ** 1.5; Nc = self.D["Nc"] * (T / 300) ** 1.5
        ep = np.log(p / Nv) + (np.array([eta_inv(u) - math.log(u) for u in p / Nv]) if self.fd else 0)
        en = np.log(n / Nc) + (np.array([eta_inv(u) - math.log(u) for u in n / Nc]) if self.fd else 0)
        Cv = self.Cv
        return dict(phip=psi - Cv + Vt * ep, phin=psi - Cv - self.Eg - Vt * en, psi=psi)

    # ------------------------------------------------------------------ solving
    def _solve(self, abs_err=1e30, rel_err=1e-9, it=60):
        ds.solve(type="dc", absolute_error=abs_err, relative_error=rel_err, maximum_iterations=it)

    def _ramp(self, pname, attr, V, step, min_step=1e-4, max_step=0.1):
        d = self.name
        while abs(getattr(self, attr) - V) > 1e-12:
            s = math.copysign(min(step, abs(V - getattr(self, attr))), V - getattr(self, attr))
            target = getattr(self, attr) + s
            ds.set_parameter(device=d, name=pname, value=target)
            try:
                self._solve()
                self.update_tunnel(iters=3, tol=1e-2)
                setattr(self, attr, target)
                step = min(step * 1.5, max_step)
            except Exception:
                ds.set_parameter(device=d, name=pname, value=getattr(self, attr))
                step /= 2
                if step < min_step:
                    raise RuntimeError(f"ramp of {pname} failed at {getattr(self, attr)}")
                try:
                    self._solve()
                except Exception:
                    pass
        self.update_tunnel()

    def set_bias(self, V, step=0.05):
        self._ramp("drain_bias", "bias", V, step)

    def set_gate(self, VG, step=0.1):
        self._ramp("gate_bias", "vg", VG, step, max_step=0.25)

    def current(self):
        d = self.name
        Ip = ds.get_contact_current(device=d, contact="drain", equation="HCE")
        In = ds.get_contact_current(device=d, contact="drain", equation="ECE")
        It = self.tunnel_current()["drain"]
        return Ip + In + It, Ip, In, It

    def current_source(self):
        d = self.name
        It = self.tunnel_current()["source"]
        return (ds.get_contact_current(device=d, contact="source", equation="HCE")
                + ds.get_contact_current(device=d, contact="source", equation="ECE") + It)

    # ------------------------------------------------------------------ extraction
    def power_balance(self):
        d = self.name
        I = self.current()[0]
        out = dict(P_el=abs(I * self.bias))
        if not self.thermal:
            return out
        Heat = np.array(ds.get_edge_model_values(device=d, region=SEMI, name="Heat"))
        vol = np.array(ds.get_edge_model_values(device=d, region=SEMI, name="EdgeNodeVol"))
        out["P_bulk"] = float(np.sum(Heat * 2 * vol))
        T = self.nodev("T")
        csa = self.nodev("ContactSurfaceArea")
        Qs = self.junction_heat()
        out["P_junction"] = Qs["source"] + Qs["drain"]
        out["P_tunnel"] = self.tunnel_heat()
        P_ct = 0.0
        if self.P["h_contact"] > 0:
            for c in ("source", "drain"):
                ids = self.cnodes[c]
                P_ct += float(np.sum(self.P["h_contact"] * (T[ids] - self.P["T0"]) * csa[ids]))
        if self.gate:
            To = self.nodev("T", OX); yo = self.nodev("y", OX); co = self.nodev("ContactSurfaceArea", OX)
            g = np.where(np.abs(yo - (self.t + self.gate["t_ox"])) < 1e-12)[0]
            P_sink = float(np.sum((To[g] - self.P["T0"]) / self.P["R_TBR"] * co[g]))
        else:
            sub = np.where(np.abs(self.y - self.t) < 1e-12)[0]
            P_sink = float(np.sum(self.P["h_sub"] * (T[sub] - self.P["T0"]) * csa[sub]))
        out.update(P_sink=P_sink, P_contacts=P_ct, P_out=P_sink + P_ct, Tmax=float(T.max()),
                   T_source=float(T[self.cnodes["source"]].mean()), T_drain=float(T[self.cnodes["drain"]].mean()))
        return out

    def junction_heat(self):
        """Interface Joule + Peltier heat deposited at each Schottky contact (W/cm), evaluated
        with the same expressions as the thermal boundary condition."""
        P, Dd = self.P, self.D
        T = self.nodev("T") if self.thermal else np.full(self.nn, P["T0"])
        Tc = T if self.noniso else np.full(self.nn, P["T0"])
        p, n = self.nodev("Holes"), self.nodev("Electrons")
        csa = self.nodev("ContactSurfaceArea")
        Gp, Gn = self.nodev("Gp"), self.nodev("Gn")
        out = {}
        for c in ("source", "drain"):
            out[c] = 0.0
            if c not in self.cexpr:
                continue
            ids = self.cnodes[c]
            Tk, Tl = Tc[ids], T[ids]
            Nv = Dd["Nv"] * (Tk / 300) ** 1.5; Nc = Dd["Nc"] * (Tk / 300) ** 1.5
            if self.v_override:
                vp = vn = np.full(len(ids), self.v_override)
            else:
                vp = Dd["A_p"] * Tk ** 2 / (q * Nv); vn = Dd["A_n"] * Tk ** 2 / (q * Nc)
            if self.fd:
                peq = Nv * np.array([F12(-self.phiBp / (kq * x)) for x in Tk])
                neq = Nc * np.array([F12(-self.phiBn / (kq * x)) for x in Tk])
                ep = np.array([eta_inv(u) for u in p[ids] / Nv]); en = np.array([eta_inv(u) for u in n[ids] / Nc])
            else:
                peq = Nv * np.exp(-self.phiBp / (kq * Tk)); neq = Nc * np.exp(-self.phiBn / (kq * Tk))
                ep = np.log(p[ids] / Nv); en = np.log(n[ids] / Nc)
            Jps = q * vp * self.T_vdW * (1 + Gp[ids]) * (p[ids] - peq)
            Jns = -q * vn * self.T_vdW * (1 + Gn[ids]) * (n[ids] - neq)
            Q = np.zeros(len(ids))
            if P["interface_heat"]:
                Q += Jps * kq * Tl * (ep + self.phiBp / (kq * Tl)) + (-Jns) * kq * Tl * (en + self.phiBn / (kq * Tl))
            if P["peltier"]:
                Q += Jps * Tl * kq * (2 - ep) + Jns * Tl * (-kq) * (2 - en)
            out[c] = float(np.sum(Q * csa[ids]))
        return out

    def centerline(self):
        if self.kind == "1d":
            idx = np.argsort(self.x)
        else:
            ys = np.unique(self.y)
            yc = ys[np.argmin(np.abs(ys - 0.5 * self.t))]
            idx = np.where(np.abs(self.y - yc) < 1e-12)[0]
            idx = idx[np.argsort(self.x[idx])]
        ph = self._phi()
        out = dict(x=self.x[idx], psi=ph["psi"][idx], phip=ph["phip"][idx], phin=ph["phin"][idx],
                   p=self.nodev("Holes")[idx], n=self.nodev("Electrons")[idx])
        if self.thermal:
            out["T"] = self.nodev("T")[idx]
        return out

    def decompose(self):
        """Exact potential-drop decomposition with the majority-current quasi-Fermi potential."""
        I, Ip, In, It = self.current()
        V = self.bias
        cl = self.centerline()
        phi = cl["phip"] if abs(Ip) >= abs(In) else cl["phin"]
        x0, x1 = self.geo["x_ch0"], self.geo["x_ch1"]
        f0 = float(np.interp(x0, cl["x"], phi)); f1 = float(np.interp(x1, cl["x"], phi))
        Rcs, Rcd, Rch = f0 / I, (V - f1) / I, (f1 - f0) / I
        return dict(R_tot=V / I, R_ch=Rch, R_c=0.5 * (Rcs + Rcd), R_c_source=Rcs, R_c_drain=Rcd,
                    frac_contact=(Rcs + Rcd) / (V / I))

    def sheet_density(self, x_frac=0.5):
        """Carrier sheet densities (cm^-2) at a channel position (fraction of channel)."""
        xc = self.geo["x_ch0"] + x_frac * (self.geo["x_ch1"] - self.geo["x_ch0"])
        xs = np.unique(self.x)
        xn = xs[np.argmin(np.abs(xs - xc))]
        col = np.where(np.abs(self.x - xn) < 1e-12)[0]
        col = col[np.argsort(self.y[col])]
        p = self.nodev("Holes")[col]; n = self.nodev("Electrons")[col]; yy = self.y[col]
        return float(np.trapezoid(p, yy)), float(np.trapezoid(n, yy))

    # ------------------------------------------------------------------ transient
    def transient_step(self, V, t_end, nsteps=40, t_ramp=1e-11, nramp=40):
        """Ramp the drain from the present DC state to V over t_ramp (nramp sub-steps), then hold
        and integrate the lattice heat equation to t_end (carriers quasi-static, tunnelling
        sources frozen). Returns rows (t, Tmax, P_el, P_out)."""
        d = self.name
        ds.solve(type="transient_dc", absolute_error=1e30, relative_error=1e-9, maximum_iterations=60)
        V0 = self.bias
        res = []
        tr = np.linspace(t_ramp / nramp, t_ramp, nramp)
        t_prev = 0.0
        for tk in tr:
            ds.set_parameter(device=d, name="drain_bias", value=V0 + (V - V0) * tk / t_ramp)
            ds.solve(type="transient_bdf1", absolute_error=1e30, relative_error=1e-8,
                     maximum_iterations=80, tdelta=tk - t_prev, charge_error=1e30)
            t_prev = tk
        self.bias = V
        pb = self.power_balance()
        res.append((t_prev, pb["Tmax"], pb["P_el"], pb["P_out"]))
        for tk in np.geomspace(t_ramp * 1.5, t_end, nsteps):
            ds.solve(type="transient_bdf1", absolute_error=1e30, relative_error=1e-8,
                     maximum_iterations=60, tdelta=tk - t_prev, charge_error=1e30)
            t_prev = tk
            pb = self.power_balance()
            res.append((tk, pb["Tmax"], pb["P_el"], pb["P_out"]))
        return np.array(res)
