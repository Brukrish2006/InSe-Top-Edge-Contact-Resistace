"""
inse2d.py -- corrected (v2) DEVSIM model for metal / few-layer gamma-InSe contacts.

This module replaces the v1 pipeline in ../src. Changes relative to v1, each tied to a
reviewer comment (R1 = Reviewer 1, R3 = Reviewer 3, R2 = Reviewer 2):

 * Poisson equation is solved with the *non-equilibrium* carrier densities
   rho = q (p - n + N_D - N_A)  (v1 kept the equilibrium Boltzmann charge
   "IntrinsicCharge" after drift-diffusion was switched on, so injected carriers
   never fed back into the electrostatics).                                   [R1.2]
 * Background doping is p-type: NetDoping = N_D - N_A = -N_A (v1 used +1e16,
   i.e. n-type, contradicting the text).                                       [R1.2]
 * Metal/InSe carrier boundary condition is a genuine flux (Robin) condition,
       J_p,s = q v_p,eff (p_s - p_eq),
   balanced against the drift-diffusion current arriving from the semiconductor.
   In v1 the expression q v (1+Gamma)(p - p_eq) was supplied as the *residual*
   of the DEVSIM contact equation, which enforces p_s = p_eq for any nonzero
   prefactor (Dirichlet), so Gamma could never influence the solution.        [R1.1]
 * v_p,eff = v_R * T_vdW * (1 + Gamma(F)):  v_R = A* T^2 /(q N_V) (Richardson
   velocity); T_vdW is a WKB transmission through the van der Waals gap (top
   contacts only); Gamma(F) is the field-dependent thermionic-field-emission
   enhancement computed from the local interface field by WKB through a
   triangular barrier and updated self-consistently.                          [R1.1]
 * Contact potential uses the correct intrinsic-level offset
   (+0.5 V_t ln(N_C/N_V); v1 had the opposite sign, so the potential boundary
   condition and p_eq were inconsistent by ~34 meV).
 * Optional Fermi-level pinning: Phi_eff = Phi_CNL + S (Phi_M - Phi_CNL).      [R3.2]
 * Joule heating uses H = -J_p . grad(phi_p) - J_n . grad(phi_n) (quasi-Fermi
   potentials, non-negative), and the interface (contact) dissipation is
   tracked separately so that the energy balance I*V = P_bulk + P_interface
   can be checked.                                                            [R3.6]
 * Thermal boundary at the substrate is a Robin condition with an explicit
   effective heat-transfer coefficient h_sub (thermal boundary resistance +
   oxide); metal contacts are adiabatic by default or heat-sinking (option).   [R3.6]
 * Resistance, temperature and power are always extracted at the *same* bias,
   and the contact / channel voltage drops are separated from the hole
   quasi-Fermi potential along the channel.                              [R1.3, R3.1, R3.5]

Units: DEVSIM cgs-like (cm, V, A, s). 2D simulations are per unit width, so
currents are A/cm and resistances Ohm*cm.
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
h_planck = 6.62607015e-34
hbar = h_planck / (2 * math.pi)
m0 = 9.1093837015e-31
eps0 = 8.8541878128e-14  # F/cm

REGION = "InSe"

# ---------------------------------------------------------------------------
# Parameters (see params_table() for provenance categories, R3.7)
# ---------------------------------------------------------------------------
DEFAULTS = dict(
    T0=300.0,
    Eg=1.30,          # eV
    chi=3.70,         # eV  -> E_V = chi + Eg = 5.00 eV below vacuum
    me_dos=0.07,
    mh_dos=0.17,
    mh_tun=0.17,      # hole tunnelling mass (assumed = DOS mass)
    eps_r=6.4,
    mu_n=1000.0,      # cm^2/Vs
    mu_p=50.0,        # cm^2/Vs
    mu_T_exp=1.5,     # mu ~ (T/300)^-mu_T_exp
    tau=1e-9,         # s (SRH)
    NA=1e16,          # cm^-3 (p-type background)
    ND=0.0,
    kappa=0.085,      # W/(cm K)
    h_sub=1.35e3,     # W/(cm^2 K): 1/(R_TBR + t_ox/k_ox), R_TBR=1e-8 m^2K/W, 90 nm SiO2
    h_contact=0.0,    # W/(cm^2 K) heat sinking through metal contacts (0 = adiabatic)
    S_pin=1.0,        # pinning factor (1 = Schottky-Mott)
    Phi_CNL=4.35,     # eV, charge-neutrality level (assumed mid-gap) used only if S_pin<1
    tunneling=True,   # field-dependent TFE enhancement Gamma(F)
    T_vdW_top=None,   # None -> computed from WKB with d_vdW, Phi_vdW below
    d_vdW=0.30e-7,    # cm, vdW gap width (top contacts)
    Phi_vdW=2.0,      # eV, effective tunnel-barrier height of the vdW gap
    thermal=True,
)


def derived(P):
    T = P["T0"]
    Vt = kB * T / q
    fac = 2 * (2 * math.pi * kB * T / h_planck ** 2) ** 1.5 * 1e-6
    Nc = fac * (P["me_dos"] * m0) ** 1.5
    Nv = fac * (P["mh_dos"] * m0) ** 1.5
    ni = math.sqrt(Nc * Nv) * math.exp(-P["Eg"] / (2 * Vt))
    A_p = 120.17 * P["mh_dos"]
    A_n = 120.17 * P["me_dos"]
    vR_p = A_p * T ** 2 / (q * Nv)
    vR_n = A_n * T ** 2 / (q * Nc)
    return dict(Vt=Vt, Nc=Nc, Nv=Nv, ni=ni, vR_p=vR_p, vR_n=vR_n, A_p=A_p, A_n=A_n)


def phi_eff(P, Phi_M):
    """Effective metal work function seen by the semiconductor (pinning)."""
    return P["Phi_CNL"] + P["S_pin"] * (Phi_M - P["Phi_CNL"])


def barrier_heights(P, Phi_M):
    """Schottky barrier heights (eV) for holes and electrons. Negative hole barriers
    (Phi_eff above E_V) are clipped to zero because Boltzmann drift-diffusion cannot
    describe a degenerate interface; the clipping is reported explicitly."""
    Pe = phi_eff(P, Phi_M)
    Ev = P["chi"] + P["Eg"]
    phiBp_raw = Ev - Pe
    phiBp = max(phiBp_raw, 0.0)
    Pe_used = Ev - phiBp
    phiBn = Pe_used - P["chi"]
    return dict(Phi_eff=Pe, Phi_used=Pe_used, phiBp=phiBp, phiBp_raw=phiBp_raw, phiBn=phiBn,
                clipped=phiBp_raw < 0)


def wkb_vdw(P):
    if P["T_vdW_top"] is not None:
        return P["T_vdW_top"]
    if P["d_vdW"] <= 0:
        return 1.0
    kappa = math.sqrt(2 * m0 * P["Phi_vdW"] * q) / hbar  # 1/m
    return math.exp(-2 * kappa * P["d_vdW"] * 1e-2)


def gamma_tfe(F_Vcm, phiB_eV, m_tun, T=300.0):
    """Thermionic-field-emission enhancement factor (dimensionless).

    Gamma = (1/kT) * int_0^phiB exp(E/kT) * exp(-(4/3) sqrt(2 m*) E^{3/2} / (q hbar F)) dE,
    E measured downward from the barrier top, WKB through a triangular barrier with
    the local interface field F. Gamma -> 0 for F -> 0 (pure thermionic emission).
    """
    if phiB_eV <= 0 or F_Vcm <= 0:
        return 0.0
    kT = kB * T
    F = F_Vcm * 100.0  # V/m
    b = (4.0 / 3.0) * math.sqrt(2 * m_tun * m0) / (q * hbar * F)
    E = np.linspace(0.0, phiB_eV * q, 4001)
    expo = E / kT - b * E ** 1.5
    g = np.trapezoid(np.exp(expo), E) / kT
    return float(g)


# ---------------------------------------------------------------------------
# Meshes
# ---------------------------------------------------------------------------
def mesh_top(name, Lc, Lch, t, mf=1.0):
    """Top contacts of length Lc on the upper InSe surface, separated by Lch."""
    ds.create_2d_mesh(mesh=name)
    xs = [(0.0, 2e-7), (Lc, 0.25e-7), (Lc + 0.5 * Lch, 2e-7), (Lc + Lch, 0.25e-7),
          (2 * Lc + Lch, 2e-7)]
    for x, ps in xs:
        ds.add_2d_mesh_line(mesh=name, dir="x", pos=x, ps=ps / mf)
    for y, ps in [(-1e-7, 0.5e-7), (0.0, 0.1e-7), (0.5 * t, 0.5e-7), (t, 0.25e-7), (t + 1e-7, 0.5e-7)]:
        ds.add_2d_mesh_line(mesh=name, dir="y", pos=y, ps=ps / mf)
    L = 2 * Lc + Lch
    ds.add_2d_region(mesh=name, material="InSe", region=REGION, xl=0, xh=L, yl=0, yh=t)
    # inert neighbour regions (DEVSIM needs a region on the far side of a 2D contact)
    ds.add_2d_region(mesh=name, material="metal", region="m_src", xl=0, xh=Lc, yl=-1e-7, yh=0)
    ds.add_2d_region(mesh=name, material="metal", region="m_drn", xl=Lc + Lch, xh=L, yl=-1e-7, yh=0)
    ds.add_2d_region(mesh=name, material="metal", region="m_sub", xl=0, xh=L, yl=t, yh=t + 1e-7)
    b = 1e-10
    ds.add_2d_contact(mesh=name, name="source", material="metal", region=REGION,
                      xl=0, xh=Lc, yl=0, yh=0, bloat=b)
    ds.add_2d_contact(mesh=name, name="drain", material="metal", region=REGION,
                      xl=Lc + Lch, xh=L, yl=0, yh=0, bloat=b)
    ds.add_2d_contact(mesh=name, name="substrate", material="metal", region=REGION,
                      xl=0, xh=L, yl=t, yh=t, bloat=b)
    ds.finalize_mesh(mesh=name)
    ds.create_device(mesh=name, device=name)
    return dict(x_ch0=Lc, x_ch1=Lc + Lch, normal="y", Lc=Lc, Lch=Lch, t=t, geom="top")


def mesh_edge(name, Lch, t, mf=1.0):
    """Edge contacts on the two vertical faces of a channel of length Lch."""
    ds.create_2d_mesh(mesh=name)
    for x, ps in [(-1e-7, 0.5e-7), (0.0, 0.25e-7), (0.5 * Lch, 2e-7), (Lch, 0.25e-7), (Lch + 1e-7, 0.5e-7)]:
        ds.add_2d_mesh_line(mesh=name, dir="x", pos=x, ps=ps / mf)
    for y, ps in [(0.0, 0.25e-7), (0.5 * t, 0.5e-7), (t, 0.25e-7), (t + 1e-7, 0.5e-7)]:
        ds.add_2d_mesh_line(mesh=name, dir="y", pos=y, ps=ps / mf)
    ds.add_2d_region(mesh=name, material="InSe", region=REGION, xl=0, xh=Lch, yl=0, yh=t)
    ds.add_2d_region(mesh=name, material="metal", region="m_src", xl=-1e-7, xh=0, yl=0, yh=t)
    ds.add_2d_region(mesh=name, material="metal", region="m_drn", xl=Lch, xh=Lch + 1e-7, yl=0, yh=t)
    ds.add_2d_region(mesh=name, material="metal", region="m_sub", xl=0, xh=Lch, yl=t, yh=t + 1e-7)
    b = 1e-10
    ds.add_2d_contact(mesh=name, name="source", material="metal", region=REGION,
                      xl=0, xh=0, yl=0, yh=t, bloat=b)
    ds.add_2d_contact(mesh=name, name="drain", material="metal", region=REGION,
                      xl=Lch, xh=Lch, yl=0, yh=t, bloat=b)
    ds.add_2d_contact(mesh=name, name="substrate", material="metal", region=REGION,
                      xl=0, xh=Lch, yl=t, yh=t, bloat=b)
    ds.finalize_mesh(mesh=name)
    ds.create_device(mesh=name, device=name)
    return dict(x_ch0=0.0, x_ch1=Lch, normal="x", Lc=0.0, Lch=Lch, t=t, geom="edge")


def mesh_1d(name, L, ps=1e-7, ps_c=None):
    ps_c = ps_c or ps
    ds.create_1d_mesh(mesh=name)
    ds.add_1d_mesh_line(mesh=name, pos=0.0, ps=ps_c, tag="l")
    ds.add_1d_mesh_line(mesh=name, pos=0.5 * L, ps=ps)
    ds.add_1d_mesh_line(mesh=name, pos=L, ps=ps_c, tag="h")
    ds.add_1d_contact(mesh=name, name="source", tag="l", material="metal")
    ds.add_1d_contact(mesh=name, name="drain", tag="h", material="metal")
    ds.add_1d_region(mesh=name, material="InSe", region=REGION, tag1="l", tag2="h")
    ds.finalize_mesh(mesh=name)
    ds.create_device(mesh=name, device=name)
    return dict(x_ch0=0.0, x_ch1=L, normal="x", Lc=0.0, Lch=L, t=None, geom="1d")


# ---------------------------------------------------------------------------
# Model helpers
# ---------------------------------------------------------------------------
def _nm(d, name, expr, deps=()):
    ds.node_model(device=d, region=REGION, name=name, equation=expr)
    for v in deps:
        ds.node_model(device=d, region=REGION, name=f"{name}:{v}", equation=f"diff({expr},{v})")


def _em(d, name, expr, deps=()):
    ds.edge_model(device=d, region=REGION, name=name, equation=expr)
    for v in deps:
        for k in ("n0", "n1"):
            ds.edge_model(device=d, region=REGION, name=f"{name}:{v}@{k}",
                          equation=f"diff({expr},{v}@{k})")


def _cnm(d, c, name, expr, deps=()):
    ds.contact_node_model(device=d, contact=c, name=name, equation=expr)
    for v in deps:
        ds.contact_node_model(device=d, contact=c, name=f"{name}:{v}", equation=f"diff({expr},{v})")


class Device:
    """One simulated device. Build with Device(kind, P, **geometry)."""

    def __init__(self, kind, P=None, Phi_M=5.0, name="dev", mf=1.0, Lc=50e-7, Lch=50e-7,
                 t=5e-7, L1d=1e-4, contact_bc=("schottky", "schottky"), ps1d=1e-7,
                 v_override=None):
        ds.reset_devsim()
        self.P = dict(DEFAULTS)
        if P:
            self.P.update(P)
        self.D = derived(self.P)
        self.kind = kind
        self.name = name
        self.Phi_M = Phi_M
        self.bh = barrier_heights(self.P, Phi_M)
        if kind == "top":
            self.geo = mesh_top(name, Lc, Lch, t, mf)
        elif kind == "edge":
            self.geo = mesh_edge(name, Lch, t, mf)
        elif kind == "1d":
            self.geo = mesh_1d(name, L1d, ps=ps1d)
        else:
            raise ValueError(kind)
        self.contact_bc = dict(zip(("source", "drain"), contact_bc))
        self.v_override = v_override
        self.T_vdW = wkb_vdw(self.P) if kind == "top" else 1.0
        self.bias = 0.0
        self._build()

    # ---------------- build ----------------
    def _build(self):
        d, P, D = self.name, self.P, self.D
        thermal = P["thermal"] and self.kind != "1d"
        self.thermal = thermal
        par = dict(q=q, Vt=D["Vt"], eps=P["eps_r"] * eps0, ni=D["ni"], Nc=D["Nc"], Nv=D["Nv"],
                   NetDoping=P["ND"] - P["NA"], mu_n0=P["mu_n"], mu_p0=P["mu_p"],
                   aT=P["mu_T_exp"], taun=P["tau"], taup=P["tau"], kappa=P["kappa"],
                   T0=P["T0"], h_sub=P["h_sub"], h_ct=P["h_contact"])
        for k, v in par.items():
            ds.set_parameter(device=d, name=k, value=v)
        for s in ("Potential", "Electrons", "Holes", "Gp", "Gn"):
            ds.node_solution(device=d, region=REGION, name=s)
            ds.edge_from_node_model(device=d, region=REGION, node_model=s)
        nn = len(ds.get_node_model_values(device=d, region=REGION, name="x"))
        self.nn = nn
        # ---- contact built-in potential (consistent with p_eq / n_eq)
        self.psi_bi = (P["chi"] - self.bh["Phi_used"] + P["Eg"] / 2
                       + 0.5 * D["Vt"] * math.log(D["Nc"] / D["Nv"]))
        NA = P["NA"] - P["ND"]
        psi_bulk = -D["Vt"] * math.log(NA / D["ni"]) if NA > 0 else D["Vt"] * math.log(-NA / D["ni"])
        ds.set_node_values(device=d, region=REGION, name="Potential", values=[psi_bulk] * nn)
        for s in ("Gp", "Gn"):
            ds.set_node_values(device=d, region=REGION, name=s, values=[0.0] * nn)
        # ---- Poisson, equilibrium (potential-only) first
        _em(d, "EField", "(Potential@n0-Potential@n1)*EdgeInverseLength", ["Potential"])
        _em(d, "DField", "eps*EField", ["Potential"])
        _nm(d, "p_i", "ni*exp(-Potential/Vt)", ["Potential"])
        _nm(d, "n_i", "ni*exp(Potential/Vt)", ["Potential"])
        _nm(d, "RhoEq", "-q*(p_i-n_i+NetDoping)", ["Potential"])
        ds.equation(device=d, region=REGION, name="PE", variable_name="Potential",
                    node_model="RhoEq", edge_model="DField", variable_update="log_damp")
        for c in ("source", "drain"):
            ds.set_parameter(device=d, name=f"{c}_bias", value=0.0)
            if self.contact_bc[c] == "ohmic":
                psi_c = psi_bulk
            else:
                psi_c = self.psi_bi
            ds.set_parameter(device=d, name=f"{c}_psi", value=psi_c)
            _cnm(d, c, f"{c}_pbc", f"Potential-{c}_bias-{c}_psi", ["Potential"])
            ds.contact_equation(device=d, contact=c, name="PE", node_model=f"{c}_pbc",
                                edge_charge_model="DField")
        self._solve(abs_err=1e30, rel_err=1e-10, it=100)
        # ---- drift-diffusion
        ds.set_node_values(device=d, region=REGION, name="Holes", init_from="p_i")
        ds.set_node_values(device=d, region=REGION, name="Electrons", init_from="n_i")
        if thermal:
            ds.node_solution(device=d, region=REGION, name="T")
            ds.edge_from_node_model(device=d, region=REGION, node_model="T")
            ds.set_node_values(device=d, region=REGION, name="T", values=[P["T0"]] * nn)
            mup = "mu_p0*pow(0.5*(T@n0+T@n1)/300.0,-aT)"
            mun = "mu_n0*pow(0.5*(T@n0+T@n1)/300.0,-aT)"
            tdeps = ["T"]
        else:
            mup, mun, tdeps = "mu_p0", "mu_n0", []
        _em(d, "vd", "(Potential@n0-Potential@n1)/Vt", ["Potential"])
        Jp = f"-q*({mup})*EdgeInverseLength*Vt*(Holes@n1*B(vd)-Holes@n0*B(-vd))"
        Jn = f"q*({mun})*EdgeInverseLength*Vt*(Electrons@n1*B(-vd)-Electrons@n0*B(vd))"
        sub = lambda e: e.replace("B(vd)", "B((Potential@n0-Potential@n1)/Vt)").replace(
            "B(-vd)", "B(-(Potential@n0-Potential@n1)/Vt)")
        self.Jp_expr, self.Jn_expr = sub(Jp), sub(Jn)
        _em(d, "Jp", self.Jp_expr, ["Holes", "Potential"] + tdeps)
        _em(d, "Jn", self.Jn_expr, ["Electrons", "Potential"] + tdeps)
        U = "(Electrons*Holes-ni^2)/(taup*(Electrons+ni)+taun*(Holes+ni))"
        _nm(d, "U", U, ["Electrons", "Holes"])
        _nm(d, "Gen_p", f"q*{U}", ["Electrons", "Holes"])
        _nm(d, "Gen_n", f"-q*{U}", ["Electrons", "Holes"])
        ds.equation(device=d, region=REGION, name="HCE", variable_name="Holes",
                    edge_model="Jp", node_model="Gen_p", variable_update="positive")
        ds.equation(device=d, region=REGION, name="ECE", variable_name="Electrons",
                    edge_model="Jn", node_model="Gen_n", variable_update="positive")
        # full Poisson with non-equilibrium carriers (fix of R1.2)
        _nm(d, "Rho", "-q*(Holes-Electrons+NetDoping)", ["Holes", "Electrons"])
        ds.equation(device=d, region=REGION, name="PE", variable_name="Potential",
                    node_model="Rho", edge_model="DField", variable_update="log_damp")
        # ---- carrier contacts
        for c in ("source", "drain"):
            if self.contact_bc[c] == "ohmic":
                pe = f"ni*exp(-(Potential-{c}_bias)/Vt)"
                ne = f"ni*exp((Potential-{c}_bias)/Vt)"
                _cnm(d, c, f"{c}_hbc", f"Holes-{pe}", ["Holes", "Potential"])
                _cnm(d, c, f"{c}_ebc", f"Electrons-{ne}", ["Electrons", "Potential"])
                ds.contact_equation(device=d, contact=c, name="HCE", node_model=f"{c}_hbc",
                                    edge_current_model="Jp")
                ds.contact_equation(device=d, contact=c, name="ECE", node_model=f"{c}_ebc",
                                    edge_current_model="Jn")
                continue
            vp = self.v_override if self.v_override else self.D["vR_p"]
            vn = self.v_override if self.v_override else self.D["vR_n"]
            ds.set_parameter(device=d, name=f"{c}_vp", value=vp * self.T_vdW)
            ds.set_parameter(device=d, name=f"{c}_vn", value=vn * self.T_vdW)
            ds.set_parameter(device=d, name="p_eq", value=D["ni"] * math.exp(-self.psi_bi / D["Vt"]))
            ds.set_parameter(device=d, name="n_eq", value=D["ni"] * math.exp(self.psi_bi / D["Vt"]))
            fp = f"q*{c}_vp*(1+Gp)*(Holes-p_eq)*ContactSurfaceArea/NodeVolume"
            fn = f"-q*{c}_vn*(1+Gn)*(Electrons-n_eq)*ContactSurfaceArea/NodeVolume"
            _cnm(d, c, f"{c}_hbc", fp, ["Holes"])
            _cnm(d, c, f"{c}_ebc", fn, ["Electrons"])
            ds.contact_equation(device=d, contact=c, name="HCE", node_model=f"{c}_hbc",
                                edge_model="Jp", edge_current_model="Jp")
            ds.contact_equation(device=d, contact=c, name="ECE", node_model=f"{c}_ebc",
                                edge_model="Jn", edge_current_model="Jn")
        # ---- thermal
        if thermal:
            _em(d, "Qflux", "kappa*(T@n0-T@n1)*EdgeInverseLength", ["T"])
            H = (f"({self.Jp_expr})*((Potential@n0+Vt*log(Holes@n0/ni))-(Potential@n1+Vt*log(Holes@n1/ni)))*EdgeInverseLength"
                 f"+({self.Jn_expr})*((Potential@n0-Vt*log(Electrons@n0/ni))-(Potential@n1-Vt*log(Electrons@n1/ni)))*EdgeInverseLength")
            ds.edge_model(device=d, region=REGION, name="Heat", equation=H)
            _em(d, "NegHeat", f"-({H})", ["Potential", "Holes", "Electrons", "T"])
            _em(d, "EdgeNodeVol", "0.5*EdgeCouple*EdgeLength")
            ds.set_parameter(name="edge_node0_volume_model", value="EdgeNodeVol")
            ds.set_parameter(name="edge_node1_volume_model", value="EdgeNodeVol")
            ds.equation(device=d, region=REGION, name="TE", variable_name="T",
                        edge_model="Qflux", edge_volume_model="NegHeat", variable_update="default")
            _cnm(d, "substrate", "sub_tbc", "h_sub*(T-T0)*ContactSurfaceArea/NodeVolume", ["T"])
            ds.contact_equation(device=d, contact="substrate", name="TE", node_model="sub_tbc",
                                edge_model="Qflux", edge_volume_model="NegHeat",
                                edge_current_model="Qflux")
            if self.P["h_contact"] > 0:
                for c in ("source", "drain"):
                    _cnm(d, c, f"{c}_tbc", "h_ct*(T-T0)*ContactSurfaceArea/NodeVolume", ["T"])
                    ds.contact_equation(device=d, contact=c, name="TE", node_model=f"{c}_tbc",
                                        edge_model="Qflux", edge_volume_model="NegHeat",
                                        edge_current_model="Qflux")
        self._contact_nodes()
        self._solve()
        if self.P["tunneling"]:
            self.update_gamma(iters=3)

    # ---------------- geometry helpers ----------------
    def _contact_nodes(self):
        d = self.name
        x = np.array(ds.get_node_model_values(device=d, region=REGION, name="x"))
        y = np.array(ds.get_node_model_values(device=d, region=REGION, name="y")) if self.kind != "1d" \
            else np.zeros_like(x)
        self.x, self.y = x, y
        g = self.geo
        tol = 1e-12
        self.cnodes = {}
        self.cnbr = {}
        if self.kind == "top":
            top = np.where(np.abs(y) < tol)[0]
            src = top[x[top] <= g["Lc"] + tol]
            drn = top[x[top] >= g["Lc"] + g["Lch"] - tol]
            ys = np.unique(y)
            y1 = ys[1]
            for c, ids in (("source", src), ("drain", drn)):
                nb = []
                for i in ids:
                    j = np.where((np.abs(x - x[i]) < tol) & (np.abs(y - y1) < tol))[0][0]
                    nb.append(j)
                self.cnodes[c], self.cnbr[c] = ids, np.array(nb)
        else:
            xs = np.unique(x)
            for c, xc, xn in (("source", xs[0], xs[1]), ("drain", xs[-1], xs[-2])):
                ids = np.where(np.abs(x - xc) < tol)[0]
                nb = []
                for i in ids:
                    j = np.where((np.abs(x - xn) < tol) & (np.abs(y - y[i]) < tol))[0][0]
                    nb.append(j)
                self.cnodes[c], self.cnbr[c] = ids, np.array(nb)

    def nodev(self, name):
        return np.array(ds.get_node_model_values(device=self.name, region=REGION, name=name))

    # ---------------- tunneling update ----------------
    def update_gamma(self, iters=3, tol=1e-3):
        if self.contact_bc["source"] == "ohmic" and self.contact_bc["drain"] == "ohmic":
            return
        d = self.name
        Tlat = self.P["T0"]
        hist = []
        for _ in range(iters):
            psi = self.nodev("Potential")
            Gp = np.zeros(self.nn)
            Gn = np.zeros(self.nn)
            for c in ("source", "drain"):
                if self.contact_bc[c] != "schottky":
                    continue
                ids, nb = self.cnodes[c], self.cnbr[c]
                dist = np.hypot(self.x[ids] - self.x[nb], self.y[ids] - self.y[nb])
                F = (psi[nb] - psi[ids]) / dist  # V/cm, >0 when field points from metal into semic.
                for k, i in enumerate(ids):
                    # holes tunnel when the band bending pushes holes away from the interface
                    # (E points into semiconductor at a hole-depleting contact) -> use |F|
                    Gp[i] = gamma_tfe(abs(F[k]), self.bh["phiBp"], self.P["mh_tun"], Tlat)
                    Gn[i] = gamma_tfe(abs(F[k]), self.bh["phiBn"], self.P["me_dos"], Tlat)
            old = self.nodev("Gp")
            ds.set_node_values(device=d, region=REGION, name="Gp", values=list(Gp))
            ds.set_node_values(device=d, region=REGION, name="Gn", values=list(Gn))
            self._solve()
            ch = np.max(np.abs(Gp - old) / (1 + np.abs(Gp)))
            hist.append(ch)
            if ch < tol:
                break
        self.gamma_hist = hist
        return hist

    # ---------------- solving ----------------
    def _solve(self, abs_err=1e30, rel_err=1e-9, it=60):
        ds.solve(type="dc", absolute_error=abs_err, relative_error=rel_err, maximum_iterations=it)

    def set_bias(self, V, step=0.05, min_step=1e-4):
        """Ramp the drain bias from the current value to V with adaptive steps."""
        d = self.name
        while abs(self.bias - V) > 1e-12:
            s = math.copysign(min(step, abs(V - self.bias)), V - self.bias)
            target = self.bias + s
            ds.set_parameter(device=d, name="drain_bias", value=target)
            try:
                self._solve()
                if self.P["tunneling"]:
                    self.update_gamma(iters=2)
                self.bias = target
                step = min(step * 1.5, 0.1)
            except Exception:
                ds.set_parameter(device=d, name="drain_bias", value=self.bias)
                step /= 2
                if step < min_step:
                    raise RuntimeError(f"bias ramp failed at {self.bias} V")
                try:
                    self._solve()
                except Exception:
                    pass

    def current(self):
        d = self.name
        Ip = ds.get_contact_current(device=d, contact="drain", equation="HCE")
        In = ds.get_contact_current(device=d, contact="drain", equation="ECE")
        return Ip + In, Ip, In

    def current_source(self):
        d = self.name
        return (ds.get_contact_current(device=d, contact="source", equation="HCE")
                + ds.get_contact_current(device=d, contact="source", equation="ECE"))

    # ---------------- extraction ----------------
    def power_balance(self):
        """Return dict with I*V, bulk Joule power, interface power and heat leaving through
        the thermal boundaries (all W/cm in 2D)."""
        d = self.name
        I, _, _ = self.current()
        V = self.bias
        P_el = abs(I * V)
        out = dict(P_el=P_el)
        out.update(self._interface_power())
        if not self.thermal:
            return out
        Heat = np.array(ds.get_edge_model_values(device=d, region=REGION, name="Heat"))
        vol = np.array(ds.get_edge_model_values(device=d, region=REGION, name="EdgeNodeVol"))
        P_bulk = float(np.sum(Heat * 2 * vol))
        T = self.nodev("T")
        csa = self.nodev("ContactSurfaceArea")
        sub = np.where(np.abs(self.y - self.geo["t"]) < 1e-12)[0]
        P_sub = float(np.sum(self.P["h_sub"] * (T[sub] - self.P["T0"]) * csa[sub]))
        P_ct = 0.0
        if self.P["h_contact"] > 0:
            for c in ("source", "drain"):
                ids = self.cnodes[c]
                P_ct += float(np.sum(self.P["h_contact"] * (T[ids] - self.P["T0"]) * csa[ids]))
        out.update(P_bulk=P_bulk, P_sub=P_sub, P_contacts=P_ct,
                   P_out=P_sub + P_ct, Tmax=float(T.max()))
        return out

    def _interface_power(self):
        """J_s * (quasi-Fermi jump) summed over both Schottky contacts (W/cm)."""
        if not hasattr(self, "cnodes") or self.kind == "1d":
            return {}
        p = self.nodev("Holes"); n = self.nodev("Electrons")
        csa = self.nodev("ContactSurfaceArea")
        Vt = self.D["Vt"]
        peq = self.D["ni"] * math.exp(-self.psi_bi / Vt)
        neq = self.D["ni"] * math.exp(self.psi_bi / Vt)
        Gp = self.nodev("Gp"); Gn = self.nodev("Gn")
        P_if = 0.0
        for c in ("source", "drain"):
            if self.contact_bc[c] != "schottky":
                continue
            ids = self.cnodes[c]
            vp = self.D["vR_p"] * self.T_vdW * (1 + Gp[ids])
            vn = self.D["vR_n"] * self.T_vdW * (1 + Gn[ids])
            Jps = q * vp * (p[ids] - peq)
            Jns = q * vn * (n[ids] - neq)
            P_if += float(np.sum((Jps * Vt * np.log(p[ids] / peq)
                                  + Jns * Vt * np.log(n[ids] / neq)) * csa[ids]))
        return dict(P_interface=P_if)

    def centerline(self):
        """Profiles along y = t/2 (2D) or the whole bar (1D)."""
        if self.kind == "1d":
            idx = np.argsort(self.x)
        else:
            ys = np.unique(self.y)
            yc = ys[np.argmin(np.abs(ys - 0.5 * self.geo["t"]))]
            idx = np.where(np.abs(self.y - yc) < 1e-12)[0]
            idx = idx[np.argsort(self.x[idx])]
        out = dict(x=self.x[idx])
        for nmn in ("Potential", "Holes", "Electrons"):
            out[nmn] = self.nodev(nmn)[idx]
        out["phip"] = out["Potential"] + self.D["Vt"] * np.log(out["Holes"] / self.D["ni"])
        if self.thermal:
            out["T"] = self.nodev("T")[idx]
        return out

    def decompose(self):
        """Decompose R_total using the hole quasi-Fermi potential phi_p on the centre line.

        Exact potential-drop decomposition (primary):
            R_c,source = (phi_p(x_ch0) - V_source)/I,  R_c,drain = (V_drain - phi_p(x_ch1))/I,
            R_ch = (phi_p(x_ch1) - phi_p(x_ch0))/I,     R_tot = R_c,s + R_ch + R_c,d (exactly),
        where x_ch0/x_ch1 are the inner contact edges (top) or the contact faces (edge). R_c
        therefore contains the interface quasi-Fermi jump plus the under-contact access path.
        R_if = P_interface / I^2 isolates the pure interface (thermionic/tunnel) part.

        Slope method (diagnostic only): R_ch,slope = |dphi_p/dx| * L_ch from the central
        50 % of the channel. In the injection (barrier-free) regime the hole density is
        strongly non-uniform along the channel and this estimate exceeds R_tot, i.e. a
        contact/channel split by linear extrapolation is not meaningful there."""
        I, Ip, In = self.current()
        V = self.bias
        cl = self.centerline()
        x0, x1 = self.geo["x_ch0"], self.geo["x_ch1"]
        L = x1 - x0
        phi = cl["phip"]
        f0 = float(np.interp(x0, cl["x"], phi))
        f1 = float(np.interp(x1, cl["x"], phi))
        R_tot = V / I
        Rcs = (f0 - 0.0) / I
        Rcd = (V - f1) / I
        Rch = (f1 - f0) / I
        mm = (cl["x"] >= x0 + 0.25 * L) & (cl["x"] <= x0 + 0.75 * L)
        slope = np.polyfit(cl["x"][mm], phi[mm], 1)[0]
        pi = self._interface_power()
        R_if = pi["P_interface"] / I ** 2 if "P_interface" in pi else float("nan")
        return dict(R_tot=R_tot, R_ch=Rch, R_c=0.5 * (Rcs + Rcd), R_c_source=Rcs,
                    R_c_drain=Rcd, frac_contact=(Rcs + Rcd) / R_tot, R_if=R_if,
                    R_ch_slope=abs(slope) * L / abs(I))


def params_table():
    """Parameter provenance table (R3.7). Category: M = experimentally measured on bulk
    or thick flakes, L = literature-derived (calculation / compilation), R = representative
    value chosen within a literature range, A = assumed (no direct source)."""
    return [
        ("E_g", "1.30 eV", "M (bulk 1.26 eV) -> R for 5 nm", "+ (larger E_g at 5 nm raises E_V, shifts threshold up)"),
        ("chi", "3.70 eV", "L", "+/- shifts E_V and the barrier-collapse threshold 1:1"),
        ("m_e,dos", "0.07 m0", "L", "small (minority carriers)"),
        ("m_h,dos", "0.17 m0", "L", "N_V, v_R, Richardson constant"),
        ("m_h,tun", "0.17 m0", "A (= DOS mass)", "heavier mass -> smaller Gamma"),
        ("eps_r", "6.4", "L", "few-layer screening may lower it"),
        ("mu_n", "1000 cm2/Vs", "M (encapsulated thin flakes)", "minor (holes dominate)"),
        ("mu_p", "50 cm2/Vs", "A (order of magnitude)", "R_ch ~ 1/mu_p"),
        ("mu(T) exponent", "1.5", "A (phonon-limited)", "self-heating feedback"),
        ("tau_n, tau_p", "1 ns", "A", "minor"),
        ("N_A", "1e16 cm-3", "A", "R_ch, depletion, Gamma"),
        ("kappa", "0.085 W/cmK", "R (few-layer range)", "T_max"),
        ("h_sub", "1.35e3 W/cm2K", "A (R_TBR=1e-8 m2K/W + 90 nm SiO2)", "T_max ~ 1/h_sub"),
        ("d_vdW, Phi_vdW", "0.30 nm, 2.0 eV", "A (typical vdW-gap tunnel barrier)", "top-contact T_vdW"),
        ("S, Phi_CNL", "1.0, 4.35 eV", "A (Schottky-Mott baseline; pinning varied)", "barrier height"),
    ]
