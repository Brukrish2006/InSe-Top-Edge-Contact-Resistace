import math
from devsim import set_parameter, contact_equation, get_node_model_values
from simple_physics import (
    CreateContactNodeModel, GetContactBiasName, GetContactNodeModelName, 
    InEdgeModelList, CreateEdgeModel, CreateEdgeModelDerivatives,
    q, k, eps_0
)

def SetInSeParameters(device, region, T=300):
    eps_inse = 6.4
    Eg = 1.3 # eV
    EA = 3.7 # eV
    m_e_dos = 0.07
    m_h_dos = 0.17
    mu_n = 1000.0
    mu_p = 50.0
    m0 = 9.10938356e-31
    h = 6.62607015e-34
    
    # Nc, Nv in cm^-3
    factor = 2 * ( (2 * math.pi * k * T / (h**2))**1.5 ) * 1e-6
    Nc = factor * (m_e_dos * m0)**1.5
    Nv = factor * (m_h_dos * m0)**1.5
    
    ni = math.sqrt(Nc * Nv) * math.exp(-Eg * q / (2 * k * T))
    
    set_parameter(device=device, region=region, name="Permittivity", value=eps_inse * eps_0)
    set_parameter(device=device, region=region, name="ElectronCharge", value=q)
    set_parameter(device=device, region=region, name="n_i", value=ni)
    set_parameter(device=device, region=region, name="T", value=T)
    set_parameter(device=device, region=region, name="kT", value=k*T)
    set_parameter(device=device, region=region, name="V_t", value=k*T/q)
    # We will make mobility temperature dependent in the equations, 
    # but store mu_n_300 and mu_p_300 as parameters
    set_parameter(device=device, region=region, name="mu_n_300", value=mu_n)
    set_parameter(device=device, region=region, name="mu_p_300", value=mu_p)
    
    # Thermal parameters
    set_parameter(device=device, region=region, name="ThermalConductivity", value=0.085) # W/cmK (8.5 W/mK)
    set_parameter(device=device, region=region, name="LatticeTemperature_Ambient", value=300.0)
    set_parameter(device=device, region=region, name="n1", value=ni)
    set_parameter(device=device, region=region, name="p1", value=ni)
    set_parameter(device=device, region=region, name="taun", value=1e-9)
    set_parameter(device=device, region=region, name="taup", value=1e-9)
    
    set_parameter(device=device, region=region, name="Nc", value=Nc)
    set_parameter(device=device, region=region, name="Nv", value=Nv)
    set_parameter(device=device, region=region, name="Eg", value=Eg)
    set_parameter(device=device, region=region, name="ElectronAffinity", value=EA)
    print(f"InSe Params Set: Nc={Nc:e}, Nv={Nv:e}, ni={ni:e}")

def CreateSchottkyContactPotential(device, region, contact, work_function, is_circuit=False):
    if not InEdgeModelList(device, region, "contactcharge_edge"):
        CreateEdgeModel(device, region, "contactcharge_edge", "Permittivity*ElectricField")
        CreateEdgeModelDerivatives(device, region, "contactcharge_edge", "Permittivity*ElectricField", "Potential")

    set_parameter(device=device, region=region, name=f"{contact}_work_function", value=work_function)

    # psi_contact = V_bias + ElectronAffinity - WorkFunction + Eg/2 - 0.5*V_t*log(Nc/Nv)
    contact_model = f"Potential - {GetContactBiasName(contact)} - (ElectronAffinity - {contact}_work_function + Eg/2 - 0.5*V_t*log(Nc/Nv))"
    
    contact_model_name = GetContactNodeModelName(contact)
    CreateContactNodeModel(device, contact, contact_model_name, contact_model)
    CreateContactNodeModel(device, contact, f"{contact_model_name}:Potential", "1")
    if is_circuit:
        CreateContactNodeModel(device, contact, f"{contact_model_name}:{GetContactBiasName(contact)}", "-1")

    if is_circuit:
        contact_equation(device=device, contact=contact, name="PotentialEquation", node_model=contact_model_name, edge_charge_model="contactcharge_edge", circuit_node=GetContactBiasName(contact))
    else:
        contact_equation(device=device, contact=contact, name="PotentialEquation", node_model=contact_model_name, edge_charge_model="contactcharge_edge")


def CreateSchottkyContactDriftDiffusion(device, region, contact, is_circuit=False, Gamma=10.0, T_vdW=1.0):
    from devsim import set_parameter
    
    # Effective mass based Richardson constants for InSe: A* = 120 * (m*/m0) A/cm^2 K^2
    set_parameter(device=device, region=region, name=f"{contact}_A_star_n", value=8.4)  # m_e = 0.07
    set_parameter(device=device, region=region, name=f"{contact}_A_star_p", value=20.4) # m_h = 0.17
    
    # Tunneling enhancement factor (Gamma) for Thermionic Field Emission
    set_parameter(device=device, region=region, name=f"{contact}_Gamma_n", value=Gamma) 
    set_parameter(device=device, region=region, name=f"{contact}_Gamma_p", value=Gamma)
    
    # van der Waals transmission coefficient
    set_parameter(device=device, region=region, name=f"{contact}_T_vdW", value=T_vdW)
    
    # Thermionic emission velocities
    v_th_n = f"({contact}_A_star_n * T^2) / (ElectronCharge * Nc)"
    v_th_p = f"({contact}_A_star_p * T^2) / (ElectronCharge * Nv)"
    
    # Equilibrium concentrations at the contact
    n_eq = f"Nc * exp( -( {contact}_work_function - ElectronAffinity ) / V_t )"
    p_eq = f"Nv * exp( -( ElectronAffinity + Eg - {contact}_work_function ) / V_t )"
    
    # Net Thermionic + Tunneling Current (from bulk into metal)
    contact_electrons_model = f"-ElectronCharge * ({v_th_n}) * {contact}_T_vdW * (1 + {contact}_Gamma_n) * (Electrons - ({n_eq}))"
    contact_holes_model = f"ElectronCharge * ({v_th_p}) * {contact}_T_vdW * (1 + {contact}_Gamma_p) * (Holes - ({p_eq}))"
    
    # Derivatives for the Newton solver
    d_electrons = f"-ElectronCharge * ({v_th_n}) * {contact}_T_vdW * (1 + {contact}_Gamma_n)"
    d_holes = f"ElectronCharge * ({v_th_p}) * {contact}_T_vdW * (1 + {contact}_Gamma_p)"

    contact_electrons_name = f"{contact}nodeelectrons"
    contact_holes_name = f"{contact}nodeholes"

    CreateContactNodeModel(device, contact, contact_electrons_name, contact_electrons_model)
    CreateContactNodeModel(device, contact, f"{contact_electrons_name}:Electrons", d_electrons)
    CreateContactNodeModel(device, contact, contact_holes_name, contact_holes_model)
    CreateContactNodeModel(device, contact, f"{contact_holes_name}:Holes", d_holes)

    if is_circuit:
        contact_equation(device=device, contact=contact, name="ElectronContinuityEquation", node_model=contact_electrons_name, edge_current_model="ElectronCurrent", circuit_node=GetContactBiasName(contact))
        contact_equation(device=device, contact=contact, name="HoleContinuityEquation", node_model=contact_holes_name, edge_current_model="HoleCurrent", circuit_node=GetContactBiasName(contact))
    else:
        contact_equation(device=device, contact=contact, name="ElectronContinuityEquation", node_model=contact_electrons_name, edge_current_model="ElectronCurrent")
        contact_equation(device=device, contact=contact, name="HoleContinuityEquation", node_model=contact_holes_name, edge_current_model="HoleCurrent")

def CreateLatticeTemperature(device, region):
    from devsim import equation
    from model_create import CreateSolution, CreateNodeModel, CreateNodeModelDerivative, CreateEdgeModel, CreateEdgeModelDerivatives, InNodeModelList
    
    # Create solution variable if not exists
    if not InNodeModelList(device, region, "LatticeTemperature"):
        CreateSolution(device, region, "LatticeTemperature")
        # Initialize to 300K
        from devsim import set_node_values
        CreateNodeModel(device, region, "T_init", "LatticeTemperature_Ambient")
        set_node_values(device=device, region=region, name="LatticeTemperature", init_from="T_init")

    # Temp-dependent mobility: mu(T) = mu_300 * (T/300)^-1.5
    mu_n_T = "mu_n_300 * pow(pow(LatticeTemperature / 300.0, 2.0), -0.75)"
    mu_p_T = "mu_p_300 * pow(pow(LatticeTemperature / 300.0, 2.0), -0.75)"
    CreateNodeModel(device, region, "mu_n", mu_n_T)
    CreateNodeModelDerivative(device, region, "mu_n", mu_n_T, "LatticeTemperature")
    CreateNodeModel(device, region, "mu_p", mu_p_T)
    CreateNodeModelDerivative(device, region, "mu_p", mu_p_T, "LatticeTemperature")
    
    from devsim import edge_from_node_model
    edge_from_node_model(device=device, region=region, node_model="mu_n")
    edge_from_node_model(device=device, region=region, node_model="mu_p")

    # Thermal flux edge model: F_th = kappa * grad(T)
    if not InEdgeModelList(device, region, "ElectronCurrent"):
        CreateEdgeModel(device, region, "ElectronCurrent", "0.0")
    if not InEdgeModelList(device, region, "HoleCurrent"):
        CreateEdgeModel(device, region, "HoleCurrent", "0.0")
        
    heat_flux = "ThermalConductivity * (LatticeTemperature@n0 - LatticeTemperature@n1) * EdgeInverseLength"
    CreateEdgeModel(device, region, "ThermalFlux", heat_flux)
    for v in ["LatticeTemperature"]:
        CreateEdgeModelDerivatives(device, region, "ThermalFlux", heat_flux, v)

    # Proper symmetric Joule Heating via edge_volume_model
    # Power generation per cm^3: J \cdot E = J \cdot \Delta V / L
    # We negate it so it acts as a heat source in the DEVSIM residual (which evaluates sum(flux) + node_models = 0)
    joule_heat_vol = "-(ElectronCurrent + HoleCurrent) * (Potential@n0 - Potential@n1) * EdgeInverseLength"
    CreateEdgeModel(device, region, "JouleHeat_vol", joule_heat_vol)
    for v in ["Potential", "Electrons", "Holes"]:
        CreateEdgeModelDerivatives(device, region, "JouleHeat_vol", joule_heat_vol, v)

    # Set up DEVSIM to integrate edge_volume_model symmetrically over the nodes
    CreateEdgeModel(device, region, "EdgeNodeVolume", "0.5 * EdgeCouple * EdgeLength")
    from devsim import set_parameter
    set_parameter(name="edge_node0_volume_model", value="EdgeNodeVolume")
    set_parameter(name="edge_node1_volume_model", value="EdgeNodeVolume")

    # Equation: div(ThermalFlux) + JouleHeat_vol = 0
    equation(device=device, region=region, name="LatticeTemperatureEquation", variable_name="LatticeTemperature",
             edge_model="ThermalFlux", edge_volume_model="JouleHeat_vol", variable_update="default")

def CreateThermalContact(device, region, contact, temperature="LatticeTemperature_Ambient"):
    # Fixed temperature boundary condition
    contact_model = f"LatticeTemperature - {temperature}"
    contact_model_name = f"{contact}_T_nodemodel"
    
    from simple_physics import CreateContactNodeModel
    CreateContactNodeModel(device, contact, contact_model_name, contact_model)
    CreateContactNodeModel(device, contact, f"{contact_model_name}:LatticeTemperature", "1")
    
    contact_equation(device=device, contact=contact, name="LatticeTemperatureEquation", node_model=contact_model_name)


