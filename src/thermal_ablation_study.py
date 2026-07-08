import sys
import math
import csv
import os

mkl_path = r"C:\Python314\Library\bin\mkl_rt.3.dll"
os.environ["DEVSIM_MATH_LIBS"] = mkl_path

import devsim
from mesh_generator import create_edge_contact_mesh
from inse_physics import SetInSeParameters, CreateSchottkyContactDriftDiffusion
from simple_physics import CreateSiliconPotentialOnly, CreateSRH, CreateECE, CreateHCE
import simple_dd
from inse_physics import CreateLatticeTemperature, CreateSchottkyContactPotential, CreateThermalContact

def run_ablation():
    device_name = "dev_ablation"
    region = "bulk"
    L_gap = 50e-7
    thickness = 5e-7
    wf = 5.0
    mu = 50.0
    kappa = 0.085
    
    create_edge_contact_mesh(device_name, device_name, L_gap, thickness, mesh_factor=1.0)
    
    devsim.set_parameter(device=device_name, region=region, name="NetDoping", value=1e16)
    SetInSeParameters(device_name, region)
    devsim.set_parameter(device=device_name, region=region, name="mu_p_300", value=mu)
    devsim.set_parameter(device=device_name, region=region, name="ThermalConductivity", value=kappa)
    
    CreateSiliconPotentialOnly(device_name, region)
    
    # Enable Thermal Equations but we'll see what happens
    CreateLatticeTemperature(device_name, region)
    CreateSchottkyContactPotential(device_name, region, "source", wf)
    CreateSchottkyContactPotential(device_name, region, "drain", wf)
    CreateThermalContact(device_name, region, "substrate")
    
    # ABLATION: Overwrite temperature-dependent mobility with constant mobility
    # The normal CreateECE/CreateHCE uses mu_n/mu_p which depends on T.
    # We will redefine mu_p and mu_n to ignore T.
    devsim.node_model(device=device_name, region=region, name="mu_p", equation=f"{mu}")
    devsim.node_model(device=device_name, region=region, name="mu_n", equation="mu_n_300")
    
    # Initial solve
    devsim.set_parameter(device=device_name, name="source_bias", value=0.0)
    devsim.set_parameter(device=device_name, name="drain_bias", value=0.0)
    devsim.solve(type="dc", absolute_error=1.0, relative_error=1e-5, maximum_iterations=100)
    
    from devsim.python_packages.model_create import CreateSolution
    CreateSolution(device_name, region, "Electrons")
    CreateSolution(device_name, region, "Holes")
    devsim.set_node_values(device=device_name, region=region, name="Electrons", init_from="IntrinsicElectrons")
    devsim.set_node_values(device=device_name, region=region, name="Holes", init_from="IntrinsicHoles")
    
    simple_dd.CreateBernoulli(device=device_name, region=region)
    simple_dd.CreateElectronCurrent(device=device_name, region=region, mu_n="(mu_n@n0 + mu_n@n1)/2")
    simple_dd.CreateHoleCurrent(device=device_name, region=region, mu_p="(mu_p@n0 + mu_p@n1)/2")
    CreateSRH(device_name, region)
    
    # We must use our custom mu strings here too
    CreateECE(device_name, region, mu_n="(mu_n@n0 + mu_n@n1)/2")
    CreateHCE(device_name, region, mu_p="(mu_p@n0 + mu_p@n1)/2")
    
    CreateSchottkyContactDriftDiffusion(device=device_name, region=region, contact="source")
    CreateSchottkyContactDriftDiffusion(device=device_name, region=region, contact="drain")
    
    devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=100)
    
    devsim.set_parameter(device=device_name, name="drain_bias", value=0.1)
    devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=100)
    
    elec_I = devsim.get_contact_current(device=device_name, contact="drain", equation="ElectronContinuityEquation")
    hole_I = devsim.get_contact_current(device=device_name, contact="drain", equation="HoleContinuityEquation")
    total_I = elec_I + hole_I
    
    R = 0.1 / total_I if total_I != 0 else float('inf')
    print(f"Ablated Resistance at 0.1V: {abs(R)} Ohm-cm")

if __name__ == "__main__":
    run_ablation()
