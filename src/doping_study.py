import sys
import math
import csv
import matplotlib.pyplot as plt
import os

mkl_path = r"C:\Python314\Library\bin\mkl_rt.3.dll"
os.environ["DEVSIM_MATH_LIBS"] = mkl_path

import devsim
from mesh_generator import create_top_contact_mesh, create_edge_contact_mesh
from inse_physics import SetInSeParameters, CreateSchottkyContactDriftDiffusion

def run_simulation(device_name, geom, doping, wf):
    L_c = 50e-7
    L_gap = 50e-7
    thickness = 5e-7
    
    if geom == "top":
        create_top_contact_mesh(device_name, device_name, L_c, L_gap, thickness)
    else:
        create_edge_contact_mesh(device_name, device_name, L_gap, thickness)
        
    region = "bulk"
    
    devsim.set_parameter(device=device_name, region=region, name="NetDoping", value=doping)
    SetInSeParameters(device_name, region)
    
    from simple_physics import CreateSiliconPotentialOnly, CreateSRH, CreateECE, CreateHCE
    import simple_dd
    CreateSiliconPotentialOnly(device_name, region)
    
    from inse_physics import CreateLatticeTemperature, CreateSchottkyContactPotential, CreateThermalContact
    CreateLatticeTemperature(device_name, region)
    CreateSchottkyContactPotential(device_name, region, "source", wf)
    CreateSchottkyContactPotential(device_name, region, "drain", wf)
    CreateThermalContact(device_name, region, "substrate")
    
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
    return abs(R)

def main():
    wf_values = [4.8, 5.0, 5.2]
    doping_values = [1e15, 1e16, 1e17]
    results = []
    
    for wf in wf_values:
        for dop in doping_values:
            devsim.reset_devsim()
            print(f"Top Contact: dop={dop}, WF={wf}")
            R_top = run_simulation(f"dev_top_wf{wf}_dop{dop}", "top", dop, wf)
            
            devsim.reset_devsim()
            print(f"Edge Contact: dop={dop}, WF={wf}")
            R_edge = run_simulation(f"dev_edge_wf{wf}_dop{dop}", "edge", dop, wf)
            
            results.append((wf, dop, R_top, R_edge))
            
    with open("doping_study.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["WF_eV", "Doping_cm3", "R_Top", "R_Edge"])
        writer.writerows(results)
        
    print("Doping study completed.")

if __name__ == "__main__":
    main()
