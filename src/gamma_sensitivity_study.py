import sys
import math
import csv
import matplotlib.pyplot as plt
import os

mkl_path = r"C:\Python314\Library\bin\mkl_rt.3.dll"
os.environ["DEVSIM_MATH_LIBS"] = mkl_path

import devsim

from mesh_generator import create_top_contact_mesh, create_edge_contact_mesh
from simple_physics import SetSiliconParameters
from inse_physics import SetInSeParameters, CreateSchottkyContactDriftDiffusion

def run_simulation(device_name, geom, L_c, L_gap, thickness, wf, Gamma, T_vdW):
    if geom == "top":
        create_top_contact_mesh(device_name, device_name, L_c, L_gap, thickness)
    else:
        create_edge_contact_mesh(device_name, device_name, L_gap, thickness)
        
    region = "bulk"
    
    # Physics setup
    devsim.set_parameter(device=device_name, region=region, name="NetDoping", value=1e16)
    SetInSeParameters(device_name, region)
    
    from simple_physics import CreateSiliconPotentialOnly, CreateSRH, CreateECE, CreateHCE
    import simple_dd
    
    CreateSiliconPotentialOnly(device_name, region)
    
    from inse_physics import CreateLatticeTemperature, CreateSchottkyContactPotential, CreateThermalContact
    CreateLatticeTemperature(device_name, region)
    
    CreateSchottkyContactPotential(device_name, region, "source", wf)
    CreateSchottkyContactPotential(device_name, region, "drain", wf)
    CreateThermalContact(device_name, region, "substrate")
    
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
    CreateECE(device_name, region, mu_n="(mu_n@n0 + mu_n@n1)/2")
    CreateHCE(device_name, region, mu_p="(mu_p@n0 + mu_p@n1)/2")
    
    CreateSchottkyContactDriftDiffusion(device=device_name, region=region, contact="source", Gamma=Gamma, T_vdW=T_vdW)
    CreateSchottkyContactDriftDiffusion(device=device_name, region=region, contact="drain", Gamma=Gamma, T_vdW=T_vdW)
    
    devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=100)
    
    # Voltage Sweep
    v_ds_target = 0.1
    v_ds = 0.0
    v_step = 0.1
    
    currents = []
    voltages = []
    
    devsim.set_parameter(device=device_name, name="source_bias", value=0.0)
    
    while v_ds <= v_ds_target + 1e-9:
        devsim.set_parameter(device=device_name, name="drain_bias", value=v_ds)
        devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=100)
        
        elec_I = devsim.get_contact_current(device=device_name, contact="drain", equation="ElectronContinuityEquation")
        hole_I = devsim.get_contact_current(device=device_name, contact="drain", equation="HoleContinuityEquation")
        total_I = elec_I + hole_I
        
        currents.append(total_I)
        voltages.append(v_ds)
        v_ds += v_step
        
    R = voltages[-1] / currents[-1] if currents[-1] != 0 else float('inf')
    return abs(R)

def main():
    mkl_path = r"C:\Python314\Library\bin\mkl_rt.3.dll"
    os.environ["DEVSIM_MATH_LIBS"] = mkl_path
    
    L_c = 50e-7      
    L_gap = 50e-7    
    thickness = 5e-7 
    
    wfs = [4.33, 4.6, 4.8, 5.0, 5.2]
    gammas = [1.0, 10.0, 100.0]
    
    results_top = {g: [] for g in gammas}
    results_edge = {g: [] for g in gammas}
    
    for Gamma in gammas:
        for wf in wfs:
            devsim.reset_devsim()
            print(f"Top Contact: Gamma={Gamma}, WF={wf}")
            try:
                R_top = run_simulation(f"dev_top_g{Gamma}_wf{wf}", "top", L_c, L_gap, thickness, wf, Gamma, 1.0)
                results_top[Gamma].append(R_top)
            except Exception as e:
                print(f"Failed Top Gamma={Gamma}, WF={wf}: {e}")
                results_top[Gamma].append(float('inf'))
                
            devsim.reset_devsim()
            print(f"Edge Contact: Gamma={Gamma}, WF={wf}")
            try:
                R_edge = run_simulation(f"dev_edge_g{Gamma}_wf{wf}", "edge", thickness, L_gap, thickness, wf, Gamma, 1.0)
                results_edge[Gamma].append(R_edge)
            except Exception as e:
                print(f"Failed Edge Gamma={Gamma}, WF={wf}: {e}")
                results_edge[Gamma].append(float('inf'))
                
    with open("gamma_sensitivity.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["WF"] + [f"Top_Gamma={g}" for g in gammas] + [f"Edge_Gamma={g}" for g in gammas])
        for i, wf in enumerate(wfs):
            writer.writerow([wf] + [results_top[g][i] for g in gammas] + [results_edge[g][i] for g in gammas])
            
    # Plotting
    plt.figure(figsize=(8, 5))
    colors = ['r', 'g'