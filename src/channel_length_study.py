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

def run_simulation(device_name, geom, L_c, L_gap, thickness, wf, mu_p, kappa):
    if geom == "top":
        create_top_contact_mesh(device_name, device_name, L_c, L_gap, thickness)
    else:
        create_edge_contact_mesh(device_name, device_name, L_gap, thickness)
        
    region = "bulk"
    
    # Physics setup
    devsim.set_parameter(device=device_name, region=region, name="NetDoping", value=1e16)
    SetInSeParameters(device_name, region)
    
    # Override defaults
    devsim.set_parameter(device=device_name, region=region, name="mu_p_300", value=mu_p)
    devsim.set_parameter(device=device_name, region=region, name="ThermalConductivity", value=kappa)
    
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
    
    CreateSchottkyContactDriftDiffusion(device=device_name, region=region, contact="source")
    CreateSchottkyContactDriftDiffusion(device=device_name, region=region, contact="drain")
    
    devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=100)
    
    # Voltage Sweep up to 2.0V to capture thermal effects
    v_ds_target = 2.0
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
        
    T_nodes = devsim.get_node_model_values(device=device_name, region=region, name="LatticeTemperature")
    max_T = max(T_nodes) if T_nodes else 300.0
    
    # Resistance at lowest non-zero step (V = 0.1V)
    R = voltages[1] / currents[1] if len(currents) > 1 and currents[1] != 0 else float('inf')
    
    return abs(R), max_T

def main():
    L_c = 50e-7      
    thickness = 5e-7 
    wf = 5.0 # Near-Ohmic regime 
    mu = 50.0
    kappa = 0.085
    
    L_gap_values_nm = [10, 30, 50, 100]
    
    results = []
    
    for L_gap_nm in L_gap_values_nm:
        L_gap = L_gap_nm * 1e-7
        
        devsim.reset_devsim()
        print(f"Top Contact: L_gap={L_gap_nm} nm")
        R_top, T_top = run_simulation(f"dev_top_Lgap{L_gap_nm}", "top", L_c, L_gap, thickness, wf, mu, kappa)
        
        devsim.reset_devsim()
        print(f"Edge Contact: L_gap={L_gap_nm} nm")
        R_edge, T_edge = run_simulation(f"dev_edge_Lgap{L_gap_nm}", "edge", thickness, L_gap, thickness, wf, mu, kappa)
        
        results.append((L_gap_nm, R_top, R_edge, T_top, T_edge))
        
    with open("channel_length_study.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["L_gap_nm", "R_Top", "R_Edge", "MaxT_Top", "MaxT_Edge"])
        writer.writerows(results)
        
    # Plotting
    L_gaps = [r[0] for r in results]
    R_top_vals = [r[1] for r in results]
    R_edge_vals = [r[2] for r in results]
    
    plt.figure(figsize=(7, 5))
    plt.plot(L_gaps, R_top_vals, '^--', color='red', label='Top Contact (50 nm)')
    plt.plot(L_gaps, R_edge_vals, 'o-', color='blue', label='Edge Contact')
    plt.xlabel('Channel Length $L_{ch}$ (nm)')
    plt.ylabel(r'Total Resistance ($\Omega\cdot$cm)')
    plt.title(r'Total Resistance vs Channel Length ($\Phi_M$ = 5.0 eV)')
    plt.grid(True, ls="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig('fig_channel_length.png', dpi=300)
    
    print("Channel length study completed.")

if __name__ == "__main__":
    main()
