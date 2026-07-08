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

def run_simulation(device_name, geom, L_c, L_gap, thickness, wf, mesh_factor):
    # Create Mesh
    if geom == "top":
        create_top_contact_mesh(device_name, device_name, L_c, L_gap, thickness, mesh_factor)
    else:
        create_edge_contact_mesh(device_name, device_name, L_gap, thickness, mesh_factor)
        
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
    
    CreateSchottkyContactDriftDiffusion(device=device_name, region=region, contact="source")
    CreateSchottkyContactDriftDiffusion(device=device_name, region=region, contact="drain")
    devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=100)
    
    # Voltage Sweep
    v_ds_target = 0.1
    v_ds = 0.0
    v_step = 0.1
    
    currents = []
    voltages = []
    
    devsim.set_parameter(device=device_name, name="source_bias", value=0.0)
    # devsim.set_parameter(device=device_name, name="source_workfunction", value=wf)
    # devsim.set_parameter(device=device_name, name="drain_workfunction", value=wf)
    
    while v_ds <= v_ds_target + 1e-9:
        devsim.set_parameter(device=device_name, name="drain_bias", value=v_ds)
        devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=100)
        
        elec_I = devsim.get_contact_current(device=device_name, contact="drain", equation="ElectronContinuityEquation")
        hole_I = devsim.get_contact_current(device=device_name, contact="drain", equation="HoleContinuityEquation")
        total_I = elec_I + hole_I
        
        currents.append(total_I)
        voltages.append(v_ds)
        v_ds += v_step
        
    # Extract Resistance at V_DS = 0.5V
    R = voltages[-1] / currents[-1] if currents[-1] != 0 else float('inf')
    
    node_count = len(devsim.get_node_model_values(device=device_name, region=region, name="Potential"))
    
    return R, node_count

def main():
    mkl_path = r"C:\Python314\Library\bin\mkl_rt.3.dll"
    os.environ["DEVSIM_MATH_LIBS"] = mkl_path
    
    # Study params
    L_c = 50e-7      # 50 nm
    L_gap = 50e-7    # 50 nm channel
    thickness = 5e-7 # 5 nm thick
    wf = 5.0         # 5.0 eV Work function (Ohmic regime)
    
    mesh_factors = [0.5, 1.0, 2.0, 4.0]
    
    results = []
    
    for mf in mesh_factors:
        devsim.reset_devsim()
        print(f"Running Top Contact, mesh_factor={mf}")
        R_top, nodes_top = run_simulation("dev_top", "top", L_c, L_gap, thickness, wf, mf)
        
        devsim.reset_devsim()
        print(f"Running Edge Contact, mesh_factor={mf}")
        R_edge, nodes_edge = run_simulation("dev_edge", "edge", L_c, L_gap, thickness, wf, mf)
        
        results.append((mf, nodes_top, R_top, nodes_edge, R_edge))
        
    with open("mesh_convergence.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["MeshFactor", "Nodes_Top", "R_Top", "Nodes_Edge", "R_Edge"])
        writer.writerows(results)
        
    # Plotting
    nodes_t = [r[1] for r in results]
    r_t = [r[2] for r in results]
    nodes_e = [r[3] for r in results]
    r_e = [r[4] for r in results]
    
    plt.figure(figsize=(8, 5))
    plt.plot(nodes_t, r_t, 'o-', label='Top Contact (Lc=50nm)')
    plt.plot(nodes_e, r_e, 's-', label='Edge Contact')
    plt.xlabel('Number of Mesh Nodes')
    plt.ylabel(r'Total Resistance ($\Omega\cdot$cm)')
    plt.title(r'Mesh Convergence Study at $\Phi_M = 5.0$ eV')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('fig_mesh_convergence.png', dpi=300)
    print("Mesh convergence study completed.")

if __name__ == "__main__":
    main()
