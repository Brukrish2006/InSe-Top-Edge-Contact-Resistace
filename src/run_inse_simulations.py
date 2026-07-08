import sys
import os
import csv
import math

# Automatically fix the DEVSIM Math library missing DLL error for VS Code
os.environ['DEVSIM_MATH_LIBS'] = r"C:\Python314\Library\bin\mkl_rt.3.dll"

from devsim import *
from mesh_generator import create_top_contact_mesh, create_edge_contact_mesh
import inse_physics
import simple_physics
import simple_dd

def run_simulation(device_name, geometry, L_c, L_gap, thickness, work_function, filename):
    region = "bulk"
    
    # 1. Mesh Generation
    if geometry == 'top':
        create_top_contact_mesh(device_name, device_name, L_c, L_gap, thickness)
    elif geometry == 'edge':
        create_edge_contact_mesh(device_name, device_name, L_gap, thickness)
        
    # 2. Physics setup
    set_parameter(device=device_name, region=region, name="NetDoping", value=1e16)
    
    inse_physics.SetInSeParameters(device=device_name, region=region, T=300)
    
    # Create Potential, Electrons, Holes variables and basic equations
    simple_physics.CreateSiliconPotentialOnly(device=device_name, region=region)
    
    # Create Thermal variables and equations
    inse_physics.CreateLatticeTemperature(device=device_name, region=region)
    
    # Create Contacts (Potential only for equilibrium)
    inse_physics.CreateSchottkyContactPotential(device=device_name, region=region, contact="source", work_function=work_function)
    inse_physics.CreateSchottkyContactPotential(device=device_name, region=region, contact="drain", work_function=work_function)
    inse_physics.CreateThermalContact(device=device_name, region=region, contact="substrate")
    
    # Initial Solution (Equilibrium)
    set_parameter(device=device_name, name="source_bias", value=0.0)
    set_parameter(device=device_name, name="drain_bias", value=0.0)
    
    solve(type="dc", absolute_error=1.0, relative_error=1e-5, maximum_iterations=100)
    
    # Initialize Drift-Diffusion Variables
    from devsim.python_packages.model_create import CreateSolution
    CreateSolution(device_name, region, "Electrons")
    CreateSolution(device_name, region, "Holes")
    set_node_values(device=device_name, region=region, name="Electrons", init_from="IntrinsicElectrons")
    set_node_values(device=device_name, region=region, name="Holes", init_from="IntrinsicHoles")
    
    # Add Drift-Diffusion
    simple_dd.CreateBernoulli(device=device_name, region=region)
    simple_dd.CreateElectronCurrent(device=device_name, region=region, mu_n="(mu_n@n0 + mu_n@n1)/2")
    simple_dd.CreateHoleCurrent(device=device_name, region=region, mu_p="(mu_p@n0 + mu_p@n1)/2")
    simple_physics.CreateSRH(device=device_name, region=region)
    simple_physics.CreateECE(device=device_name, region=region, mu_n="(mu_n@n0 + mu_n@n1)/2")
    simple_physics.CreateHCE(device=device_name, region=region, mu_p="(mu_p@n0 + mu_p@n1)/2")
    
    # Add Drift-Diffusion to Contacts
    inse_physics.CreateSchottkyContactDriftDiffusion(device=device_name, region=region, contact="source")
    inse_physics.CreateSchottkyContactDriftDiffusion(device=device_name, region=region, contact="drain")
    
    import devsim
    print("Available edge models:", devsim.get_edge_model_list(device=device_name, region="bulk"))
    solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=100)
    
    # Voltage Sweep
    v_step = 0.1
    v_max = 2.0  # Increased to 2.0V to induce more self-heating
    v_drain = 0.0
    
    results = []
    
    while v_drain <= v_max + 1e-9:
        set_parameter(device=device_name, name="drain_bias", value=v_drain)
        solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=100)
        
        # Get currents
        i_drain_e = get_contact_current(device=device_name, contact="drain", equation="ElectronContinuityEquation")
        i_drain_h = get_contact_current(device=device_name, contact="drain", equation="HoleContinuityEquation")
        i_drain = i_drain_e + i_drain_h
        
        # Get max temperature
        T_nodes = get_node_model_values(device=device_name, region=region, name="LatticeTemperature")
        max_T = max(T_nodes) if T_nodes else 300.0
        
        results.append((v_drain, i_drain, max_T))
        
        v_drain += v_step
        
    # Write to CSV
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['V_drain (V)', 'I_drain (A/cm)', 'Max_T (K)'])
        for r in results:
            writer.writerow(r)
            
    # Calculate Resistance at V=0.1V (linear region)
    if len(results) > 1 and results[1][1] != 0:
        R_total = abs(results[1][0] / results[1][1])
    else:
        R_total = float('inf')
        
    final_max_T = results[-1][2] if results else 300.0
        
    return R_total, final_max_T

def main():
    L_gap = 50e-7 # 50 nm
    thickness = 5e-7 # 5 nm
    
    geometries = ['top', 'edge']
    L_c_values = [10e-7, 20e-7, 30e-7, 50e-7]
    work_functions = [4.33, 4.8, 5.0, 5.2]
    
    summary_file = 'summary_results_thermal.csv'
    
    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Geometry', 'L_c (cm)', 'Work Function (eV)', 'R_total (Ohm*cm)', 'Max_T at 2.0V (K)'])
        
    counter = 0
    for geom in geometries:
        for wf in work_functions:
            if geom == 'top':
                for L_c in L_c_values:
                    counter += 1
                    device_name = f"dev_{counter}"
                    print(f"Running: geom={geom}, L_c={L_c}, WF={wf}")
                    fname = f"iv_{geom}_Lc{L_c*1e7:.0f}nm_WF{wf}eV.csv"
                    try:
                        reset_devsim()
                        R, max_T = run_simulation(device_name, geom, L_c, L_gap, thickness, wf, fname)
                        with open(summary_file, 'a', newline='') as f:
                            csv.writer(f).writerow([geom, L_c, wf, R, max_T])
                    except Exception as e:
                        import traceback
                        traceback.print_exc()
                        print(f"Failed: {e}")
            elif geom == 'edge':
                counter += 1
                device_name = f"dev_{counter}"
                L_c = thickness
                print(f"Running: geom={geom}, L_c=thickness, WF={wf}")
                fname = f"iv_{geom}_WF{wf}eV.csv"
                try:
                    reset_devsim()
                    R, max_T = run_simulation(device_name, geom, L_c, L_gap, thickness, wf, fname)
                    with open(summary_file, 'a', newline='') as f:
                        csv.writer(f).writerow([geom, 'N/A', wf, R, max_T])
                except Exception as e:
                    import traceback
                    traceback.print_exc()
                    print(f"Failed: {e}")

if __name__ == "__main__":
    main()
