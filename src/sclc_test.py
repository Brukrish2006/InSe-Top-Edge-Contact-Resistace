import sys
import os
import csv
mkl_path = r"C:\Python314\Library\bin\mkl_rt.3.dll"
os.environ["DEVSIM_MATH_LIBS"] = mkl_path
import devsim
from mesh_generator import create_top_contact_mesh, create_edge_contact_mesh
from inse_physics import SetInSeParameters, CreateSchottkyContactDriftDiffusion

def run_simulation(device_name, geom, L_c, L_gap, thickness, wf, mu_p, kappa, mesh_factor, v_start, v_end, v_step):
    if geom == "top":
        create_top_contact_mesh(device_name, device_name, L_c, L_gap, thickness, mesh_factor=mesh_factor)
    else:
        create_edge_contact_mesh(device_name, device_name, L_gap, thickness, mesh_factor=mesh_factor)
        
    region = "bulk"
    
    devsim.set_parameter(device=device_name, region=region, name="NetDoping", value=1e16)
    SetInSeParameters(device_name, region)
    
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
    
    v_ds = v_start
    currents = []
    voltages = []
    resistances = []
    
    devsim.set_parameter(device=device_name, name="source_bias", value=0.0)
    
    while v_ds <= v_end + 1e-9:
        devsim.set_parameter(device=device_name, name="drain_bias", value=v_ds)
        devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=100)
        
        elec_I = devsim.get_contact_current(device=device_name, contact="drain", equation="ElectronContinuityEquation")
        hole_I = devsim.get_contact_current(device=device_name, contact="drain", equation="HoleContinuityEquation")
        total_I = elec_I + hole_I
        
        currents.append(abs(total_I))
        voltages.append(v_ds)
        
        if total_I != 0:
            resistances.append(abs(v_ds / total_I))
        else:
            resistances.append(float('inf'))
            
        v_ds += v_step
        
    return voltages, currents, resistances

def main():
    L_c = 50e-7      
    thickness = 5e-7 
    wf = 5.0
    mu = 50.0
    kappa = 0.085
    L_gap = 100e-7
    
    # 1. SCLC test: V_DS sweep for Top contact (L=100nm, mesh=1.0)
    print("Running SCLC V_DS sweep test for Top...")
    devsim.reset_devsim()
    v, i_top, r_top_sclc = run_simulation("dev_top_sclc", "top", L_c, L_gap, thickness, wf, mu, kappa, 1.0, 0.02, 0.5, 0.02)
    
    # 1b. SCLC test: V_DS sweep for Edge contact (L=100nm, mesh=1.0)
    print("Running SCLC V_DS sweep test for Edge...")
    devsim.reset_devsim()
    v, i_edge, r_edge_sclc = run_simulation("dev_edge_sclc", "edge", thickness, L_gap, thickness, wf, mu, kappa, 1.0, 0.02, 0.5, 0.02)
    
    with open("sclc_vds_sweep.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["V_DS", "R_Top", "R_Edge"])
        for volt, r_top, r_edge in zip(v, r_top_sclc, r_edge_sclc):
            writer.writerow([volt, r_top, r_edge])
            
    # 2. Mesh refinement test at L=100nm, V=0.1V, for both Top and Edge
    print("Running Fine Mesh Test (mesh=2.0)...")
    devsim.reset_devsim()
    _, _, r_top_fine = run_simulation("dev_top_fine", "top", L_c, L_gap, thickness, wf, mu, kappa, 2.0, 0.1, 0.1, 0.1)
    
    devsim.reset_devsim()
    _, _, r_edge_fine = run_simulation("dev_edge_fine", "edge", thickness, L_gap, thickness, wf, mu, kappa, 2.0, 0.1, 0.1, 0.1)
    
    print(f"Fine Mesh Results (L=100nm, V_DS=0.1V, Mesh Factor=2.0):")
    print(f"R_top: {r_top_fine[0]}")
    print(f"R_edge: {r_edge_fine[0]}")

    # 3. SCLC test at standard L=50nm length, wf=5.2 eV to check 17-28% claim
    print("Running SCLC V_DS sweep test for Top (L_c=50nm, L=50nm, wf=5.2)...")
    L_gap_50 = 50e-7
    wf_52 = 5.2
    devsim.reset_devsim()
    v_50, _, r_top_sclc_50 = run_simulation("dev_top_sclc_50_52", "top", L_c, L_gap_50, thickness, wf_52, mu, kappa, 1.0, 0.02, 0.5, 0.02)
    
    print("Running SCLC V_DS sweep test for Edge (L=50nm, wf=5.2)...")
    devsim.reset_devsim()
    _, _, r_edge_sclc_50 = run_simulation("dev_edge_sclc_50_52", "edge", thickness, L_gap_50, thickness, wf_52, mu, kappa, 1.0, 0.02, 0.5, 0.02)
    
    print("Running SCLC V_DS sweep test for Top (L_c=10nm, L=50nm, wf=5.2)...")
    L_c_10 = 10e-7
    devsim.reset_devsim()
    _, _, r_top_sclc_50_lc10 = run_simulation("dev_top_sclc_50_lc10_52", "top", L_c_10, L_gap_50, thickness, wf_52, mu, kappa, 1.0, 0.02, 0.5, 0.02)
    
    with open("sclc_vds_sweep_50nm_wf52.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["V_DS", "R_Top_Lc50", "R_Edge", "R_Top_Lc10"])
        for volt, r_t_50, r_e, r_t_10 in zip(v_50, r_top_sclc_50, r_edge_sclc_50, r_top_sclc_50_lc10):
            writer.writerow([volt, r_t_50, r_e, r_t_10])

if __name__ == "__main__":
    main()
