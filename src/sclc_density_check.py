import sys
import os
import csv
mkl_path = r"C:\Python314\Library\bin\mkl_rt.3.dll"
os.environ["DEVSIM_MATH_LIBS"] = mkl_path
import devsim
from mesh_generator import create_edge_contact_mesh
from inse_physics import SetInSeParameters, CreateSchottkyContactDriftDiffusion

def run_density_check():
    device_name = "dev_density_check"
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
    
    biases = [0.02, 0.1, 0.3, 0.5]
    results = {}
    
    for v_ds in biases:
        devsim.set_parameter(device=device_name, name="drain_bias", value=v_ds)
        devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=100)
        
        x_coords = devsim.get_node_model_values(device=device_name, region=region, name="x")
        y_coords = devsim.get_node_model_values(device=device_name, region=region, name="y")
        holes = devsim.get_node_model_values(device=device_name, region=region, name="Holes")
        
        center_holes = []
        for x, y, h in zip(x_coords, y_coords, holes):
            # Check near center of channel
            if abs(x - L_gap/2) < 2e-7:
                center_holes.append(h)
                
        avg_center_hole = sum(center_holes) / len(center_holes) if center_holes else 0
        results[v_ds] = avg_center_hole
        
    for v_ds, h_dens in results.items():
        print(f"V_DS = {v_ds} V -> Avg center hole density: {h_dens:.2e} cm^-3")

    with open("sclc_density_check.csv", "w", newlin