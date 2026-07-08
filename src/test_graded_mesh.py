import os
import csv
import sys
import time

mkl_path = r"C:\Python314\Library\bin\mkl_rt.3.dll"
os.environ["DEVSIM_MATH_LIBS"] = mkl_path
import devsim

from mesh_generator import create_top_contact_mesh
from inse_physics import SetInSeParameters, CreateSchottkyContactDriftDiffusion
from simple_physics import CreateSiliconPotentialOnly, CreateSRH, CreateECE, CreateHCE
import simple_dd
from devsim.python_packages.model_create import CreateSolution

def run_graded_mesh_test():
    """
    Tests whether refining the metal-semiconductor-vacuum triple-point corner
    (graded_corner=True in mesh_generator.py) changes the top-contact
    mesh-convergence behavior reported in Section 3.6's uniform-mesh sweep.

    NOTE: current is extracted as the sum of electron + hole contact current
    (ElectronContinuityEquation + HoleContinuityEquation), matching the
    methodology used everywhere else in this repository
    (mesh_convergence_study.py, run_inse_simulations.py). An earlier version
    of this script extracted hole current alone; at Phi_M=5.0 eV the electron
    contribution turns out to be numerically negligible here, so this fix
    does not change the resulting R values, but the corrected form is kept
    for methodological consistency with the rest of the pipeline.

    Extended to mesh_factor=4.0 and 8.0 (beyond the original 0.5-2.0 range)
    because the first three points alone are insufficient to establish
    whether the graded mesh actually plateaus.
    """
    geom = "top"
    L_c = 50e-7
    L_gap = 50e-7
    thickness = 5e-7
    wf = 5.0
    mu = 50.0
    kappa = 0.085
    v_ds = 0.1
    mesh_factors = [0.5, 1.0, 2.0, 4.0, 8.0]

    results = []

    for mf in mesh_factors:
        devsim.reset_devsim()
        device_name = f"dev_{geom}_{mf}"
        region = "bulk"

        # USE GRADED MESH
        create_top_contact_mesh(device_name, device_name, L_c, L_gap, thickness, mesh_factor=mf, graded_corner=True)

        devsim.set_parameter(device=device_name, region=region, name="NetDoping", value=1e16)
        SetInSeParameters(device_name, region)
        devsim.set_parameter(device=device_name, region=region, name="mu_p_300", value=mu)
        devsim.set_parameter(device=device_name, region=region, name="ThermalConductivity", value=kappa)

        CreateSiliconPotentialOnly(device_name, region)

        from inse_physics import CreateLatticeTemperature, CreateSchottkyContactPotential, CreateThermalContact
        CreateLatticeTemperature(device_name, region)
        CreateSchottkyContactPotential(device_name, region, "source", wf)
        CreateSchottkyContactPotential(device_name, region, "drain", wf)
        CreateThermalContact(device_name, region, "substrate")

        devsim.set_parameter(device=device_name, name="source_bias", value=0.0)
        devsim.set_parameter(device=device_name, name="drain_bias", value=0.0)

        devsim.solve(type="dc", absolute_error=1.0, relative_error=1e-5, maximum_iterations=100)

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

        devsim.set_parameter(device=device_name, name="drain_bias", value=v_ds)
        devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=100)

        elec_current = devsim.get_contact_current(device=device_name, contact="drain", equation="ElectronContinuityEquation")
        hole_current = devsim.get_contact_current(device=device_name, contact="drain", equation="HoleContinuityEquation")
        total_current = elec_current + hole_current
        abs_current = abs(total_current)
        if abs_current == 0:
            resistance_ohm_cm = float('inf')
        else:
            resistance_ohm_cm = (v_ds / abs_current)

        nodes = devsim.get_node_model_values(device=device_name, region=region, name="x")
        num_nodes = len(nodes)

        print(f"Mesh Factor {mf}: Nodes = {num_nodes}, R = {resistance_ohm_cm:.4f} Ohm-cm")
        results.append({"mesh_factor": mf, "nodes": num_nodes, "R_ohm_cm": resistance_ohm_cm})

    with open("graded_mesh_test.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["mesh_factor", "nodes", "R_ohm_cm"])
        writer.writeheader()
        for row in results:
            writer.writerow(row)

if __name__ == "__main__":
    run_graded_mesh_test()
