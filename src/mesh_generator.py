from devsim import *

def create_top_contact_mesh(mesh_name, device_name, L_c, L_gap, thickness, mesh_factor=1.0, graded_corner=False):
    """
    Creates a 2D mesh for a top contact structure.
    L_c: Contact length (horizontal) in cm
    L_gap: Gap between contacts in cm
    thickness: InSe thickness in cm
    mesh_factor: Multiplier for mesh density (higher = denser)
    graded_corner: Whether to refine the metal-semiconductor-vacuum triple points
    """
    create_2d_mesh(mesh=mesh_name)
    
    total_length = 2 * L_c + L_gap
    
    # X coordinates (horizontal)
    gap_spacing = 2.5e-7 / mesh_factor # Fixed 2.5 nm absolute spacing for channel gap
    
    add_2d_mesh_line(mesh=mesh_name, dir="x", pos=0.0, ps=(L_c/10)/mesh_factor)
    
    if graded_corner:
        corner_ps = 1e-8 / mesh_factor # 0.1 nm spacing at the corner
        # Left contact inner corner
        add_2d_mesh_line(mesh=mesh_name, dir="x", pos=L_c - 1e-7, ps=gap_spacing)
        add_2d_mesh_line(mesh=mesh_name, dir="x", pos=L_c, ps=corner_ps)
        add_2d_mesh_line(mesh=mesh_name, dir="x", pos=L_c + 1e-7, ps=gap_spacing)
        
        # Right contact inner corner
        add_2d_mesh_line(mesh=mesh_name, dir="x", pos=L_c + L_gap - 1e-7, ps=gap_spacing)
        add_2d_mesh_line(mesh=mesh_name, dir="x", pos=L_c + L_gap, ps=corner_ps)
        add_2d_mesh_line(mesh=mesh_name, dir="x", pos=L_c + L_gap + 1e-7, ps=gap_spacing)
    else:
        add_2d_mesh_line(mesh=mesh_name, dir="x", pos=L_c, ps=gap_spacing)
        add_2d_mesh_line(mesh=mesh_name, dir="x", pos=L_c + L_gap, ps=gap_spacing)
        
    add_2d_mesh_line(mesh=mesh_name, dir="x", pos=total_length, ps=(L_c/10)/mesh_factor)
    
    add_2d_mesh_line(mesh=mesh_name, dir="y", pos=-thickness, ps=(thickness/5)/mesh_factor) # Dummy space above
    if graded_corner:
        add_2d_mesh_line(mesh=mesh_name, dir="y", pos=0.0, ps=1e-8/mesh_factor)
        add_2d_mesh_line(mesh=mesh_name, dir="y", pos=1e-7, ps=(thickness/5)/mesh_factor)
    else:
        add_2d_mesh_line(mesh=mesh_name, dir="y", pos=0.0, ps=(thickness/5)/mesh_factor)
    add_2d_mesh_line(mesh=mesh_name, dir="y", pos=thickness, ps=(thickness/5)/mesh_factor)
    add_2d_mesh_line(mesh=mesh_name, dir="y", pos=2*thickness, ps=(thickness/5)/mesh_factor) # Dummy space below
    
    # Regions
    add_2d_region(mesh=mesh_name, material="InSe", region="bulk", xl=0.0, xh=total_length, yl=0.0, yh=thickness)
    add_2d_region(mesh=mesh_name, material="gas", region="air", xl=0.0, xh=total_length, yl=-thickness, yh=0.0)
    add_2d_region(mesh=mesh_name, material="metal", region="source_metal", xl=0.0, xh=L_c, yl=-thickness, yh=0.0)
    add_2d_region(mesh=mesh_name, material="metal", region="drain_metal", xl=L_c + L_gap, xh=total_length, yl=-thickness, yh=0.0)
    add_2d_region(mesh=mesh_name, material="metal", region="substrate_metal", xl=0.0, xh=total_length, yl=thickness, yh=2*thickness)
    
    # Add Contacts on the top surface (y=0)
    add_2d_contact(mesh=mesh_name, name="source", material="metal", region="bulk", xl=0.0, xh=L_c, yl=0.0, yh=0.0, bloat=1e-10)
    add_2d_contact(mesh=mesh_name, name="drain", material="metal", region="bulk", xl=L_c + L_gap, xh=total_length, yl=0.0, yh=0.0, bloat=1e-10)
    
    # Add Thermal Contact on the bottom surface (y=thickness)
    add_2d_contact(mesh=mesh_name, name="substrate", material="metal", region="bulk", xl=0.0, xh=total_length, yl=thickness, yh=thickness, bloat=1e-10)
    
    finalize_mesh(mesh=mesh_name)
    create_device(mesh=mesh_name, device=device_name)


def create_edge_contact_mesh(mesh_name, device_name, L_gap, thickness, mesh_factor=1.0):
    """
    Creates a 2D mesh for an edge contact structure.
    L_gap: Gap between contacts (channel length) in cm
    thickness: InSe thickness in cm
    mesh_factor: Multiplier for mesh density (higher = denser)
    For edge contact, the contacts are on the left and right vertical edges.
    """
    create_2d_mesh(mesh=mesh_name)
    
    # X coordinates (horizontal)
    gap_spacing = 2.5e-7 / mesh_factor # Fixed 2.5 nm absolute spacing for channel gap
    
    L_dummy = thickness
    add_2d_mesh_line(mesh=mesh_name, dir="x", pos=-L_dummy, ps=(L_dummy/5)/mesh_factor)
    add_2d_mesh_line(mesh=mesh_name, dir="x", pos=0.0, ps=gap_spacing)
    add_2d_mesh_line(mesh=mesh_name, dir="x", pos=L_gap, ps=gap_spacing)
    add_2d_mesh_line(mesh=mesh_name, dir="x", pos=L_gap+L_dummy, ps=(L_dummy/5)/mesh_factor)
    
