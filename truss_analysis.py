###############################################################################
# 2D TRUSS ANALYSIS USING THE FINITE ELEMENT METHOD (FEM)
# Author: Ayesha Ahmed
# Date: 2025
# Description:
#   This program performs static structural analysis of 2D truss structures
#   using the Finite Element Method. It calculates nodal displacements,
#   reaction forces, and axial member forces, and visualises both undeformed
#   and deformed truss geometries.
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# USER INPUT SECTION – DEFINE NODES AND GEOMETRY
###############################################################################
how_many_nodes = int(input("How many nodes? "))

if how_many_nodes <= 1:
    print("Insufficient nodes")
    exit()  # Program must have at least two nodes for an element

# Create a dictionary to store node coordinates
nodes = {1: (0.0, 0.0)}  # Node 1 fixed at origin by convention
print("Node 1 assumed at (0,0).")

# Request coordinates for remaining nodes
for this_node in range(2, how_many_nodes + 1):
    node_x = float(input(f"Enter x coordinate of node {this_node}: "))
    node_y = float(input(f"Enter y coordinate of node {this_node}: "))
    nodes[this_node] = (node_x, node_y)

# Display node coordinate library
print("\nNode Coordinate Library:")
for node, coord in nodes.items():
    print(f"Node {node}: ({coord[0]}, {coord[1]})")

###############################################################################
# DEFINE ELEMENT CONNECTIVITY, SUPPORTS, AND LOADS
###############################################################################
elements = {}
num_elements = int(input("How many elements?"))

# Define each element with start and end nodes, cross-section, and material
for e in range(1, num_elements + 1):
    i = int(input(f"Enter start node for element {e}: "))
    j = int(input(f"Enter end node for element {e}: "))
    A = float(input(f"Enter cross-sectional area (m^2) for element {e}: "))
    E = float(input(f"Enter Young's Modulus (Pa) for element {e}: "))
    elements[e] = (i, j, A, E)

# Define supports (fixed nodes)
supports = {}
print("Enter supports as: node direction(s)")
print("Example: 1 x,y  or  2 y  or  3 x")
print("Use commas between multiple directions.\n")

num_supports = int(input("How many supported nodes? "))

for n in range(num_supports):
    support_input = input("Enter support (e.g. '1 x,y' or '2 y'): ").split()
    node_id = int(support_input[0])
    dirs = tuple(d.strip() for d in support_input[1].split(","))
    supports[node_id] = dirs

# Display support summary
print("\nSupport summary:")
for node, dirs in supports.items():
    print(f"Node {node}: fixed in {', '.join(dirs)} direction(s)")


# Define applied loads
loads = {}
num_loads = int(input("\nHow many loaded nodes? "))

for load in range(num_loads):
    node_id = int(input("Enter node number where load is applied: "))
    Fx = float(input(f"Enter Fx (N) at node {node_id}: "))
    Fy = float(input(f"Enter Fy (N) at node {node_id}: "))
    loads[node_id] = (Fx, Fy)

# Display load summary
print("\nLoad summary:")
for node, (Fx, Fy) in loads.items():
    print(f"Node {node}: Fx = {Fx} N, Fy = {Fy} N")

###############################################################################
# FUNCTION: ELEMENT STIFFNESS MATRIX FORMULATION
###############################################################################
def element_stiffness(x1, y1, x2, y2, A, E):
    """
    Compute the element stiffness matrix for a 2D truss element.

    Inputs:
        x1, y1, x2, y2 : coordinates of the two nodes (m)
        A : cross-sectional area (m^2)
        E : Young's Modulus (Pa)

    Returns:
        k : 4x4 element stiffness matrix in global coordinates
    """
    L = ((x2 - x1)**2 + (y2 - y1)**2)**0.5     # Element length
    c = (x2 - x1) / L                          # cos(theta)
    s = (y2 - y1) / L                          # sin(theta)

    # Element stiffness matrix in global coordinates
    k = (A * E / L) * np.array([
        [ c*c,  c*s, -c*c, -c*s],
        [ c*s,  s*s, -c*s, -s*s],
        [-c*c, -c*s,  c*c,  c*s],
        [-c*s, -s*s,  c*s,  s*s]
    ])
    return k

###############################################################################
# GLOBAL STIFFNESS MATRIX ASSEMBLY
###############################################################################
ndof = 2 * len(nodes)                 # 2 DOFs per node (x and y)
K = np.zeros((ndof, ndof))            # Initialise global stiffness matrix

# Assemble element stiffness matrices into global matrix
for e, (i, j, A, E) in elements.items():
    x1, y1 = nodes[i]
    x2, y2 = nodes[j]
    k = element_stiffness(x1, y1, x2, y2, A, E)

    # Define DOF mapping for element nodes i and j
    dof = [2*(i-1), 2*(i-1)+1, 2*(j-1), 2*(j-1)+1]

    # Add element stiffness to global stiffness matrix
    for m in range(4):
        for n in range(4):
            K[dof[m], dof[n]] += k[m, n]

###############################################################################
# LOAD VECTOR FORMATION
###############################################################################
F = np.zeros(ndof)                    # Global load vector initialised to zero

# Populate load vector with applied nodal forces
for n, (Fx, Fy) in loads.items():
    F[2*(n-1)] = Fx                   # x-direction force
    F[2*(n-1)+1] = Fy                 # y-direction force

###############################################################################
# APPLY BOUNDARY CONDITIONS
###############################################################################
fixed_dofs = []

# Identify restrained degrees of freedom (x=0, y=1 per node)
for n, dirs in supports.items():
    for d in dirs:
        fixed_dofs.append(2*(n-1) + (0 if d == 'x' else 1))

# Free DOFs are all DOFs not in the fixed list
free_dofs = list(set(range(ndof)) - set(fixed_dofs))

###############################################################################
# SOLVE FOR NODAL DISPLACEMENTS
###############################################################################
# Partition global stiffness and load vector for free DOFs
K_ff = K[np.ix_(free_dofs, free_dofs)]
F_f = F[free_dofs]

# Solve for displacements at free DOFs
U = np.zeros(ndof)
U[free_dofs] = np.linalg.solve(K_ff, F_f)

###############################################################################
# POST-PROCESSING: CALCULATE ELEMENT AXIAL FORCES
###############################################################################
for e, (i, j, A, E) in elements.items():
    x1, y1 = nodes[i]
    x2, y2 = nodes[j]
    L = ((x2 - x1)**2 + (y2 - y1)**2)**0.5
    c = (x2 - x1) / L
    s = (y2 - y1) / L

    # Extract nodal displacements for this element
    dof = [2*(i-1), 2*(i-1)+1, 2*(j-1), 2*(j-1)+1]
    u_e = U[dof]

    # Transformation vector (relates global displacement to axial extension)
    T = np.array([-c, -s, c, s])

    # Compute axial force using F = (AE/L) * T * u_e
    F_axial = (A * E / L) * T.dot(u_e)
    print(f"Element {e}: Axial force = {F_axial:.2f} N")

###############################################################################
# VISUALISATION: UNDEFORMED VS DEFORMED STRUCTURE
###############################################################################
scale = 1000  # Scale factor to exaggerate deformation visually
plt.figure()

for e, (i, j, _, _) in elements.items():
    # Original (undeformed) coordinates
    x1, y1 = nodes[i]
    x2, y2 = nodes[j]
    plt.plot([x1, x2], [y1, y2], 'b--', label='Original' if e == 1 else "")  

    # Deformed coordinates
    u = U * scale
    x1d = x1 + u[2*(i-1)]
    y1d = y1 + u[2*(i-1)+1]
    x2d = x2 + u[2*(j-1)]
    y2d = y2 + u[2*(j-1)+1]
    plt.plot([x1d, x2d], [y1d, y2d], 'r-', label='Deformed' if e == 1 else "")

# Plot settings
plt.axis('equal')
plt.title("2D Truss Deformation (Blue = Original, Red = Deformed)")
plt.xlabel("X-coordinate (m)")
plt.ylabel("Y-coordinate (m)")
plt.legend()
plt.show()
