# Truss-Analysis

This project demonstrates the implementation of a 2D truss solver using Finite Element Analysis in Python. It computes nodal displacement, reaction forces and axial element forces for truss structures under external loads. The solution also provides a visual comparison between undeformed and deformed configurations.

Each element may be inclined at an angle theta.
To transform stiffness into the global coordinate system, we use the transformation matrix, T:
T= [c s 0 0]
   [0 0 c s]
c = cos (theta), s = sin (theta)
Global stiffness:
K = T^T  k  T
This yields:
k = AE/L  [c^2    cs     -c^2    -cs ]
          [cs     s^2    -cs     -s^2]
          [c^2    -cs    c^2     cs  ]
          [cs     s^2    cs      s^2 ]

Global System Assembly
All element stiffness matrices are assembled into the global stiffness matrix ‘K’ which satisfies:
K/U = F
U = global displacement vector  
F = global load vector  
Boundary conditions are applied by removing fixed degrees of freedom.
The reduced system for free DOFs is solved as K_f_f, U_f = F_f

--Program Workflow
Define Geometry
   - User inputs nodal coordinates and element connectivity.
Material and Section Properties
   - ‘A’ and ‘E’ are defined for each element.
Boundary Conditions & Loads
   - Supports and external nodal loads are specified.
Matrix Assembly  
   - Element stiffness matrices combined into global ‘K’.
   - Displacements solved via `numpy.linalg.solve’
Post-Processing
   - Axial forces and deformations calculated and plotted.

Example Input:

How many nodes? 3

Node 1 assumed at (0,0).

Enter x coordinate of node 2: 3

Enter y coordinate of node 2: 0

Enter x coordinate of node 3: 1.5

Enter y coordinate of node 3: 2.6

Enter start node for element 1: 1

Enter end node for element 1: 2

Enter cross-sectional area (m^2) for element 1: 0.01

Enter Young's Modulus (Pa) for element 1: 210e9

Example Output:

Undeformed and Deformed Shape plotted on Matplotlib 

Blue = Original, Red = Deformed  
