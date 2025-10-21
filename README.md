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
Enter x coordinate of node 3: 3
Enter y coordinate of node 3: 2

Node Coordinate Library:
Node 1: (0.0, 0.0)
Node 2: (3.0, 0.0)
Node 3: (3.0, 2.0)
How many elements? 3
Enter start node for element 1: 1
Enter end node for element 1: 2
Enter cross-sectional area (m^2) for element 1: 0.001
Enter Young's Modulus (Pa) for element 1: 200e9
Enter start node for element 2: 1
Enter end node for element 2: 3
Enter cross-sectional area (m^2) for element 2: 0.001
Enter Young's Modulus (Pa) for element 2: 200e9
Enter start node for element 3: 2
Enter end node for element 3: 3
Enter cross-sectional area (m^2) for element 3: 0.001
Enter Young's Modulus (Pa) for element 3: 200e9
Enter supports as: node direction(s)
Example: 1 x,y  or  2 y  or  3 x
Use commas between multiple directions.

How many supported nodes? 2
Enter support (e.g. '1 x,y' or '2 y'): 1 x,y 2 y
Enter support (e.g. '1 x,y' or '2 y'): 2 y

Support summary:
Node 1: fixed in x, y direction(s)
Node 2: fixed in y direction(s)

How many loaded nodes? 1
Enter node number where load is applied: 3
Enter Fx (N) at node 3: 0
Enter Fy (N) at node 3: -5000

Load summary:
Node 3: Fx = 0.0 N, Fy = -5000.0 N

Nodal Displacements:
Node 1: ux = 0.0000 mm, uy = 0.0000 mm
Node 2: ux = 0.0000 mm, uy = 0.0000 mm
Node 3: ux = 0.0333 mm, uy = -0.0500 mm

Element Axial Forces:
Element 1: Axial force = 0.00 N
Element 2: Axial force = 0.00 N
Element 3: Axial force = -5000.00 N

Example Output:
Undeformed and Deformed Shape plotted on Matplotlib 
Blue = Original, Red = Deformed  

Example Output:

Undeformed and Deformed Shape plotted on Matplotlib 

Blue = Original, Red = Deformed  
