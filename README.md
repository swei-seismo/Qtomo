# Qtomo
Q tomography from t*
This package is still under development

Procedures:
1. t* measurements obtained from a seperate program
2. Create a grid of nodes for tomography: currently a 3D grid over a 1D grid (creatgrid.m)
3. Construct the G matrix (genGsameQpQs.bash calling tomosrc/Gatten3D1D)
4. (Optional) Plot ray tracing results (pl_raypath.m)
5. Tomography inversion

5.1. Invert for Q with a large inversin problem, use iterative LSQR (invertQp.m)

5.2. Invert for Q with a small inversin problem, use non-negative least square (invertsameQpQsnnls.m)

5.3. Invert for Qp/Qs with a small inversion problem, use SVD (invertsameQps3D1Dsvd_vard.m)
6. (Optional) Calculate Qk based on Qp/Qs (calQk_Qps.m)
7. Plot tomography results (plot*.bash calling atten4gmt.py and hits4gmt.py)
8. Create a grid on nodes for checkerboard tests (creat_checker.m)
9. Compute synthetic t* (runsynchecksameQpQs.bash calling tomosrc/syn_tstar and addnoise.py)
10. Repeat 3-7 for the checkerboard tests

tomosrc includes all Fortran codes. Please read 'Makefile' for descriptions. Note most programs require ifort (Intel Fortran Compiler).
