Example for Qp tomography

Inputs:
Q.grid         3D grid of Q
Q.prem         PREM Q model
Q1D.grid       1D grid of Q beneath Q.grid
event.lst      Event list: Longitude, Latitude, Depth (km)
lau.vel        3D grid of V, has to fully include Q.grid
p_tstar.dat    P-wave t* measurements
               station_name event_latitude event_longitude event_depth(km) t*(P) t*(P)_error not_used average_1/Qp
s_tstar.dat    S-wave t* measurements (not necessary for Qp tomography)
               station_name event_latitude event_longitude event_depth(km) t*(S) t*(S)_error t*(P) t*(P)_error average_1/Qs average_Qp/Qs
station.lst    List of station
               First line: number of stations
               Following lines: station_name longitude latitude elevation(km) not_used not_used

Important terminal outputs:
skip           Number of t* measurements skipped by various reasons, should be small
nbigmis10      Number of bad rays whose measured traveltime is 10-s different from the traveltime from ray tracing
nbigmis5       Number of bad rays whose measured traveltime is 5-s different from the traveltime from ray tracing

Outputs:
Qinv_*.p          Main output file of 1/Q, use atten4gmt.py to plot
Qmisfit.all       Overall inversion misfit for different sets of parameters
cartesian.Qgrid   Q grid converted to cartesian
cartesian.Vgrid   V grid converted to cartesian
cartesian.sta     Station locations converted to cartesian
hitsP             Hit counts at Q.grid nodes
hitsP1D           Hit counts at Q1D.grid nodes
raypaths.p        Raypaths of P waves, use pl_raypath.m to plot
raypaths.bad      Raypaths of bad waves that are counted in nbigmis10

Subdirectories (all files following the same format):
checker   Checkerboard test shown in Figs. 7a and 7b of Wei & Wiens (2018, EPSL)
synmod    Synthetic test shown in Figs. 7c and 7d of Wei & Wiens (2018, EPSL)
