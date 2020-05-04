#FCMP=gfortran
FCMP=ifort -fpe0
# module load intel/2015a

all: Gatten3D1D atten4gmt hits4gmt syn_tstar syn_ratio momcalc

# Main program for 3D1D attenuation inversion
Gatten3D1D: Gatten3D1D.f tomosubs.f raybend3D.f getgrowQ3D1D.f smoothg.f layertt.f
	$(FCMP) -mcmodel medium Gatten3D1D.f tomosubs.f raybend3D.f getgrowQ3D1D.f smoothg.f layertt.f -o Gatten3D1D

atten4gmt: atten4gmt.f
	$(FCMP) atten4gmt.f -o atten4gmt

hits4gmt: hits4gmt.f
	$(FCMP) hits4gmt.f -o hits4gmt

# Compute synthetic t*
syn_tstar: syn_tstar.f sumtstar.f tomosubs.f raybend3D.f mt19937ar.f layertt.f
	$(FCMP) syn_tstar.f sumtstar.f tomosubs.f raybend3D.f mt19937ar.f layertt.f -o syn_tstar

# Compute synthetic t*(S)/t*(P)
syn_ratio: syn_ratio.f sumtstar.f tomosubs.f raybend3D.f mt19937ar.f layertt.f
	$(FCMP) syn_ratio.f sumtstar.f tomosubs.f raybend3D.f mt19937ar.f layertt.f -o syn_ratio
