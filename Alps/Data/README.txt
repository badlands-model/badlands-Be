Data File for Southern French Alps and Mediterranean Sea models (C. Petit et al., 2022).

Here is the description of the data file you need to run the models:
:
eetHR.csv: grid of spatially variable Effective Elastic Thickness based on yield stress enveloppes (after Burov and Diament, 1995). Otherwise a constant EET value of 20 km is used.
epica_temp.csv: global temperatures from the EPICA-Dome C database
ice_conv.csv: ice ratio (1 is for full ice, 0 for no ice) between -100 ka and now (in model time, between 0 and +100 ka).
ice_conv200.csv: ice ratio (1 is for full ice, 0 for no ice) between -200 ka and now (in model time, between 0 and +200 ka).
paca_ice2.csv: maximum ice thickness map (in m)
pacaqtz.csv: quartz content map
pacasmooth5.csv: initial topography and bathymetry map
sealevel100ka.txt: sea level variations between -100 ka and now (in model time, between 0 and +100 ka).
sealevel200ka.txt: sea level variations between -200 ka and now (in model time, between 0 and +200 ka).

Here are the files that you need to compare with model outputs:

river_name.dat : xy file with age (x, in ka) and height above riverbed (y, in m)
Var_Pts.txt and Var_Pts5.txt: 2 lists of point IDs where to extract 10Be in submarine sediments
River_name_16_x.txt: 3 column files (extracted from selection in paraview) where the 3rd column is the pointID corresponding to the
upstream catchment of the given site (see sites in Mariotti et al., 2019). The 2 first columns are not used.