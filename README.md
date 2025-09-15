



For near-cloud sounding analyses:

1) Hugh_env_profiles_oneprof_allsondes.m: generates sounding files for each updraft for passing into Morrison et al 2022's equations to predict cloud top height.

2) 





_____________________________________________________________
Final analyses:

1) LASSO_dBZcell_wID_mesoW_matchup_COMMONCORE_testbed100m_v13f.m: primary analysis script to concatinate sub-LFC mesoscale ascent objects, cloudy updraft objects, and near-cloud environmental data from all 7 cases.

2) All_sims_z300m_D4_examples.m:  generates low-level convergence plots for each case.

3) momentum_profiles.m : Makes x-z cross section plots of u-w vectors, qv & qcloud plots.

4) plot_terr_lassodomains.m: for making figure 1 - lasso domains.

5) GOESvisualizer.m: used to make the obs & LASSSo cloud top maps relative to terrain and the sub-LFC ascent object figure.


_____________________________________________________________
Data: 

1) profs8_theory_ztop_L04R_v3.dat: cloud top height results from Morrison et al. 2022 equations.
2) 
