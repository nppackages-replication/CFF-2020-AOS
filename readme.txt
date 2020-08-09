
DATE: 19-MAY-2019
AUTHORS: Matias D. Cattaneo, Max H. Farrell and Yingjie Feng

This folder contains all files needed to replicate the simulation study in SA:

1. preamble_rep.R: supporting functions

2. simul_pointwise_rep.R: generate pointwise results

3. simul_uniform_rep.R: generate uniform results

4. model.csv: 19 specifications used in the simulation

Note: 
1. Please install R package "lspartition" first. 
2. To replicate the simulation study, make sure all directories used in these files exist. You may change them if you want.
3. Uniform inference results may be slightly different from that in the SA, since the simulation of SA relies on a modified version of lspartition in order to speed up computation and random draws used in uniform inference are slightly different. 