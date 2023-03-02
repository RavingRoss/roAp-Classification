<p align="center">
  <img width="600" height="200" src="[https://www.python.org/python-.png](https://github.com/RavingRoss/roAp-Classification/blob/main/Sphar_l3_m0_RB.gif)">
</p>
# roAp-Classification
The goal of this project is to find rapidly oscillating Ap stars (roAps) candidates in star clusters using known data. 
## Description
The classification of roAps is done by using data from the TESS observations, compiled/analyzed in the java based runtime environment TOPCAT from the gaia archives, of known roAps (~ 190) to find roAp candidates. These candidates are picked by analyzing clusters in TOPCAT and choosing which are the best options. For this project we have data from NGC 2264, Christmas Tree cluster, and will be normalizing the dataset to the known roAps dataset. The plot we will be using to do this process is known as a Hertzsprung-Russell diagram (HRD), with axes 'absolute magnitude in G filter' vs. 'color index of B-R' (essentially Luminosity vs. Temperature). The reason for using the HRD to analyze the plots using regression techniques and not using pure data techniques is so we do not have to worry about error analysis. The HRD is efficient in the sense there is little to no error, the only error being dust extinction (which is accounted for in the datasets using Starhorse archives through TOPCAT) and syhstematic error from the gaia archives which is assumed to be negligible. 

The length of the project is only a few months, though we will keep it updated due to the need of finding more candidates in other clusters, not only NGC 2264.
