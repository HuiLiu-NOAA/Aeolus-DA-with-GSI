Brief description of the codes with Aeolus DA in NOAA/NCEP Global_workflow v15.3:

These Aeolus DA codes are based on the original codes developed by Will Mccarty at NASA. The codes are modifdied to accormodate ESA Aeolus space-lidar winds in the NOAA GSI system. 

Specifically, the Aeolus winds are assimilted as wind observations as radiosonde winds with only one component, i.e., along the Horizontal-Line-Of-Sight (HLOS). The observation error estimated for the Aeolus winds for use in GSI are derived using the OmB and the Hollingsworth and Lonnberg method. The Aeolus winds are assimilated at height level. The Mie winds are thinned (not superobed) to 90km resolution in setupdw.f90. An additional Total-Least-Squres (TLS) bias correction is developed and applied to the Aeolus innovations to remove impact of the biases in either the FV3GFS background or Aeolus winds. The TLS bias corrections are calculated based on the Aeolus winds over previous 24 hours, which is set in tls_bc.f90. 
 

Contributors: Hui Liu, Ross N. Hoffman, Kayo Ide, Kevin Garrett, Katherine Lukens

Archive locations:    

HPSS:  /NCEPDEV/nesdis-drt/5year/Hui.Liu/HERA/scratch/globalaeo19r_onlineBC_archive20230911.tar,  

https://github.com/HuiLiu-NOAA/Aeolus-DA-with-GSI


Aeolus related codes and data:

read_lidar.f90:   
  read in Aeolus winds.

setupdw.f90:  
  process QCs of Aeolus winds, thinning of Mie winds, and calculate OmB of Aeolus winds.

 tls_bc.f90 (in TLSonlineBC.tar):   
  a small code to calculate the additional TLS bias correction over the previous 24-hours.

DATA_TLS_BC.tar:   
  Aeolus TLS bias correction coefficients used in the TLS bias correction applied to the OmB in setupdw.f90. 

exglobal_analysis_fv3gfs.sh.ecf:   
  main GSI analysis script modified to read Aeolus data and call the TLS bias correction.

setup2021b.csh:   
  the script to submit Aeolus OSE.


----- References and how to cite these codes ----
1.	Garrett K, Hui Liu, K. Ide, R.N. Hoffman, and K. Lukens, 2022: Optimization and Impact Assessment of Aeolus HLOS Wind Data Assimilation in NOAA’s Global Forecast System, Quart. J. of Royal Meteor. Soc., v148, p.2703-2716, doi: 10.1002/qj.4331
2.	Liu Hui, K. Garrett, K. Ide, and R.N. Hoffman, 2023: On the Use of Consistent Bias Corrections to Enhance the Impact of Aeolus Level-2B Rayleigh Winds on NOAA Global Forecast Skill, Quart. J. of Royal Meteor. Soc. (accepted).
3.	Liu Hui, K. Garrett, K. Ide, R.N. Hoffman, and K. Lukens, 2022: A Statistically Optimal Analysis of Systematic Differences between Aeolus HLOS Winds and NOAA’s Global Forecast System, Atmos. Meas. Tech., 15, 3925-3940, 2022. doi:10.5194/amt-15-3925-2022.
