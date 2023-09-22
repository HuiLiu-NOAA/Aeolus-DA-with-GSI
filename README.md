Brief description of the Global_workflow with Aeolus DA:

Archive locations:
HPSS:  /NCEPDEV/nesdis-drt/5year/Hui.Liu/HERA/scratch/globalaeo19r_onlineBC_archive20230911.tar
S4: /data/users/huiliu/globalaeo19r_onlineBC_archive20230911.tar


Aeolus related codes and data:

read_lidar.f90 and setupdw.f90:
  main codes to read in and process Aeolus winds

TLSonlineBC/
  a small code and script to do bias correction

DATA_TLS_BC/
  Aeolus bias correction coefficients

exglobal_analysis_fv3gfs.sh.ecf:
  main analysis script to read Aeolus data and call bias correction
