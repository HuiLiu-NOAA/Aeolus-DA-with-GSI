#!/bin/csh
set echo

cd ush/rocoto

 set idate0 = 2020081000
 set edate0 = 2020111800

# ------ experiment name -------
 set pslot0 = aeo20rBC

# ------- set up run directory, also in config.base --------------------
 set conf0 = /scratch2/NESDIS/nesdis-rdo1/Hui.Liu/CODE/globalaeo19r_onlineBC/parm/config
 set  exp0 = /scratch2/NESDIS/nesdis-rdo1/Hui.Liu/ccc/work
 set  com0 = /scratch2/NESDIS/nesdis-rdo1/Hui.Liu/ccc/out
 set  diag = /scratch2/NESDIS/nesdis-rdo1/Hui.Liu/ccc/diagaeo20rBC
   mkdir $diag

 set   ics0 = /scratch2/NESDIS/nesdis-rdo1/Hui.Liu/ICs

 ./setup_expt.py --pslot $pslot0 --configdir $conf0  --idate $idate0  --edate $edate0 --icsdir $ics0  --comrot $com0  --expdir $exp0 --resdet 384  --resens 192 --nens 80 --gfs_cyc 1

 ./setup_workflow.py --expdir $exp0/$pslot0

 cd  $exp0/$pslot0

#  rocotorun -d $pslot0.db  -w $pslot0.xml

  crontab  $pslot0.crontab

