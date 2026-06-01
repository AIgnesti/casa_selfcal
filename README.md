Save ma_script.py in auto_selfcal/bin

Then
1) Store every calibrated MS downloaded from NAS in one folder
2) Set imagename, phase_center, central V and deltaV, and imsize for cleaning in the autoselfcal_par.txt file and save it in the same folder with the MS files
3) From the folder, run the script with:
casa -c [path to auto_selfcal]/auto_selfcal/bin/ma_script.py

SCRIPT STEPS:
1) Split observations in targets
2) split targets in pre spw (25 channels)+main spw centered on central V + post spw (25 channels)
3) auto_selfcal
4) uvcontsub
5) imaging
  a) dirty image
   b) measure expected noise threshold and identify correct channels from dirty image
   c) Deep imaging
7) impbcor
8) exportfits

   use at your own risk.
