SCRIPTS AND TEMPLATE FILES USED IN ALMA JELLY PROCESSING:

1) ma_script.py: auto-selfcal and Robust=0 imaging with Sofia

Save ma_script.py in auto_selfcal/bin

Then
1) Store every calibrated MS downloaded from NAS in one folder
2) Store sofia .par files in your home. Otherwise, specifiy the path in the code
3) Set imagename, phase_center, central V and deltaV, and imsize for cleaning in the autoselfcal_par.txt file and save it in the same folder with the MS files
4) From the folder, run the script with:
casa -c [path to auto_selfcal]/auto_selfcal/bin/ma_script.py

SCRIPT STEPS:
1) Split observations in targets
2) split targets in pre spw (25 channels)+main spw centered on central V + post spw (25 channels)
3) auto_selfcal [auto-multithresh masking parameters require fine tuning]
4) uvcontsub: continuum is fitted in spw 0,2
5) imaging spw 1 only (aka the line)
   1) dirty image
   2) measure expected noise threshold and identify correct channels from dirty image
   3) Deep imaging down to 10 sigma sigma with high-thresh mask
   4) SOFIA mask at 6 sigma+dilation on the final image
   6) Deeper cleaning down to 4 sigma withion SOFIA mask
   7) SOFIA mask at 4 sigma+reliability+dilation on the final image
   8) Deeper cleaning down to 1.5 sigma withion SOFIA mask

7) impbcor
8) exportfits
9) Moment0,1 and 2 with SOFIA at 3.5 sigma w/o reliability

   use at your own risk.

2) script_clean.py: Robust=2 imaging
   Save script in auto_selfcal/bin
   MUST be run after ma_script.py
   
