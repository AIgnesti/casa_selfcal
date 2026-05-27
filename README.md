save ma_script.py in autoselfcal/bin
Set imagename, phase_center, central V and deltaV, and imsize for cleaning in the autoselfcal_par.txt file in the same folder as the .ms
SCRIPT STEPS:
1) Split observations in targets
2) split targets in pre spw (25 channels)+main spw centered on central V + post spw (25 channels)
3) auto_selfcal
4) imaging
5) imcontsub
6) impbcor
7) exportfits

   use at your own risk.
