#/usr/bash/python
# a python script
import sys
sys.path.append(os.path.dirname(__file__)+"/..")
from auto_selfcal import auto_selfcal, split_calibrated_final
import glob
import casatasks
from casatools import msmetadata

os.system('cp ~/bin/casa/auto_selfcal/bin/ma_script.py .')
os.system('mv ma_script.py ma_script_run.txt')
msmda = msmetadata()

# Mac builds of CASA lack MPI and error without this try/except
try:
   from casampi.MPIEnvironment import MPIEnvironment
   parallel=MPIEnvironment.is_mpi_enabled
except:
   parallel=False
img_name='GMP3779_7_12_autoselfcal'
phase_c='J2000 12:59:02.5974 +027.38.25.374'
imsz=[1280,1152]
v_c= 5400 #5700 #km/s
deltaV=800
channel_wv=3.8
nchan=int(deltaV/channel_wv)
v0=int(v_c-deltaV/2.)
v1=int(v_c+deltaV/2.)
f_mid=230.538*(1.-v_c/3e5)
f0=230.538*(1.-v0/3e5) #GHz
f1=230.538*(1.-(v0+deltaV)/3e5) #GHz

self_list=[]
for fname in glob.glob('./*'): # change directory as needed
    if fname.endswith('.ms') and not '.selfcal.' in fname and not '_targets.ms' in fname:
        file=fname.replace('./','')
        self_list.append(file)
#print(self_list)
#self_list=['uid___A002_Xf96bbc_Xfdee.ms']

for file in self_list:

    #split_calibrated_final(file, overwrite=True)
    targets=str(file.replace('.ms','_targets.ms'))
    msmda.open(targets)
    spw_list=-1
    vmin=0
    vmax=0
    for spw_id in [0,1,2,3]:
        chan_freqs = msmda.chanfreqs(spw=spw_id, unit="Hz")/1e9

        if chan_freqs[-1]<f_mid and chan_freqs[0]>f_mid:
            spw_list=spw_id
            print(chan_freqs[0],f0,f_mid,f1,chan_freqs[-1],len(chan_freqs))
            vmin=3e5*(1.-chan_freqs[0]/230.538)
            vmax=3e5*(1.-chan_freqs[-1]/230.538)
    #   Always close the tool when finished
    msmda.done()
    print(vmin,v0,v1,vmax)
    nchan_pre=min(int(abs(v0-vmin)/channel_wv),25)
    nchan_post=min(int(abs(v1-vmax)/channel_wv),25)
    mstransform(targets,outputvis=targets.replace('_targets.ms','_targets_pre.ms'),spw=str(spw_list),datacolumn='data',regridms=True,nchan=nchan_pre,start=str(v0-nchan_pre*channel_wv)+'km/s',width=str(channel_wv)+'km/s',phasecenter=phase_c,restfreq='230.538GHz',mode='velocity',nspw=1)
    mstransform(targets,outputvis=targets.replace('_targets.ms','_targets_post.ms'),spw=str(spw_list),datacolumn='data',regridms=True,nchan=nchan_post,start=str(v1)+'km/s',width=str(channel_wv)+'km/s',phasecenter=phase_c,restfreq='230.538GHz',mode='velocity',nspw=1)
    mstransform(targets,outputvis=targets.replace('_targets.ms','_targets_mid.ms'),spw=str(spw_list),datacolumn='data',regridms=True,nchan=int(nchan),start=str(v0)+'km/s',width=str(channel_wv)+'km/s',phasecenter=phase_c,restfreq='230.538GHz',mode='velocity',nspw=1)
    concat([targets.replace('_targets.ms','_targets_pre.ms'),targets.replace('_targets.ms','_targets_mid.ms'),targets.replace('_targets.ms','_targets_post.ms')],concatvis=targets.replace('_targets.ms','_targets_vel.ms'))


    auto_selfcal(targets.replace('_targets.ms','_targets_vel.ms'), parallel=parallel,spectral_average=False,optimize_spw_combine=True,minsnr_to_proceed=2.0,allow_gain_interpolation=True,guess_scan_combine=True,allow_cocal=False,delta_beam_thresh=10,apply_to_target_ms=False,check_all_spws=False,apply_cal_mode_default='calflag',inf_EB_gaintype='G',inf_EB_gaincal_combine='scan')
    os.system('rm -r *.tt0')
    os.system('rm -r *.mask')
    os.system('rm -r '+targets.replace('_targets.ms','_targets_pre.ms'))
    os.system('rm -r '+targets.replace('_targets.ms','_targets_post.ms'))
    os.system('rm -r '+targets.replace('_targets.ms','_targets_mid.ms'))

    os.system('mv weblog weblog_'+file.replace('.ms',''))

#uvc_list=[]
#for fname in glob.glob('./*'): # change directory as needed
#    if fname.endswith('.selfcal.ms'):# and fname.startswith('Target'):
#        file=fname.replace('./','')
#        uvc_list.append(file)
#print(uvc_list)

#for vis in uvc_list:
    #uvcontsub(vis=vis,spw='',fitspec='0,2',fitorder=0,outputvis=vis.replace('.ms','_line.ms'),datacolumn='data') #CHECK SPW!!!
#    uvcontsub(vis=vis,spw='',fitspec='1:0~10;'+str(nchan-10)+'~'+str(nchan-1),fitorder=1,outputvis=vis.replace('.ms','_line.ms'),datacolumn='data') #CHECK SPW!!!



##imaging
img_list=[]
for fname in glob.glob('./*'): # change directory as needed
    #if fname.endswith('.selfcal_line.ms'):# and fname.startswith('Target'):
    if fname.endswith('.selfcal.ms'):
        file=fname.replace('./','')
        img_list.append(file)

print(img_list)
#vstart=v0
#for im in img_list:
#    msmda.open(im)
#    spw_list=-1
#    for spw_id in [0,1]:
#    chan_freqs = msmda.chanfreqs(spw=1, unit="Hz")/1e9
#    print(v0,vstart,3e5*(1.-chan_freqs[0]/230.538))
#    if 3e5*(1.-chan_freqs[0]/230.538)>vstart:
#        vstart=3e5*(1.-chan_freqs[0]/230.538)
#
#        print(3e5*(1.-chan_freqs[0]/230.538),3e5*(1.-chan_freqs[-1]/230.538))
#


#channel_wv=3.81
#nchan=int(deltaV/channel_wv)

tclean(vis=img_list,selectdata=True,field='',spw='1',timerange='',uvrange='',antenna='',scan='',observation='',intent='',datacolumn='data',imagename=img_name,imsize=imsz,cell='0.15arcsec',phasecenter=phase_c,stokes='I',projection='SIN',startmodel='',specmode='cube',reffreq='',outframe='',veltype='radio',restfreq='230.538GHz',interpolation='linear',perchanweightdensity=True,gridder='mosaic',facets=1,psfphasecenter='',wprojplanes=1,vptable='',mosweight=True,aterm=True,psterm=False,wbawp=True,conjbeams=False,cfcache='',usepointing=False,computepastep=360.0,rotatepastep=360.0,pointingoffsetsigdev=[],pblimit=0.2,normtype='flatnoise',deconvolver='multiscale',scales=[0, 6, 16],nterms=2,smallscalebias=0.0,fusedthreshold=0.0,largestscale=-1,restoration=True,restoringbeam='common',pbcor=False,outlierfile='',weighting='briggs',robust=0.5,npixels=0,uvtaper=[],niter=10000,gain=0.1,threshold='3.50mJy/beam',nsigma=0.0,cycleniter=100,cyclefactor=3.0,minpsffraction=0.05,maxpsffraction=0.8,interactive=False,nmajor=-1,fullsummary=False,usemask='auto-multithresh',mask='',pbmask=0.5,sidelobethreshold=2.0,noisethreshold=4.25,lownoisethreshold=1.5,negativethreshold=0.0,smoothfactor=1.0,minbeamfrac=0.3,cutthreshold=0.01,growiterations=75,dogrowprune=True,minpercentchange=-1.0,verbose=False,fastnoise=True,restart=True,savemodel='none',calcres=True,calcpsf=True,psfcutoff=0.35,parallel=True )
#nchan=nchan,start=str(int(vstart)+1)+'km/s',width=str(channel_wv)+'km/s',
#Primary beam correction
imcontsub(imagename=img_name+'.image/',linefile=img_name+'_line.image/',contfile=img_name+'_cont.image/',fitorder=0)
impbcor( imagename=img_name+'_line.image/',pbimage=img_name+'.pb/',outfile=img_name+'_line_pbcorr.image/',overwrite=False,box='',region='',chans='',stokes='I',mask='',mode='divide',cutoff=-1.0,stretch=False )
#Export
exportfits(imagename=img_name+'_line_pbcorr.image/',fitsimage=img_name+'_line_pbcorr.fits',velocity=True)
exportfits(imagename=img_name+'.pb/',fitsimage=img_name+'_pb.fits',velocity=True)

## SOFIA run
#
#os.system('sofia moments_par_file.par')
