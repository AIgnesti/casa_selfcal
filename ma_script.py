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

with open('autoselfcal_par.txt') as f:
    lines = [line.rstrip('\n') for line in f if not line.startswith('#')]


img_name=lines[0]#'D100_7_12_autoselfcal'
phase_c=lines[1]#'J2000 12:59:02.5974 +027.38.25.374'
imsz=[int(lines[2]),int(lines[3])]
v_c= int(lines[4]) #5700 #km/s
deltaV=int(lines[5])

print(img_name,phase_c,imsz,v_c,deltaV)


channel_wv=3.8
nchan=int(deltaV/channel_wv)
v0=int(v_c-deltaV/2.)
v1=int(v_c+deltaV/2.)
f_mid=230.538*(1.-v_c/3e5)
f0=230.538*(1.-v0/3e5) #GHz
f1=230.538*(1.-(v0+deltaV)/3e5) #GHz

# Identify MS

self_list=[]
for fname in glob.glob('./*'): # change directory as needed
    if fname.endswith('.ms') and not '.selfcal.' in fname and not 'targets' in fname:
        file=fname.replace('./','')
        self_list.append(file)


# AUTO-SELFCAL per MS: split into targets -> split into line -> autoselfcal

for file in self_list:

    split_calibrated_final(file, overwrite=True)
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
    msmda.done()
    print(vmin,v0,v1,vmax)
    nchan_pre=min(int(abs(v0-vmin)/channel_wv),20)
    nchan_post=min(int(abs(v1-vmax)/channel_wv),20)
    mstransform(targets,outputvis=targets.replace('_targets.ms','_targets_pre.ms'),spw=str(spw_list),datacolumn='data',regridms=True,nchan=nchan_pre,start=str(v0-nchan_pre*channel_wv)+'km/s',width=str(channel_wv)+'km/s',phasecenter=phase_c,restfreq='230.538GHz',mode='velocity',nspw=1)
    mstransform(targets,outputvis=targets.replace('_targets.ms','_targets_post.ms'),spw=str(spw_list),datacolumn='data',regridms=True,nchan=nchan_post,start=str(v1)+'km/s',width=str(channel_wv)+'km/s',phasecenter=phase_c,restfreq='230.538GHz',mode='velocity',nspw=1)
    mstransform(targets,outputvis=targets.replace('_targets.ms','_targets_mid.ms'),spw=str(spw_list),datacolumn='data',regridms=True,nchan=int(nchan),start=str(v0)+'km/s',width=str(channel_wv)+'km/s',phasecenter=phase_c,restfreq='230.538GHz',mode='velocity',nspw=1)
    concat([targets.replace('_targets.ms','_targets_pre.ms'),targets.replace('_targets.ms','_targets_mid.ms'),targets.replace('_targets.ms','_targets_post.ms')],concatvis=targets.replace('_targets.ms','_targets_vel.ms'))


    auto_selfcal(targets.replace('_targets.ms','_targets_vel.ms'), parallel=parallel,spectral_average=False,optimize_spw_combine=True,minsnr_to_proceed=3.0,gaincal_minsnr=3.0,allow_gain_interpolation=True,guess_scan_combine=True,allow_cocal=False,delta_beam_thresh=10,apply_to_target_ms=False,check_all_spws=False,apply_cal_mode_default='calflag',inf_EB_gaintype='T',inf_EB_gaincal_combine='scan')
    os.system('rm -r *.tt0')
    os.system('rm -r *.mask')
    os.system('rm -r '+targets.replace('_targets.ms','_targets_pre.ms'))
    os.system('rm -r '+targets.replace('_targets.ms','_targets_post.ms'))
    os.system('rm -r '+targets.replace('_targets.ms','_targets_mid.ms'))

    os.system('mv weblog weblog_'+file.replace('.ms',''))


# Continuum subtraction

uvc_list=[]
for fname in glob.glob('./*'): # change directory as needed
    if fname.endswith('.selfcal.ms'):# and fname.startswith('Target'):
        file=fname.replace('./','')
        uvc_list.append(file)
print(uvc_list)

for vis in uvc_list:
    #uvcontsub(vis=vis,spw='',fitspec='0,2',fitorder=0,outputvis=vis.replace('.ms','_line.ms'),datacolumn='data') #CHECK SPW!!!
    uvcontsub(vis=vis,spw='',fitspec='0,2',fitorder=1,outputvis=vis.replace('.ms','_line.ms'),datacolumn='data') #CHECK SPW!!!


##imaging
img_list=[]
for fname in glob.glob('./*'): # change directory as needed
    #if fname.endswith('.selfcal_line.ms'):# and fname.startswith('Target'):
    if fname.endswith('.selfcal_line.ms'):
        file=fname.replace('./','')
        img_list.append(file)

print(img_list)
nchan_img=10000
vimg=1000000
for im in img_list:
    msmda.open(im)
    spw_list=-1
    chan_freqs = msmda.chanfreqs(spw=1, unit="Hz")/1e9
    chan_w = msmda.chanwidths(spw=1, unit="Hz")/1e9
    print(im,chan_freqs[0],3e5*(1.-chan_freqs[-1]/230.538),len(chan_freqs))
    if len(chan_freqs)<=nchan_img:
        nchan_img=len(chan_freqs)
    if 3e5*(1.-chan_freqs[-1]/230.538)<=vimg:
        vimg=3e5*(1.-chan_freqs[-1]/230.538)
#
    msmda.done()


### Dirty image to compute goal noise
def psf_per_channel(cube_path):

    ia.open(cube_path)

    header = ia.summary()
    ia.close()
    chan_list=[]
    if 'perplanebeams' in header:
        beams_data = header['perplanebeams']['beams']


        for chan_idx in beams_data.keys():
            chan_num = int(chan_idx.replace('*', ''))

            beam = beams_data[chan_idx]['*0']

            major = beam['major']['value']
            minor = beam['minor']['value']
            if major<2. and minor<2.:
                chan_list.append(chan_num)
            #print(f"{chan_num:<10}{major:<15.4f}{minor:<15.4f}{pa:<10.2f}")
        return min(chan_list),max(chan_list)

tclean(vis=img_list,selectdata=True,field='',spw='1',timerange='',uvrange='',antenna='',scan='',observation='',intent='',datacolumn='data',imagename=img_name+'_dirty',imsize=imsz,cell='0.15arcsec',phasecenter=phase_c,stokes='I',projection='SIN',startmodel='',specmode='cube',reffreq='',outframe='',veltype='radio',restfreq='230.538GHz',interpolation='linear',perchanweightdensity=True,gridder='mosaic',facets=1,psfphasecenter='',wprojplanes=1,vptable='',mosweight=True,aterm=True,psterm=False,wbawp=True,conjbeams=False,cfcache='',usepointing=False,computepastep=360.0,rotatepastep=360.0,pointingoffsetsigdev=[],pblimit=0.2,normtype='flatnoise',deconvolver='multiscale',scales=[0, 6, 12],nterms=2,smallscalebias=0.0,fusedthreshold=0.0,largestscale=-1,restoration=True,restoringbeam='',pbcor=False,outlierfile='',weighting='briggs',robust=0.5,npixels=0,uvtaper=[],niter=0,gain=0.1,threshold='2.0mJy/beam',nsigma=0.0,cycleniter=100,cyclefactor=3.0,minpsffraction=0.05,maxpsffraction=0.8,interactive=False,nmajor=-1,fullsummary=False,usemask='auto-multithresh',mask='',pbmask=0.2,sidelobethreshold=2.0,noisethreshold=4.25,lownoisethreshold=1.5,negativethreshold=0.0,smoothfactor=1.0,minbeamfrac=0.3,cutthreshold=0.01,growiterations=75,dogrowprune=True,minpercentchange=-1.0,verbose=False,fastnoise=True,restart=True,savemodel='none',calcres=True,calcpsf=True,psfcutoff=0.35,parallel=True )


# Identify 7+12 channels
min_chan,max_chan=psf_per_channel(img_name+'_dirty.psf')


# Estimate noise
chanstat=imstat(imagename=img_name+'_dirty.image',chans=str(min_chan)+'~'+str(min_chan+4))
rms1= chanstat['rms'][0]
chanstat=imstat(imagename=img_name+'_dirty.image',chans=str(max_chan-5)+'~'+str(max_chan-1))
rms2= chanstat['rms'][0]
rms=0.5*(rms1+rms2)*1e3
print('GOAL RMS: ',rms,min_chan,max_chan)

#### Deep cleaning

tclean(vis=img_list,selectdata=True,field='',spw='1',timerange='',uvrange='',antenna='',scan='',observation='',intent='',datacolumn='data',imagename=img_name,imsize=imsz,cell='0.15arcsec',start=min_chan,nchan=int(max_chan-min_chan),phasecenter=phase_c,stokes='I',projection='SIN',startmodel='',specmode='cube',reffreq='',outframe='',veltype='radio',restfreq='230.538GHz',interpolation='linear',perchanweightdensity=True,gridder='mosaic',facets=1,psfphasecenter='',wprojplanes=1,vptable='',mosweight=True,aterm=True,psterm=False,wbawp=True,conjbeams=False,cfcache='',usepointing=False,computepastep=360.0,rotatepastep=360.0,pointingoffsetsigdev=[],pblimit=0.2,normtype='flatnoise',deconvolver='multiscale',scales=[0, 6, 12],nterms=2,smallscalebias=0.0,fusedthreshold=0.0,largestscale=-1,restoration=True,restoringbeam='common',pbcor=False,outlierfile='',weighting='briggs',robust=0.5,npixels=0,uvtaper=[],niter=500000,gain=0.1,threshold=str(rms)+'mJy/beam',nsigma=0.0,cycleniter=100,cyclefactor=3.0,minpsffraction=0.05,maxpsffraction=0.8,interactive=False,nmajor=-1,fullsummary=False,usemask='auto-multithresh',mask='',pbmask=0.2,sidelobethreshold=2.0,noisethreshold=4.25,lownoisethreshold=1.5,negativethreshold=0.0,smoothfactor=1.0,minbeamfrac=0.3,cutthreshold=0.01,growiterations=75,dogrowprune=True,minpercentchange=-1.0,verbose=False,fastnoise=True,restart=True,savemodel='none',calcres=True,calcpsf=True,psfcutoff=0.35,parallel=True )
#Primary beam correction
#imcontsub(imagename=img_name+'.image/',linefile=img_name+'_line.image/',contfile=img_name+'_cont.image/',fitorder=0)
impbcor( imagename=img_name+'.image/',pbimage=img_name+'.pb/',outfile=img_name+'_pbcorr.image/',overwrite=False,box='',region='',chans='',stokes='I',mask='',mode='divide',cutoff=-1.0,stretch=False )
#Export
exportfits(imagename=img_name+'_pbcorr.image/',fitsimage=img_name+'_pbcorr.fits',velocity=True)
exportfits(imagename=img_name+'.pb/',fitsimage=img_name+'_pb.fits',velocity=True)
