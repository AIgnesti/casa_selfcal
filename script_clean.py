import sys
sys.path.append(os.path.dirname(__file__)+"/..")
from auto_selfcal import auto_selfcal, split_calibrated_final
import glob
import casatasks
from casatools import msmetadata
path_sofia_templates='~/'
msmda = msmetadata()
os.system('ulimit -n 5000')
# Mac builds of CASA lack MPI and error without this try/except
try:
   from casampi.MPIEnvironment import MPIEnvironment
   parallel=MPIEnvironment.is_mpi_enabled
except:
   parallel=False



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
            if major<2. and minor<2. and major>0. and minor>0.:
                chan_list.append(chan_num)
            #print(f"{chan_num:<10}{major:<15.4f}{minor:<15.4f}{pa:<10.2f}")
        return min(chan_list),max(chan_list)

def make_clean_mask(img_name,mask_thr,outfile,rel):
    os.system('rm -rf final_mask.mask')
    os.system('rm -rf final_mask_I.mask')
    impbcor(imagename=img_name+'.image/',pbimage=img_name+'.pb/',outfile=img_name+'_pbcorr.image/',box='',region='',chans='',stokes='I',mask='',mode='divide',cutoff=0.5,stretch=False,overwrite=True)
    exportfits(imagename=img_name+'_pbcorr.image/',fitsimage=img_name+'_pbcorr.fits',velocity=True,overwrite=True)
    exportfits(imagename=img_name+'.pb/',fitsimage=img_name+'_pb.fits',velocity=True,overwrite=True)

    os.system('sofia '+path_sofia_templates+'clean_mask_template.par input.data='+img_name+'_pbcorr.fits input.primaryBeam='+img_name+'_pb.fits scfind.threshold='+str(mask_thr)+' reliability.enable='+rel)

    importfits(imagename='final_mask.mask',fitsimage='final_clean_mask.fits',overwrite=True)

    # Some dark magic to convince casa to use the new mask
    ia.open('final_mask.mask')
    ia.adddegaxes(outfile='final_mask_I.mask', stokes='I', overwrite=True)
    ia.close()
    #immath(imagename='final_mask_I.mask', mode='evalexpr', outfile='final_mask_I_flat.mask', expr='IM0[IM0>0.]/IM0[IM0>0.]')#, mask='final_mask_I.mask>0.')

    imtrans(imagename='final_mask_I.mask', outfile=outfile, order=['Right Ascension', 'Declination', 'Stokes', 'Frequency'])

    os.system('rm -rf '+img_name+'.mask')

with open('autoselfcal_par.txt') as f:
    lines = [line.rstrip('\n') for line in f if not line.startswith('#')]

suffix='_Robust2_dilation5'
img_name=lines[0]+suffix#'D100_7_12_autoselfcal'
phase_c=lines[1]#'J2000 12:59:02.5974 +027.38.25.374'
imsz=[int(lines[2]),int(lines[3])]
v_c= int(lines[4]) #5700 #km/s
deltaV=int(lines[5])

print(img_name,phase_c,imsz,v_c,deltaV)


channel_wv=3.8 #km
nchan=int(deltaV/channel_wv)
v0=int(v_c-deltaV/2.)
v1=int(v_c+deltaV/2.)
f_mid=230.538*(1.-v_c/3e5)
f0=230.538*(1.-v0/3e5) #GHz
f1=230.538*(1.-(v0+deltaV)/3e5) #GHz


img_list=[]
for fname in os.listdir(): # change directory as needed

    if fname.endswith('.selfcal_line.ms'):#_line
        img_list.append(fname)

print(img_list)

##Dirty image
if not img_name+'_dirty.psf' in os.listdir():
    tclean(vis=img_list,selectdata=True,field='',spw='1',timerange='',uvrange='',phasecenter=phase_c,antenna='',scan='',observation='',intent='',datacolumn='corrected',imagename=img_name+'_dirty',imsize=imsz,cell='0.15arcsec',stokes='I',projection='SIN',startmodel='',specmode='cube',reffreq='',outframe='',veltype='radio',restfreq='230.538GHz',interpolation='linear',perchanweightdensity=True,gridder='mosaic',facets=1,psfphasecenter='',wprojplanes=1,vptable='',mosweight=True,aterm=True,psterm=False,wbawp=True,conjbeams=False,cfcache='',usepointing=False,computepastep=360.0,rotatepastep=360.0,pointingoffsetsigdev=[],pblimit=0.2,normtype='flatnoise',deconvolver='multiscale',scales=[0, 6, 12],nterms=2,smallscalebias=0.0,fusedthreshold=0.0,largestscale=-1,restoration=True,restoringbeam='',pbcor=False,outlierfile='',weighting='briggs',robust=2.0,npixels=0,uvtaper=[],niter=0,gain=0.2,threshold='2.0mJy/beam',nsigma=0.0,cycleniter=100,cyclefactor=3.0,minpsffraction=0.05,maxpsffraction=0.8,interactive=False,nmajor=-1,fullsummary=False,usemask='auto-multithresh',mask='',pbmask=0.2,sidelobethreshold=2.0,noisethreshold=4.25,lownoisethreshold=1.5,negativethreshold=0.0,smoothfactor=1.0,minbeamfrac=0.3,cutthreshold=0.01,growiterations=75,dogrowprune=True,minpercentchange=-1.0,verbose=False,fastnoise=True,restart=True,savemodel='none',calcres=True,calcpsf=True,psfcutoff=0.35,parallel=True )


# Identify 7+12 channels
min_chan,max_chan=psf_per_channel(img_name+'_dirty.psf')


# Noise estimate
chanstat=imstat(imagename=img_name+'_dirty.image',chans=str(min_chan)+'~'+str(min_chan+10))
rms1= chanstat['rms'][0]
chanstat=imstat(imagename=img_name+'_dirty.image',chans=str(max_chan-11)+'~'+str(max_chan-1))
rms2= chanstat['rms'][0]
rms=(0.5*(rms1+rms2)*1e3) # Clean threshold at 1 sigma in mJy
print('GOAL RMS: ',rms,min_chan,max_chan)

#### Deep cleaning
# down to 10 sigma with high-thresh mask
tclean(vis=img_list,selectdata=True,field='',spw='1',timerange='',uvrange='',antenna='',scan='',observation='',intent='',datacolumn='corrected',imagename=img_name,imsize=imsz,cell='0.15arcsec',start=min_chan,nchan=int(max_chan-min_chan),phasecenter=phase_c,stokes='I',projection='SIN',startmodel='',specmode='cube',reffreq='',outframe='',veltype='radio',restfreq='230.538GHz',interpolation='linear',perchanweightdensity=True,gridder='mosaic',facets=1,psfphasecenter='',wprojplanes=1,vptable='',mosweight=True,aterm=True,psterm=False,wbawp=True,conjbeams=False,cfcache='',usepointing=False,computepastep=360.0,rotatepastep=360.0,pointingoffsetsigdev=[],pblimit=0.2,normtype='flatnoise',deconvolver='multiscale',scales=[0, 6, 12],nterms=2,smallscalebias=0.6,fusedthreshold=0.0,largestscale=-1,restoration=True,restoringbeam='common',pbcor=False,outlierfile='',weighting='briggs',robust=2.0,npixels=0,uvtaper=[],niter=500000,gain=0.2,threshold=str(round(10.*rms,3))+'mJy/beam',nsigma=0.0,cycleniter=50,cyclefactor=3.0,minpsffraction=0.05,maxpsffraction=0.8,interactive=False,nmajor=-1,fullsummary=False,usemask='auto-multithresh',mask='',pbmask=0.2,sidelobethreshold=3.5,noisethreshold=4.25,lownoisethreshold=2.5,negativethreshold=0.0,smoothfactor=1.0,minbeamfrac=0.3,cutthreshold=0.01,growiterations=75,dogrowprune=True,minpercentchange=-1.0,verbose=False,fastnoise=False,restart=False,savemodel='none',calcres=True,calcpsf=True,psfcutoff=0.35,parallel=True )

# Making final clean mask

make_clean_mask(img_name,9.0,'final_mask_I_sorted'+suffix+'.mask','false')


tclean(vis=img_list,selectdata=True,field='',spw='1',timerange='',uvrange='',antenna='',scan='',observation='',intent='',datacolumn='corrected',imagename=img_name,imsize=imsz,cell='0.15arcsec',start=min_chan,nchan=int(max_chan-min_chan),phasecenter=phase_c,stokes='I',projection='SIN',startmodel='',specmode='cube',reffreq='',outframe='',veltype='radio',restfreq='230.538GHz',interpolation='linear',perchanweightdensity=True,gridder='mosaic',facets=1,psfphasecenter='',wprojplanes=1,vptable='',mosweight=True,aterm=True,psterm=False,wbawp=True,conjbeams=False,cfcache='',usepointing=False,computepastep=360.0,rotatepastep=360.0,pointingoffsetsigdev=[],pblimit=0.2,normtype='flatnoise',deconvolver='multiscale',scales=[0, 6, 12],nterms=2,smallscalebias=0.6,fusedthreshold=0.0,largestscale=-1,restoration=True,restoringbeam='common',pbcor=False,outlierfile='',weighting='briggs',robust=2.0,npixels=0,uvtaper=[],niter=500000,gain=0.2,threshold=str(round(3.5*rms,4))+'mJy/beam',nsigma=0.0,cycleniter=50,cyclefactor=2.0,minpsffraction=0.05,maxpsffraction=0.8,interactive=False,nmajor=-1,fullsummary=False,usemask='user',mask='final_mask_I_sorted'+suffix+'.mask',restart=True,savemodel='none',calcres=False,calcpsf=False,psfcutoff=0.35,parallel=True )

make_clean_mask(img_name,5.0,'final_mask_I_sorted'+suffix+'_2.mask','false')


tclean(vis=img_list,selectdata=True,field='',spw='1',timerange='',uvrange='',antenna='',scan='',observation='',intent='',datacolumn='corrected',imagename=img_name,imsize=imsz,cell='0.15arcsec',start=min_chan,nchan=int(max_chan-min_chan),phasecenter=phase_c,stokes='I',projection='SIN',startmodel='',specmode='cube',reffreq='',outframe='',veltype='radio',restfreq='230.538GHz',interpolation='linear',perchanweightdensity=True,gridder='mosaic',facets=1,psfphasecenter='',wprojplanes=1,vptable='',mosweight=True,aterm=True,psterm=False,wbawp=True,conjbeams=False,cfcache='',usepointing=False,computepastep=360.0,rotatepastep=360.0,pointingoffsetsigdev=[],pblimit=0.2,normtype='flatnoise',deconvolver='multiscale',scales=[0, 6, 12],nterms=2,smallscalebias=0.6,fusedthreshold=0.0,largestscale=-1,restoration=True,restoringbeam='common',pbcor=False,outlierfile='',weighting='briggs',robust=2.0,npixels=0,uvtaper=[],niter=500000,gain=0.2,threshold=str(round(1.5*rms,4))+'mJy/beam',nsigma=0.0,cycleniter=50,cyclefactor=2.0,minpsffraction=0.05,maxpsffraction=0.8,interactive=False,nmajor=-1,fullsummary=False,usemask='user',mask='final_mask_I_sorted'+suffix+'_2.mask',restart=True,savemodel='none',calcres=False,calcpsf=False,psfcutoff=0.35,parallel=True )



#### Primary beam correction
impbcor(imagename=img_name+'.image/',pbimage=img_name+'.pb/',outfile=img_name+'_pbcorr.image/',box='',region='',chans='',stokes='I',mask='',mode='divide',cutoff=0.5,stretch=False,overwrite=True)

### Export
exportfits(imagename=img_name+'_pbcorr.image/',fitsimage=img_name+'_pbcorr.fits',velocity=True,overwrite=True)
exportfits(imagename=img_name+'.pb/',fitsimage=img_name+'_pb.fits',velocity=True,overwrite=True)

#### SOFIA moments
os.system('sofia '+path_sofia_templates+'moments_par_template.par input.data='+img_name+'_pbcorr.fits'+' input.primaryBeam='+img_name+'_pb.fits output.filename='+img_name+'_pbcorr')
