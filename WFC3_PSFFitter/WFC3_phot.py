#! /usr/bin/python

import pyfits
import numpy as num
import random 
#import sys
from scipy import interpolate as interp
#import libsorter
import os
import string

#better version of centroid comes from sourceFind
from sourceFind import centroid



def grid(fn,Filter,X,Y,m,z,detector,despace,Zpars,wx,wy,rot,r,R,star,
         xMin,xMax,numx,yMin,yMax,numy,psfSubFactor):
    """
    Evaluate the chis values on a grid that is centered on the X/Y.
    It evalues between X-xmin to X+xmax,
    and Y-ymin to Y+ymax using numx and numy evaluation points 
    in the x and y axes respectively. It returns the xvalues
    and yvalues, and the chis evaluated at each x and y value.
    """

    import WFC3_PSFFits as P

    xvals=num.linspace(X-xMin,X+xMax,numx)
    yvals=num.linspace(Y-yMin,Y+yMax,numy)
    chis=[]
    for i in range(len(xvals)):
        chis.append([])
        for j in range(len(yvals)):
            minuit=P.Min4PSFFit(fn,1,[[0,xvals[i],yvals[j],m,3.5]],
                                Filter,detector=detector,useMemory='n',
                                star=star,Zpars=Zpars,
                                psfSubSampFactor=psfSubFactor)
            minuit.initMin(minx=X-r,maxx=X+R,miny=Y-r,maxy=Y+R,Z=z,
                           wx=wy,wy=wy,rot=rot,fitMask=[],despace=despace,fixFocus=True)
            res=minuit.getMin(getWithoutMin=True)
            minuit.gMin.mnrset(0)
            chis[len(chis)-1].append(float(res[2]))
            
            print
            print res[2],xvals[i],yvals[j],'***'
            print
                  

    return (xvals,yvals,num.array(chis))

def getBestPoint(x,y,chi):
    m=num.min(chi)
    w=num.where(chi==m)
    return (x[w[0][0]],y[w[1][0]],m)


def getPSFCentroid(fni,X,Y,M,z=0.0,despace=0.0,Zpars=[0.0,0.0,0.0,0.0,0.0,0.0,0.0],wx=0.44,wy=0.44,rot=0.0,
                   r=5,R=6,star=[],
                   xMin=3.,xMax=3.,numx=21,yMin=3.,yMax=3.,numy=21,bestpsfSubFactor=55):
    """
    get a psf based centroid using a number of grid evaluations.

    If you think you are more than +-3 pixels from the true
    centroid, then use something else like a photocentroid
    first!

    First do a 21x21 grid using 0.06666 pixel subsampling.
    
    """

    if '/' in fni:
        (chdir,fn)=fni.split('/')
    else:
        fn=fni
        chdir='.'

    if '.fits' not in fn: fn+='.fits'

    os.chdir(chdir)

    #get the image data
    ha=pyfits.open(fn)
    header0=ha[0].header
    header1=ha[1].header
    ha.close()

    DETECTOR=header0['DETECTOR']
    FILTER=string.strip(header0['FILTER'])
    XCEN=int(float(header1['CENTERA1']))
    YCEN=int(float(header1['CENTERA2']))
    XSIZE=int(float(header1['SIZAXIS1']))
    YSIZE=int(float(header1['SIZAXIS2']))

    #get the pixel area map
    if DETECTOR=='UVIS':
        if YCEN<2051:
            #chip 1
            detector='uv1'
        else:
            #chip 2
            detector='uv2'
    else:
        detector='ir'

    

    #the simple way of refined gridding.
    (x,y,chi)=grid(fn,FILTER,X,Y,M,z,detector,despace,Zpars,wx,wy,rot,r,R,star,
                   xMin,xMax,numx,yMin,yMax,numy,15)
    (newx1,newy1,guessChi)=getBestPoint(x,y,chi)
    (x,y,chi)=grid(fn,FILTER,newx1,newy1,M,z,detector,despace,Zpars,wx,wy,rot,r,R,star,
                   0.34,0.34,18,0.34,0.34,18,55)
    (newx2,newy2,guessChi)=getBestPoint(x,y,chi)
    (x,y,chi)=grid(fn,FILTER,newx2,newy2,M,z,detector,despace,Zpars,wx,wy,rot,r,R,star,
                   0.07,0.07,15,0.07,0.07,15,105)
    (newx3,newy3,guessChi)=getBestPoint(x,y,chi)
    cent_x=newx3
    cent_y=newy3



    #use below for the variable centroiding precision.
    """


    while bestpsfSubFactor>currentSubSampFactor:
        oldStepSizeX=(xMax+xMin)/float(numx)
        oldStepSizeY=(yMax+yMin)/float(numy)

        currentSubSampFactor*=4
        if currentSubSampFactor>105: currentSubSampFactor=105


        newStepX=num.max(oldStepSizeX/5.,2./currentSubSampFactor) #ensure the step is the inverse of the subsample factor
        newStepY=num.max(oldStepSizeY/5.,2./currentSubSampFactor)

        (x,y,chi)=grid(fn,FILTER,newx1,newy1,M,z,detector,despace,Zpars,wx,wy,rot,r,R,star,
                       newStepX*3.,newStepX*3.,7,newStepY*3.,newStepY*3.,7,currentSubSampFactor)
        (newx1,newy1,guessChi)=getBestPoint(x,y,chi)
    """

    if chdir<>'.':
        os.chdir('..')

    return (cent_x,cent_y)


def phot(im,coords,radii,
         useOtherData=None,
         aperCorrs=[],
         refImDir='refImages',
         cent='centroid',
         xxx='',
         linearInterp=False,
         psfInterp=False):

    """Aperture photometry of WFC3 ir/uv1,2 images.

    This code uses iraf's phot package to perform
    simple aperture photometry.

    It then returns a list of arrays, where each entry 
    is the magnitudes and errors for every coordinate pair 
    given.

    It can do three different versions. The most
    basic (and default) is photometry where
    all pixels with imquality flags>0 are 
    ignored. This is the fastest, but typically 
    produces faintward biased answers.

    Another technique is the use of linear interpolation
    using the nearest four good pixels in a cross pattern
    that bound the bad pixel. This produces less biased
    answers the the simple verion. WARNING: not very well
    tested! This is set with linearInterp=True

    The last and most reliable technique (also the slowest)
    uses tinytim psfs to perform the interpolation. All
    good pixels within the aperture are used to match the
    psf flux to that observed. The bad pixels are sampled
    off the psf and scaled to the correct flux to produce
    the corrected image from which the photometry is 
    acquired. This is set with psfInterp=True

    Inputs:
    im - the image filename without the '.fits'.
    coords - an array of coordinates in form [x,y]. These
             are in iraf coordinates ie. (+1,+1) from pyfits
             coordinates!
    radii - an array of radii [r1,r2,r3,..] float acceptable.
    aperCorrs - an array of aperture corrections to apply.
    refImDir='refImages' - the directory with reference images.
                           This will actually be ../refImages
    cent='centroid' - perform an iraf centroid as default.
                      set to 'none' if no centroid desired.
                      
                      if centroiding, returns IRAF type 
                      coordinates.

    xxx='' - a string to add to the corrected image filename.
    """

    if psfInterp:
        import WFC3_PSFFits as P
        import sys

    from pyraf import iraf

    ###for pyraf
    yes=1
    no=0

    #import the necessary packages
    iraf.digiphot()
    iraf.daophot()
    
    #set local parameters
    iraf.set(uparm="./")
    
    iraf.daophot.verbose=no
    iraf.daophot.verify=no
    iraf.daophot.update=no
    
    iraf.centerpars.calgori=cent
    iraf.centerpars.maxshift=2.
    iraf.centerpars.minsnra=1.0
    #iraf.daopars.fitrad=1.0      #enough to leave set
    
    iraf.photpars.apertur=2
    iraf.photpars.zmag=0.0
    
    iraf.datapars.fwhmpsf=2.
    iraf.datapars.sigma=0.0
    iraf.datapars.datamin=-1000
    iraf.datapars.datamax=45000
    iraf.datapars.ccdread=""
    iraf.datapars.gain=""
    iraf.datapars.readnoi=0.0
    iraf.datapars.epadu=1.0
    iraf.datapars.exposure=""
    iraf.datapars.itime=1.0
    iraf.datapars.airmass=""
    iraf.datapars.filter=""
    iraf.datapars.obstime=""
    
    iraf.fitskypars.salgori='mode'
    iraf.fitskypars.annulus=max(radii)*1.2
    iraf.fitskypars.dannulu=8
    
    #save the new parameters
    iraf.datapars.saveParList()
    iraf.photpars.saveParList()
    #iraf.daopars.saveParList()
    iraf.centerpars.saveParList()
    iraf.daophot.saveParList()
    iraf.fitskypars.saveParList()



    if '.fits' in im: fn=im
    else: fn=im+'.fits'
    
    #get the image data
    ha=pyfits.open(fn)
    header0=ha[0].header
    header1=ha[1].header
    data=ha[1].data.astype('f')
    error=ha[2].data
    qual=ha[3].data

    if useOtherData<>None: data=useOtherData
    (a,b)=data.shape
        
    DETECTOR=header0['DETECTOR']
    FILTER=string.strip(header0['FILTER'])
    GAIN=float(header0['CCDGAIN'])
    EXPTIME=float(header0['EXPTIME'])
    if DETECTOR=="UVIS":
        PHOTZPT=float(header1['PHOTZPT'])
        PHOTFLAM=float(header1['PHOTFLAM'])
        data/=EXPTIME
    else:
        PHOTZPT=float(header0['PHOTZPT'])
        PHOTFLAM=float(header0['PHOTFLAM'])
    XCEN=int(float(header1['CENTERA1']))
    YCEN=int(float(header1['CENTERA2']))
    XSIZE=int(float(header1['SIZAXIS1']))
    YSIZE=int(float(header1['SIZAXIS2']))

    #headers are fucked in some of the ir images. This is a 
    #cluge to fix that
    if XSIZE==522: XSIZE=512
    if YSIZE==522: YSIZE=512


    ###get the flat
    #flathan=pyfits.open('flat_'+FILTER+'.fits')
    #flatdat=flathan[1].data
    #flathan.close()

    #get the pixel area map
    if DETECTOR=='UVIS':
        if YCEN<2051:
            #chip 1
            name='UVIS1wfc3_map.fits'
            detector='uv1'
        else:
            #chip 2
            name='UVIS1wfc3_map.fits'
            detector='uv2'
    else:
        name='ir_wfc3_map.fits'
        detector='ir'
    PAMhan=pyfits.open(refImDir+'/'+name)
    (A,B)=PAMhan[1].data.shape
    if (a==YSIZE or a==A):
        my=0
        My=a
    else:
        if detector=='uv1' or detector=='uv2':
            my=YCEN-YSIZE/2-2051-1
            My=my+YSIZE+1
        else:
            my=YCEN-YSIZE/2-1
            My=my+YSIZE+1

    if (b==XSIZE or b==B):
        mx=0
        Mx=b
    else:
        mx=XCEN-XSIZE/2-1
        Mx=mx+XSIZE+1
    PAMdata=PAMhan[1].data.astype('f')[my:My,mx:Mx]
    data*=PAMdata
    #data/=flatdat

    #get the hand-selected bad-pixel map
    badHan=open(im+'.fits.badPix')
    dunk=badHan.readlines()
    badHan.close()
    badPix=[]
    for ii in range(len(dunk)):
        s=dunk[ii].split()
        xbp=int(float(s[0]))
        ybp=int(float(s[1]))
        for yy in range(-1,2):
            for xx in range(-1,2):
                badPix.append([xbp-1+yy,ybp-1+xx])

    #fill in all the bad pixels using the median bg value around the images 
    #edges (4 pixels wide) and with a linear interpolation within the 
    #central area near each individual input coordinate
    #ldata=data.reshape(a*b)*1.0
    #ldata.sort()
    #median=ldata[a*b/2]
    
    median=num.median(data)
    noise=num.std(data)
    
    """
    for i in range(a):
        for j in range(4):
            if qual[i][j]>0:
                data[i][j]=median
            if qual[i][b-j-1]>0:
                data[i][b-j-1]=median
    for i in range(b):
        for j in range(4):
            if qual[j][i]>0:
                data[j][i]=median
            if qual[a-j-1][i]>0:
                data[a-j-1][i]=median
    """
                
    

    if psfInterp:
        for kk in range(len(coords)):
            xp=int(coords[kk][0]-1+0.5)
            yp=int(coords[kk][1]-1+0.5)
            R=int(max(radii))+2

            x=[]
            y=[]
            gx=[]
            gy=[]
            for i in range(yp-R,yp+R):
                for j in range(xp-R,xp+R):
                    if (j<0 or j>b or i<0 or i>a): continue
                    
                    if qual[i][j]>0 or ([j,i] in badPix):

                        x.append(j)
                        y.append(i)
                    else:
                        gx.append(j)
                        gy.append(i)
            
            bad_x=num.array(x)
            bad_y=num.array(y)
            good_x=num.array(gx)
            good_y=num.array(gy)
            
            if len(x)==0: continue

            [d,f]=im.split('/')
            os.chdir(d)
            
            FILTER=FILTER.replace('W','w').replace('M','m')
            obj=P.PSFFit(f+'.fits',1,useMemory='n')
            obj.initObject(0,coords[kk][0]-1,coords[kk][1]-1,1,3.5)
            try: os.remove('genImage.fits')
            except: pass
            try: os.remove('genImagePSF.fits')
            except: pass
            obj.produceImage(write='y')
            
            psfHan=pyfits.open('genImage.fits')
            psfData=(psfHan[0].data/EXPTIME)*PAMdata
            psfHan.close()
            
            
            try: os.remove('genImage.fits')
            except: pass
            try: os.remove('genImagePSF.fits')
            except: pass
            
            os.chdir('../')
            
            
            xmin=max(0,xp-R)
            xmax=min(b,xp+R)
            ymin=max(0,yp-R)
            ymax=min(a,yp+R)


            dcut=data[good_y,good_x]
            psfCut=psfData[good_y,good_x]
            w=num.where(dcut>median+2.*median**0.5)[0]
            multi=num.sum(dcut[w]/psfCut[w])/len(w)
            
            data[bad_y,bad_x]=psfData[bad_y,bad_x]*multi+median

    if linearInterp: #probably will fail if the peak of the psf is on a bad pixel!
        for kk in range(len(coords)):                    
            xp=int(coords[kk][0]-1+0.5)
            yp=int(coords[kk][1]-1+0.5)
            R=int(max(radii))+2
            for i in range(yp-R,yp+R):
                for j in range(xp-R,xp+R):
                    if (j<0 or j>b or i<0 or i>a): continue
                    
                    if qual[i][j]>0 or ([j,i] in badPix):
                        x=[]
                        y=[]
                        z=[]
                        val=-1
                        c=-1
                        while val==-1:
                            if qual[i+c][j]==0:
                                val=data[i+c][j]
                            else:
                                c-=1
                        x.append(j)
                        y.append(i+c)
                        z.append(val)
                        
                        val=-1
                        c=1
                        while val==-1:
                            if qual[i+c][j]==0:
                                val=data[i+c][j]
                            else:
                                c+=1
                        x.append(j)
                        y.append(i+c)
                        z.append(val)
                        
                        val=-1
                        c=-1
                        while val==-1:
                            if qual[i][j+c]==0:
                                val=data[i][j+c]
                            else:
                                c-=1
                        x.append(j+c)
                        y.append(i)
                        z.append(val)
                        
                        val=-1
                        c=1
                        while val==-1:
                            if qual[i][j+c]==0:
                                val=data[i][j+c]
                            else:
                                c+=1
                                
                        x.append(j+c)
                        y.append(i)
                        z.append(val)
                        
                        inte=interp.bisplrep(x,y,z,kx=1,ky=1)
                        value=interp.bisplev([i],[j],inte)
                        data[i][j]=value
                        #print "BORKED! Something about linear interpolation produces negative pixel values. This has not been checked to see if it's OK. It could just be noise. Though it might not be!!!"


    #write the coordinates file
    han=open('spot','w+')
    for ii in range(len(coords)):
        print >>han,coords[ii][0],coords[ii][1]
    han.close()

    #setup the necessary keywords
    iraf.datapars.gain='CCDGAIN'
    iraf.datapars.itime=1.0
    iraf.datapars.saveParList()
    iraf.photpars.zmag=PHOTZPT

    
    #write out the interpolated image for photometry

    #get centroid
    #only works on the non-phot-flam-ed image
    if cent=='centroid':

                            
        try:
            os.remove('junk'+xxx+'_cent.fits')
            os.remove('junk'+xxx+'_cent.fits1.ctr.1')
        except:
            pass
        for ii in range(a):
            ha[1].data[ii,:]=data[ii,:]
        ha.writeto('junk'+xxx+'_cent.fits')


        iraf.photpars.apertur=radii[0]
        iraf.photpars.saveParList()

        iraf.center(image='junk'+xxx+'_cent.fits[1]',coords='spot',verify=no,interactive=no)
        Zing= iraf.txdump(textfiles = 'junk'+xxx+'_cent.fits1.ctr.1', fields = 'XCENTER,YCENTER', expr=iraf.yes, Stdout=1)
        try:
            os.remove('junk'+xxx+'_cent.fits1.ctr.1')
        except:
            pass

        XC=[]
        YC=[]
        ZZ=Zing[0].split()
        for ii in range(len(ZZ)):
            if 'INDEF' not in ZZ[0]:
                XC.append(float(ZZ[0]))
                YC.append(float(ZZ[1]))
            else:
                XC.append(-32768.)
                YC.append(-32768.)

        #write the new centroided coordinates file
        han=open('spot','w+')
        for ii in range(len(coords)):
            print >>han,XC[ii],YC[ii]
        han.close()

        iraf.centerpars.calgori=cent
        iraf.centerpars.saveParList()



    data*=PHOTFLAM
    for ii in range(a):
        ha[1].data[ii,:]=data[ii,:]

    try:
        os.remove('junk'+xxx+'.fits')
    except:
        pass

    ha.writeto('junk'+xxx+'.fits')
    ha.close()



    MS=[]
    MES=[]
    for i in range(len(coords)):
        MS.append([])
        MES.append([])
        
    for kk in range(len(radii)):
        iraf.photpars.apertur=radii[kk]
        iraf.photpars.saveParList()


        try:
            os.remove('spottedm')
        except:
            pass
        iraf.daophot.phot(image='junk'+xxx+'.fits[1]', coords='spot', output= 'spottedm', interactive=no, verbose=no)


        ####magnitude equation is
        ####  m = Z-2.5*LOG(f_DN*PHOTFLAM)
        ####magnitude error is then
        ####  dm=1.09/sqrt(f_DN*EXPTIME*GAIN)
        #### Remember that to get f back in ADU we need to divide by PHOTFLAM
        #### 0.95 included to consider approximate qe effect
        Zing= iraf.txdump(textfiles = 'spottedm', fields = 'MAG,FLUX,SUM', expr=iraf.yes, Stdout=1)

        for i in range(len(MS)):
            ZZ=Zing[i].split()
            #print ZZ
            if 'INDEF' not in ZZ[0]:
                magVal=float(ZZ[0])
                if aperCorrs<>[]: magVal-=aperCorrs[i]
                MS[i].append(magVal)
                
                flux=float(ZZ[1])/PHOTFLAM
                dflux=num.sqrt(float(ZZ[2])/PHOTFLAM)

                #print flux,dflux
                MES[i].append(1.09*dflux/(flux*num.sqrt(GAIN*EXPTIME/0.95)))

            else:
                MS[i].append(-32768)
                MES[i].append(-32768)

    try:
        os.remove('junk'+xxx+'_cent.fits')
    except:
        pass
    try:
        os.remove('junk'+xxx+'.fits')
    except:
        pass
                    

    if cent<>'none':
        return (MS,MES,median,XC,YC)
    else:
        return (MS,MES,median)




def getAperCorr(im,coords,radii,
                despace=0.0,
                xxx='',
                refImDir='refImages'):

    import WFC3_PSFFits as P
    #import sys


    #get the image data
    ha=pyfits.open(im+'.fits')
    header0=ha[0].header
    header1=ha[1].header
    ha.close()

    DETECTOR=header0['DETECTOR']
    FILTER=string.strip(header0['FILTER']).replace('W','w').replace('M','m')
    XCEN=int(float(header1['CENTERA1']))
    YCEN=int(float(header1['CENTERA2']))
    XSIZE=int(float(header1['SIZAXIS1']))
    YSIZE=int(float(header1['SIZAXIS2']))



    if DETECTOR=='UVIS':
        if YCEN<2051:
            #chip 1
            detector='uv1'
        else:
            detector='uv2'
    else:
        detector='ir'

    [d,f]=im.split('/')
    os.chdir(d)
    psfDir='PSFS_'+FILTER+'+'
    
    #setup the object
    obj=P.PSFFit(f+'.fits',1,FILTER,detector=detector,useMemory='n')

    aperCorrs=[]
    for kk in range(len(coords)):
        aperCorrs.append([])

        obj.initObject(0,coords[kk][0]-1,coords[kk][1]-1,1,3.5)
        obj.produceImagePSF(despace=despace)
        psfName=obj.objects[0]['psfName']
        

        han=pyfits.open('../'+psfDir+'/'+psfName)
        data=han[0].data
        (a,b)=data.shape
        han.close()
        
        for rr in range(len(radii)):
            s=0.0
            for i in range(a):
                for j in range(b):
                    delta=((i-a/2)**2+(j-b/2)**2)**0.5
                    if delta<=radii[rr]*5:
                        s+=data[i,j]
            s/=num.sum(data)

            aperCorrs[kk].append(-2.5*num.log(s)/num.log(10))

        os.chdir('../')


    return aperCorrs


