#! /usr/bin/python


######
#This bit includes the code necessary to run an PSFfit on WFC3 data.
#When initiallizing the code, coordinates should be in numpy values. That is,
#the first pixel has (0,0) rather than (1,1) as in the fits format. This is 
#the case for ALL coordinates input by hand. That is, functions initObject,
#and initMin all take as input NUMPY coordinates (not image coordinates).
#getMin also returns in NUMPY coordinates not image coordinates.
#
#as an example, if one creates a fake star with coordinates 200,300 
#using initObject, produceImage, and genModelImage, then the fake 
#image will have star centered at 201,301.
#
#To remove that star using initObject,produceImage, and imageDiff,
#one would provice to initObject 200,300, the numpy coordinates.
#
#Then to fit that star using initMin, runMin, getMin,
#one would provide as input the NUMPY coordinates 200,300, and
#getMin would return NUMPY coordinates.
#
###############################################################
#Final note about produceImage and genModelImage, YOU MUST SET#
#depace and dE explicitly here if you don't want it to default# 
#to despace=0.0 and dE to the header kernel!!!                #
###############################################################
#
#####

import numpy as num
#from numpy.fft import fft2, ifft2
import commands
import os
import glob

import ROOT   
from array import array  #needed by root for the input arglists

import pyfits

from scipy.optimize import optimize
from scipy import special as spec

import random

import sys

import guass as GU

from scipy import interpolate as interp,signal

import distutils.spawn

#Kernel=num.array([[0.0125,0.05,0.0125],[0.05,0.75,0.05],[0.0125,0.05,0.0125]])
psfsInMemory={}
psfHeadersInMemory={}

def initMemory(psfDir):
    print "Loading PSFS into memory...",
    psfFiles=glob.glob(psfDir+'/*fits')

    imDat=[]
    for ii in psfFiles:
        if ii not in psfsInMemory:
            imHan=pyfits.open(ii)
            imDat=imHan[0].data
            headCards=imHan[0].header.ascardlist()
            imHan.close()
            
            psfsInMemory[ii]=imDat*1.0
            psfHeadersInMemory[ii]=headCards[:]
    print ".done."
    print
    return

def rint(x):
    return int(x+0.5)






def getBG(image,x,y,delta,close=12):
    if type(image)==type('1'):
        han=pyfits.open(image)
        image=han[1].data
        han.close()

    
    a=image[y-delta:y+delta,x-delta:x+delta]
    b=[]
    for qu in range(a.shape[0]):
        for quu in range(a.shape[1]):
            if abs(delta-qu)>close and abs(delta-quu)>close:
                b.append(a[qu][quu])
    
    av=num.sum(b)/(len(b))
    std=num.sqrt(num.sum((b-av)**2)/len(b))

    for qu in range(5):

        absB=num.abs(b-av)
        
        bunk=num.less(absB,2.0*std)
        
        newAv=num.sum(bunk*b)/num.sum(bunk)
        newStd=num.sqrt(num.sum((bunk*(b-av))**2)/num.sum(bunk))

        av=newAv
        std=newStd

    
    
    vals=[]
    for qu in range(4*delta):
        rx=int(random.random()*(len(b)-1))

        vals.append(b[rx])
    vals.sort()

    return (av,vals[delta])


def convolve2d(image1, kernel):
    (y1,x1)=image1.shape
    (y2,x2)=kernel.shape

 
    padImage=num.zeros((y1+2*y2,x1+2*x2),dtype=image1.dtype)
    padImage[y2:y1+y2,x2:x1+x2]=image1[:y1,:x1]
    #for ci in range(y1):
    #    for cj in range(x1):
    #        padImage[ci+y2,cj+x2]=image1[ci,cj]


    convImage=padImage*0.0

    #rev=kernel*0.0
    for ci in range(y2,y1+y2):
        for cj in range(x2,x1+x2):
            
            #print ci,cj,padImage[ci:ci-y2:-1,cj:cj-x2:-1]
            convImage[ci,cj]=num.sum(kernel*padImage[ci+1:ci-y2+1:-1,cj+1:cj-x2+1:-1])
    return convImage[y2:y1+y2,x2:x2+x1]



#this is the tinytim execution command
# It defaults to a blackbody with Solar Temperature (5770 K) and an non-over
# sampled psf. It also assumes it is generating for observations after the 
# cooler was installed. Output is nameX_Y.fits
#
#it takes as input:
#    filt  - strig containing the filter name eg. 'F110W'
#    x     - float x pixel location
#    y     - float y pixel location
#    name  - output psf name without the .fits ending
#            This will also generate a 'name'.par file for tiny2
#            and an input file 'name'.inp for tiny1.
#    detector - takes uv1, uv2, or ir according to which detector is considered
#    Optional:
#    psfRad - float radius of the generated psf. Defaults to 3.0 Arcseconds
#
# Returns the name of the file just writen.
#NOTE: Pixel locations for this command are in fits array form, not the numpy array form. 
def TinyTim(filt,x,y,name,workDir,detector,psfRad=3.0,despace=0.0,Zpars=False):

    """
    this is the tinytim execution command
    It defaults to a blackbody with Solar Temperature (5770 K) and an non-over
    sampled psf. It also assumes it is generating for observations after the 
    cooler was installed. Output is nameX_Y.fits
    
    it takes as input:
        filt  - strig containing the filter name eg. 'F110W'
        x     - float x pixel location
        y     - float y pixel location
        name  - output psf name without the .fits ending
        This will also generate a 'name'.par file for tiny2
        and an input file 'name'.inp for tiny1.
        detector - takes uv1, uv2, or ir according to which detector is considered
        Optional:
        psfRad - float radius of the generated psf. Defaults to 3.0 Arcseconds

    Returns the name of the file just writen.
    NOTE: Pixel locations for this command are in fits array form, not the numpy array form. 
    """

    s=distutils.spawn.find_executable('tiny1')
    if s==None:
        print "Cannot find the Tinytim executables directory!"
        sys.exit()
    else:
        tinyPath=s.split('tiny1')[0]

    currentDir=os.getcwd()
    os.chdir(workDir)

    rdespace=round(despace,2)

    imName='psf_'+str(rint(x))+'_'+str(rint(y))+'_'+str(rdespace)+'_'+filt+'_'
    imName+=detector
    imName+='.fits'
    if Zpars:
        imName='d'+imName

    
    #print '\n\n\n                                ',rdespace

    files=glob.glob('psf*fits')
    if imName in files: #file already exists
        if not Zpars:
            os.chdir(currentDir)
            return imName
        else:#delete the old distortion image to ensure we don't open up the wrong distortion file!
            os.remove(imName)

##########
    #To be sure that despace = 0.0 doesn't occur
    #this is a hack for an unexplained behaviour in the way Tim3 interprets despace values.
    if rdespace == 0.0:
        rdespace+=0.01
##########

    if detector=='ir':
        han=open(name+'.inp','w+')
        print >>han,23
        print >>han,str(rint(x))+' '+str(rint(y))
        print >>han,filt
        print >>han,2
        print >>han,'5770'
        print >>han,str(psfRad)
        print >>han,str(rdespace)
        print >>han,name
        han.close()
    elif detector=='uv1':
        han=open(name+'.inp','w+')
        print >>han,22
        print >>han,'1'
        print >>han,str(rint(x))+' '+str(rint(y))
        print >>han,filt
        print >>han,2
        print >>han,'5770'
        print >>han,str(psfRad)
        print >>han,str(rdespace)
        print >>han,name
        han.close()
    elif detector=='uv2':
        han=open(name+'.inp','w+')
        print >>han,22
        print >>han,'2'
        print >>han,str(rint(x))+' '+str(rint(y))
        print >>han,filt
        print >>han,2
        print >>han,'5770'
        print >>han,str(psfRad)
        print >>han,str(rdespace)
        print >>han,name
        han.close()
    command=tinyPath+'tiny1 '+name+'.par <'+name+'.inp' #Canfar
    #command='/usr/local/bin/tt/tiny1 '+name+'.par <'+name+'.inp' #Local
    c1=command
    #print command
    comStr=commands.getoutput(command)

    if Zpars:
        addDistort(name+'.par',Zpars)
        
    #commands.getoutput('rm '+name+'.inp')
    
    #subsample the psf with tiny2
    command=tinyPath+'tiny2 '+name+'.par' #Canfar
    #command='/usr/local/bin/tt/tiny2 '+name+'.par' #Local
    c2=command
    comStr2=commands.getoutput(command).split('\n')

    #distort using tiny3
    command=tinyPath+'tiny3 '+name+'.par'+' SUB=5' #Canfar
    #command='/usr/local/bin/tt/tiny3 '+name+'.par'+' SUB=5' #Local
    c3=command
    comStr3=commands.getoutput(command).split('\n')

    #delete the paramater and input files
    #try:
        #os.remove(name+'.par')
        #os.remove(name+'.inp')
        #os.remove(name+'_psf.fits') ##don't think this is needed
    #except: pass
    
    
    for i in range(len(comStr3)):
        if 'Writing' in comStr3[i]:
            s=comStr3[i].split()
            fileName=s[4]
        #else:
        #    print comStr[i]
    try:        
        mvCommand='mv '+fileName+' '+imName
        commands.getoutput(mvCommand)
        #print 'Wrote '+name+'psf_'+str(x)+'_'+str(y)+'.fits'
        
        os.chdir(currentDir)
        return imName
    except:
        print c1,'&'
        print comStr
        print c2,'&'
        print comStr2
        print c3,'&'
        print comStr3
        sys.exit()

def addDistort(parFile,Zvals):
    print Zvals
    fileHan=open(parFile)
    fileData=fileHan.readlines()
    fileHan.close()

    
    for ii in range(len(fileData)):
        #if 'Z4-Z8' in fileData[ii]:continue#not sure why this line is here
        if 'Z5' in fileData[ii] and Zvals[0]<>0.:
            fileData[ii]='  %2.5f%s'%(Zvals[0],fileData[ii][9:])
            #fileData[ii]='  '+format(Zvals[0],"2.5f")+fileData[ii][9:]
        elif 'Z6' in fileData[ii] and Zvals[1]<>0.:
            fileData[ii]='  %2.5f%s'%(Zvals[1],fileData[ii][9:])
            #fileData[ii]='  '+format(Zvals[1],"2.5f")+fileData[ii][9:]
        elif 'Z7' in fileData[ii] and Zvals[2]<>0.: 
            fileData[ii]='  %2.5f%s'%(Zvals[2],fileData[ii][9:])
            #fileData[ii]='  '+format(Zvals[2],"2.5f")+fileData[ii][9:]
        elif ('Z8' in fileData[ii] and 'Set' not in fileData[ii]) and Zvals[3]<>0.:
            fileData[ii]='  %2.5f%s'%(Zvals[3],fileData[ii][9:])
            #fileData[ii]='  '+format(Zvals[3],"2.5f")+fileData[ii][9:]
        elif 'Z9' in fileData[ii] and Zvals[4]<>0.:
            fileData[ii]='  %2.5f%s'%(Zvals[4],fileData[ii][9:])
            #fileData[ii]='  '+format(Zvals[4],"2.5f")+fileData[ii][9:]
        elif 'Z10' in fileData[ii] and Zvals[5]<>0.:
            fileData[ii]='  %2.5f%s'%(Zvals[5],fileData[ii][9:])
            #fileData[ii]='  '+format(Zvals[5],"2.5f")+fileData[ii][9:]
        elif 'Z11' in fileData[ii] and Zvals[6]<>0.:
            fileData[ii]='  %2.5f%s'%(Zvals[6],fileData[ii][9:])
            #fileData[ii]='  '+format(Zvals[6],"2.5f")+fileData[ii][9:]
            
        
    fileHan=open(parFile,'w+')
    for ii in range(len(fileData)):
        print>>fileHan,fileData[ii],
    fileHan.close()

def superPSFResample(subPSFName,xCen,yCen,dE=[0.0,0.0,0.0],subSampFactor=105,useMemory='n',useHeaderKernel=False,star=False):
    """
    subSampFactor has to be a value such that subSampFactor/5 is an ODD number. eg 5, 15 (not 10), 25 (not 20) etc. 55.... 105... screw off.

    RECOMMENDED VALUES (5,15,55,105)


    This resamples the PSF according to the sub pixel position.
    
    100x pixel subsampling as suggested by "Tinytim modeling of the WFC3/IR images (Biretta 2012)
    
    **to get the resampling to use the default kernel provided by tinytim, either
    set useHeaderKernel=True orset dE=[0.0,0.0,0.0]**


    HACK! In the IR images, the PSF produced by TinyTim is off by 2 pixels for some reason.
    So we hack it by trimming the left most and bottom most two pixels from the PSF before
    we perform any resampling.
    DOENST APPEAR TO WORK!
    """
    trueSubFactor=subSampFactor/5 ###start with 5x supersampled PSFs!
    halfSSF=subSampFactor/2 ###will be rounded. eg. 105/2 =  52
    
    if useMemory == 'n' or  subPSFName not in psfsInMemory:
        print "Getting "+subPSFName+" file from hd"
        imhan=pyfits.open(subPSFName)
        data=imhan[0].data
        header=imhan[0].header.ascardlist()
        imhan.close()
        
        psfsInMemory[subPSFName]=data*1.0
        psfHeadersInMemory[subPSFName]=header[:]
    else:
        data=psfsInMemory[subPSFName]
        header=psfHeadersInMemory[subPSFName]
    
    
    #subsample the data
    data=num.repeat(num.repeat(data,trueSubFactor,axis=0),trueSubFactor,axis=1)
    (aa,bb)=data.shape

        
    #add multiple psfs to account for trailing.
    if star:
        starShiftsX=[]
        starShiftsY=[]
        for i in range(len(star)):
            starShiftsX.append(int(subSampFactor*star[i][0]))
            starShiftsY.append(int(subSampFactor*star[i][1]))
    else:
        starShiftsX=[0.]
        starShiftsY=[0.]
        
    
    if len(starShiftsX)>1:
        origData=data*1.0
        data*=0.0
        for i in range(len(starShiftsX)):
            if starShiftsY[i]>=0:
                ymin=starShiftsY[i]
                ymax=aa
                Ymin=0
                Ymax=aa-starShiftsY[i]
            else:
                ymin=0
                ymax=aa-abs(starShiftsY[i])
                Ymin=abs(starShiftsY[i])
                Ymax=aa
            if starShiftsX[i]>=0:
                xmin=starShiftsX[i]
                xmax=bb
                Xmin=0
                Xmax=bb-starShiftsX[i]
            else:
                xmin=0
                xmax=bb-abs(starShiftsX[i])
                Xmin=abs(starShiftsX[i])
                Xmax=bb
            
            #print ymin,ymax,xmin,xmax,Ymin,Ymax,Xmin,Xmax,starShiftsX,starShiftsY
            data[ymin:ymax,xmin:xmax]+=origData[Ymin:Ymax,Xmin:Xmax]
    data/=num.sum(data)


    #trim excess edge pixels to make sure we are commensurate 
    #with the subSampFactor
    trim_y=aa%subSampFactor
    trim_x=bb%subSampFactor

    trim_bottom=int(round(float(trim_y)/2))
    trim_top=trim_y-trim_bottom
    trim_left=int(round(float(trim_x)/2))
    trim_right=trim_x-trim_left

    data=data[trim_bottom:aa-trim_top,trim_left:bb-trim_right]
    (aa,bb)=data.shape


    ####below will add a point pixel with value 1
    ####Just for testing purposes.
    #(max_i,max_j)=num.unravel_index(data.argmax(),data.shape)
    #data*=0.0
    #data[max_i-subSampFactor/2:max_i+subSampFactor/2+1,max_j-subSampFactor/2:max_j+subSampFactor/2+1]=1./(subSampFactor*subSampFactor)
    ####
    

    #now trim the edges to offset the subsampled array
    offX=int(round(divmod(xCen,1)[1]*subSampFactor))
    offY=int(round(divmod(yCen,1)[1]*subSampFactor))


    trim_left=subSampFactor-offX
    trim_right=offX
    trim_bottom=subSampFactor-offY
    trim_top=offY

    data=data[trim_bottom:aa-trim_top,trim_left:bb-trim_right]
    (aa,bb)=data.shape


    #now do the binning to get to the desired shape
    data_view=data.reshape(aa//subSampFactor, subSampFactor , bb//subSampFactor, subSampFactor)

    trimmed_data=num.sum(data_view,axis=(3,1))
    trimmed_data/=num.sum(trimmed_data)

    #now convolve the entire thing with the global Kernel
    if useHeaderKernel or dE==[0.0,0.0,0.0]:
        kernel=[]
        for iiiii in range(len(header)-3,len(header)):
             #print header[iiiii]
            s=str(header[iiiii]).split()
            kernel.append([float(s[1]),float(s[2]),float(s[3])])
        kernel=num.array(kernel)
        print "Using header kernel for CTE."
    else:
        kernel=GU.getKernel(dE[0],dE[1],dE[2])
    convData=convolve2d(trimmed_data,kernel)
    return convData



def superPSFResamplewInterp(subPSFName,xCen,yCen,dE=[0.0,0.0,0.0],subSampFactor=105,useMemory='n',useHeaderKernel=False,star=False):
     """
     HASN'T BEEN CORRECTED TO ACCOUNT FOR CORRECT BINNING AS w/out interp HAS!


     This resamples the PSF according to the sub pixel position.

     100x pixel subsampling as suggested by "Tinytim modeling of the WFC3/IR images (Biretta 2012)

     **to get the resampling to use the default kernel provided by tinytim, either
     set useHeaderKernel=True orset dE=[0.0,0.0,0.0]**
     """

     print xCen,yCen

     trueSubFactor=subSampFactor/5 ###start with 5x supersampled PSFs!
     halfSSF=subSampFactor/2 ###will be rounded. eg. 105/2 =  52

     if useMemory == 'n' or  subPSFName not in psfsInMemory:
         print "Getting "+subPSFName+" file from hd"
         imhan=pyfits.open(subPSFName)
         data=imhan[0].data
         header=imhan[0].header.ascardlist()
         imhan.close()

         psfsInMemory[subPSFName]=data*1.0
         psfHeadersInMemory[subPSFName]=header[:]
     else:
         data=psfsInMemory[subPSFName]
         header=psfHeadersInMemory[subPSFName]

     (Y,X)=data.shape
     xVals=num.linspace(0,X,X*trueSubFactor)
     yVals=num.linspace(0,Y,Y*trueSubFactor)

     resampData=data*0.0


     if star:
         starShiftsX=[]
         starShiftsY=[]
         for i in range(len(star)):
             starShiftsX.append(star[i][0])
             starShiftsY.append(star[i][1])
     else:
         starShiftsX=[0.]
         starShiftsY=[0.]


     for ss in range(len(starShiftsX)):
         sx=starShiftsX[ss]
         sy=starShiftsY[ss]

         offX=(xCen-int(xCen))*5.
         offY=(yCen-int(yCen))*5.

         print offX,offY,sx,sy

         xPoints=num.arange(len(data)) + offX + sx*5.
         yPoints=num.arange(len(data)) + offY + sy*5.

         f=interp.RectBivariateSpline(xPoints,yPoints,data)#,bbox=[-10*subSampFactor,(X+10)*subSampFactor,-10*subSampFactor,(Y+10)*subSampFactor])

         #break below into sub forloops because we don't have the memory to brute force it.
         for ii in range(trueSubFactor):
             xv=xVals[ii::trueSubFactor]
             for jj in range(trueSubFactor):
                 yv=yVals[jj::trueSubFactor]
                 a=f(xv,yv)
                 a=num.clip(a,0.0,num.max(a))
                 resampData+=a
     resampData/=num.sum(resampData)

     offX=int((xCen-int(xCen))*5)
     offY=int((yCen-int(yCen))*5)

     outData=num.zeros([Y/5,X/5],dtype=num.float32)
     for jjjjj in range(Y/5):
         for iiiii in range(X/5):
             xMin=max(0,iiiii*5-offX-2)
             xMax=min(X,iiiii*5-offX+3)
             yMin=max(0,jjjjj*5-offY-2)
             yMax=min(Y,jjjjj*5-offY+3)

             npix=(xMax-xMin)*(yMax-yMin)
             if npix<=0: continue
             outData[jjjjj][iiiii]+=num.sum(resampData[yMin:yMax,xMin:xMax])*25./npix
             #this should be commented out for speed!
             #if outData[jjjjj][iiiii]>1e14:
             #    print data[yMin:yMax,xMin:xMax]
             #    print xMin,xMax,yMin,yMax,npix,y,x,npix
             #    sys.exit()
     
     #print num.max(outData)
     #print num.sum(outData)
     #sys.exit()
     #now convolve the entire thing with the global Kernel
     if useHeaderKernel or dE==[0.0,0.0,0.0]:
         kernel=[]
         for iiiii in range(len(header)-3,len(header)):
             #print header[iiiii]
             s=str(header[iiiii]).split()
             kernel.append([float(s[1]),float(s[2]),float(s[3])])
         kernel=num.array(kernel)
         print "Using header kernel for CTE."
     else:
         kernel=GU.getKernel(dE[0],dE[1],dE[2])
     convData=convolve2d(outData,kernel)
     return convData




 #to get the resampling to use the default kernel provided by tinytim, either
 #set useHeaderKernel=True or
 #set dE=[0.0,0.0,0.0]
def resample(subPSFName,xCen,yCen,dE=[0.0,0.0,0.0],useMemory='n',useHeaderKernel=False,star=False):
    if useMemory == 'n' or  subPSFName not in psfsInMemory:
        print "Getting "+subPSFName+" file from hd"
        imhan=pyfits.open(subPSFName)
        data=imhan[0].data
        header=imhan[0].header.ascardlist()
        imhan.close()
        
        psfsInMemory[subPSFName]=data*1.0
        psfHeadersInMemory[subPSFName]=header[:]
    else:
        data=psfsInMemory[subPSFName]
        header=psfHeadersInMemory[subPSFName]

    offX=int((xCen-int(xCen))*5)
    offY=int((yCen-int(yCen))*5)
    
    #if offX>2: offX-=5
    #if offY>2: offY-=5
    (Y,X)=data.shape
    
    if star:
        starShiftsX=[]
        starShiftsY=[]
        for i in range(len(star)):
            starShiftsX.append(int(5*star[i][0]))
            starShiftsY.append(int(5*star[i][1]))
    else:
        starShiftsX=[0]
        starShiftsY=[0]
        
    if len(starShiftsX)>1:
        origData=data*1.0
        data*=0.0
        for i in range(len(starShiftsX)):
            ymin=max(0,-starShiftsY[i])
            ymax=min(Y,Y-starShiftsY[i])
            xmin=max(0,-starShiftsX[i])
            xmax=min(X,X-starShiftsX[i])
            Ymin=max(0,starShiftsY[i])
            Ymax=min(Y,Y+starShiftsY[i])
            Xmin=max(0,starShiftsX[i])
            Xmax=min(X,X+starShiftsX[i])
            data[ymin:ymax,xmin:xmax]+=origData[Ymin:Ymax,Xmin:Xmax]
    data/=num.sum(data)
    
    outData=num.zeros([Y/5,X/5],dtype=num.float32)
    for jjjjj in range(Y/5):
        for iiiii in range(X/5):
            xMin=max(0,iiiii*5-offX-2)
            xMax=min(X,iiiii*5-offX+3)
            yMin=max(0,jjjjj*5-offY-2)
            yMax=min(Y,jjjjj*5-offY+3)
            
            npix=(xMax-xMin)*(yMax-yMin)
            outData[jjjjj][iiiii]+=num.sum(data[yMin:yMax,xMin:xMax])*25./npix
            #this should be commented out for speed!
            #if outData[jjjjj][iiiii]>1e14:
            #    print data[yMin:yMax,xMin:xMax]
            #    print xMin,xMax,yMin,yMax,npix,y,x,npix
            #    sys.exit()

    #now convolve the entire thing with the global Kernel
    if useHeaderKernel or dE==[0.0,0.0,0.0]:
        kernel=[]
        for iiiii in range(len(header)-3,len(header)):
            #print header[iiiii]
            s=str(header[iiiii]).split()
            kernel.append([float(s[1]),float(s[2]),float(s[3])])
        kernel=num.array(kernel)
        print "Using header kernel for CTE."
    else:
        kernel=GU.getKernel(dE[0],dE[1],dE[2])
    convData=convolve2d(outData,kernel)
    return convData


def fcn(npar, gin, f, par, iflag):
     print "++++++++++"

     #background parameters
     Z=par[0]

     #image region to consider
     minx=int(par[1])
     maxx=int(par[2])
     miny=int(par[3])
     maxy=int(par[4])

 ##########
     #dE=[par[5],par[6],par[7]]
     dE=[par[6],par[6],0.] #Test with symetric kernel.
 ##########

     nObj=int(par[8])

     despace=par[9]
     Z5=par[10]
     Z6=par[11]
     Z7=par[12]
     Z8=par[13]
     Z9=par[14]
     Z10=par[15]
     Z11=par[16]

     #initialize the objects with the parameter guesses
     for ici in range(nObj):
         x=par[17+ici*3]
         y=par[17+ici*3+1]
         M=par[17+ici*3+2]
         target.initObject(ici,x,y,M,target.objects[ici]['radPSF'],star=target.objects[ici]['star'])
         #print x,y,M,nObj,target.objects[ici]['radPSF'],'***'

     if Z5<>0 or Z6<>0 or Z7<>0 or Z8<>0 or Z9<>0 or Z10<>0 or Z11<>0:
         target.initZpars(Z5,Z6,Z7,Z8,Z9,Z10,Z11)

     #generate the image
     target.produceImage(Z=Z,dE=dE,despace=despace)

     #open the gen'd image and the target image, and calculate the difference
     realIm=target.realImData
     tarIm=target.returnImage()
     weightIm=target.realWeightData
     killMask=target.killMask
     eKillMask=num.equal(killMask[miny:maxy,minx:maxx],0.0)

     chi=num.sum(((tarIm[miny:maxy,minx:maxx]-realIm[miny:maxy,minx:maxx])**2/(weightIm[miny:maxy,minx:maxx]**2))*eKillMask)
     npixel=num.sum(eKillMask)
     #do it by hand
     #chi=0.0
     #npixel=0.
     #for ii in range(miny,maxy):
     #    for jj in range(minx,maxx):
     #        if not killMask[ii][jj]:
     #            #print tarIm[ii][jj],realIm[ii][jj],weightIm[ii][jj]
     #            chi+=((tarIm[ii][jj]-realIm[ii][jj])**2)/((weightIm[ii][jj])**2)
     #            npixel+=1

     out=chi#/(npixel-target.freePars)

     if target.bestChi>out:
         target.bestChi=out
         outPars=[Z,minx,maxx,miny,maxy,dE[0],dE[1],dE[2],nObj,despace,Z5,Z6,Z7,Z8,Z9,Z10,Z11]
         for ici in range(nObj):
             x=par[17+ici*3]
             y=par[17+ici*3+1]
             M=par[17+ici*3+2]
             outPars.append(x)
             outPars.append(y)
             outPars.append(M)
         outPars.append('*')
         target.bestPars=outPars
     f[0]=out/(npixel-target.freePars)

     print 
     print f[0]
     print str(f[0])[:7],  str(Z)[:7], str(dE[0])[:7],str(dE[1])[:7],str(dE[2])[:7],str(despace)[:4],
     if Z5<>0 or Z6<>0 or Z7<>0 or Z8<>0 or Z9<>0 or Z10<>0 or Z11<>0:
         print
         print str(Z5)[:7],str(Z6)[:7],str(Z7)[:7],str(Z8)[:7],str(Z9)[:7],str(Z10)[:7],str(Z11)[:7],
     for boob in range(target.nObj):
         print str(par[17+boob*3])[:7],str(par[17+boob*3+1])[:7],str(par[17+boob*3+2])[:7],
     print
     #sys.exit()
     print "----------"
     print

def fcnEval(par):

     #background parameters
     Z=par[0]

     #image region to consider
     minx=int(par[1])
     maxx=int(par[2])
     miny=int(par[3])
     maxy=int(par[4])

     dE=[par[5],par[6],par[7]]

     nObj=int(par[8])

     despace=par[9]
     Z5=par[10]
     Z6=par[11]
     Z7=par[12]
     Z8=par[13]
     Z9=par[14]
     Z10=par[15]
     Z11=par[16]

     #initialize the objects with the parameter guesses
     for ai in range(nObj):
            x=par[17+ai*3]
            y=par[17+ai*3+1]
            M=par[17+ai*3+2]


            target.initObject(ai,x,y,M,target.objects[ai]['radPSF'],star=target.objects[ai]['star'])

     if Z5<>0 or Z6<>0 or Z7<>0 or Z8<>0 or Z9<>0 or Z10<>0 or Z11<>0:
         target.initZpars(Z5,Z6,Z7,Z8,Z9,Z10,Z11)

     #generate the image
     target.produceImage(Z=Z,dE=dE,despace=despace)


     #get chi squared
     chi=0.0
     realIm=target.realImData
     tarIm=target.returnImage()
     weightIm=target.realWeightData
     killMask=target.killMask
     eKillMask=num.equal(killMask[miny:maxy,minx:maxx],0.0)

     chi=num.sum(((tarIm[miny:maxy,minx:maxx]-realIm[miny:maxy,minx:maxx])**2/(weightIm[miny:maxy,minx:maxx]**2))*eKillMask)
     npixel=num.sum(eKillMask)

     return (chi/(npixel-target.freePars),npixel-target.freePars)



def fcnScipy(par,Z=0.0,minx=0,maxx=0,miny=0,maxy=0,nObj=0):
     "Not currently setup to handle aberation fitting"

     #background parameters
     #Z=args[0]
     #
     #image region to consider
     #minx=int(args[5])
     #maxx=int(args[6])
     #miny=int(args[7])
     #maxy=int(args[8])
     #nObj=int(par[9])

     dE=[par[0],par[1],par[2]]


     #initialize the objects with the parameter guesses
     for ai in range(nObj):
            x=par[3+ai*3]
            y=par[3+ai*3+1]
            M=par[3+ai*3+2]


            target.initObject(ai,x,y,M,target.objects[ai]['radPSF'])

     #generate the image
     target.produceImage(Z=Z,dE=dE)

     #get chi squared
     chi=0.0
     realIm=target.realImData
     tarIm=target.returnImage()
     weightIm=target.realWeightData
     killMask=target.killMask
     eKillMask=num.equal(killMask[miny:maxy,minx:maxx],0.0)

     chi=num.sum(((tarIm[miny:maxy,minx:maxx]-realIm[miny:maxy,minx:maxx])**2/(weightIm[miny:maxy,minx:maxx]**2))*eKillMask)
     npixel=num.sum(eKillMask)
     chi/=(npixel-target.freePars)
     print chi
     return chi

class PSFFit:

     def __init__(self,realImage,nObj=1,filter='',detector='',weightCut=10000.0,despace=0.0,refImDir='../refImages',useMemory='n',subSampFactor=55):

         #the psf subsampling factor (not the one used in TinyTim, which we leave at 5)
         self.subSampFactor=subSampFactor

         #tinytim random number
         #generate to prevent many concurrent instances from interfering
         #with one another
         self.ttRandom=int(random.random()*10000)

         if '.fits' not in realImage:
             self.realImageHan=pyfits.open(realImage+'.fits')
             self.filename=realImage+'.fits'
         else:
             self.realImageHan=pyfits.open(realImage)
             self.filename=realImage

         #for minimum keeping, because Minuit doesn't do this very well!
         self.bestChi=100000000.
         self.bestPars=[]


         #necessary header information
         self.header=self.realImageHan[0].header
         self.header1=self.realImageHan[1].header

         if detector=='':
             if 'UVIS1' in self.header['APERTURE']:
                 self.detector='uv1'
             elif 'UVIS2' in self.header['APERTURE']:
                 self.detector='uv2'
             else:
                 self.detector=self.header['DETECTOR'].lower()
         else:
             self.detector=detector

         if self.detector=='ir':
             self.photflam=float(self.header['PHOTFLAM'])
             self.zpt=self.header['PHOTZPT']
             self.photfnu=self.header['PHOTFNU']
         else:
             self.photflam=float(self.header1['PHOTFLAM'])
             self.zpt=self.header1['PHOTZPT']
             self.photfnu=self.header1['PHOTFNU']

         self.time=self.header['TIME-OBS']
         self.exptime=self.header['EXPTIME']
         self.gain=self.header['ATODGNA']
         #readnoise is ignored for the IR channel
         self.readnoise=3.1

         #pixel window center 
         self.XCEN=self.header1['CENTERA1']
         self.YCEN=self.header1['CENTERA2']
         self.XSIZE=self.header1['SIZAXIS1']
         self.YSIZE=self.header1['SIZAXIS1']
         #if uvis, make sure we are on the correct detector!
         #this doesn't actually work as the pixels per image are always less than 2050
         #if self.detector=='uv1' and self.YCEN>2050: 
         #    print "Wrong uvis detector jackass."
         #elif self.detector=='uv2' and self.YCEN<2050: 
         #    print "Wrong uvis detector jackass."

         #get the illumination pattern map
         if self.detector=='ir':
             PAMhan=pyfits.open(refImDir+'/ir_wfc3_map.fits')
             self.PAMdata=PAMhan[1].data.astype('f')
             PAMhan.close()
         elif self.detector=='uv1':
             PAMhan=pyfits.open(refImDir+'/UVIS1wfc3_map.fits')
             self.PAMdata=PAMhan[1].data.astype('f')
             PAMhan.close()
         elif self.detector=='uv2':
             PAMhan=pyfits.open(refImDir+'/UVIS2wfc3_map.fits')
             self.PAMdata=PAMhan[1].data.astype('f')
             PAMhan.close()
             #fix the YCEN coordinate
             #dont think this needs to be done
             #if self.detector=='uv1':
             #    self.YCEN-=len(self.PAMdata)


         self.realImData=self.realImageHan[1].data
         self.realWeightData=self.realImageHan[2].data

         self.maskData=self.realImageHan[3].data

         #this is the rejection threshold, above which all noisy pixels 
         #are ignored
         self.weightCut=weightCut

         self.badPix=[]

         self.killMask=self.maskData*0.0
         for i in range(len(self.maskData)):
             for j in range(len(self.maskData[i])):
                 if (self.maskData[i][j]<>0 or self.realWeightData[i][j]==0 or self.realWeightData[i][j]>=self.weightCut ) or (num.isnan(self.realWeightData[i][j])):
                     self.realWeightData[i][j]=1.0
                     self.killMask[i][j]=1.0 #one are pixels to ignore
                     self.badPix.append([j,i])

                 #check to see if any "bloom-cross" pixels are flagged
                 if self.maskData[i][j]==48:
                     killpairs=[[max(0,i-1),j],
                                [min(len(self.maskData)-1,i+1),j],
                                [i,max(0,j-1)],
                                [i,min(len(self.maskData[0])-1,j+1)]]
                     for jqrt in killpairs:
                         self.realWeightData[jqrt[0]][jqrt[1]]=1.0
                         self.killMask[jqrt[0]][jqrt[1]]=1.0 #one are pixels to ignore
                         self.badPix.append(jqrt)

                    
         if '.fits' not in realImage:
             badPixName=realImage+'.fits.badPix'
         else:
             badPixName=realImage+'.badPix'
         hang=open(badPixName)
         dang=hang.readlines()
         hang.close()
         for i in dang:
             sss=i.split()
             x=int(float(sss[0]))
             y=int(float(sss[1]))

             for qrti in range(x-1,x+2):
                 for qrtj in range(y-1,y+2):
                     self.killMask[qrtj][qrti]=1.0
                     self.realWeightData[qrtj][qrti]=1.0
                     self.badPix.append([qrtj,qrti])

         self.numBadPix=len(self.badPix)

         #the number of objects to psf fit
         self.nObj=nObj
         self.despace=despace
         if filter=='':
             self.filter=self.header['FILTER'].strip()
         else:
             self.filter=filter
         self.dE=[0.0,0.0,0.0]
         self.Zpars=False
         self.freePars=3*self.nObj

         #the object parameters
         self.objects=[]
         for i in range(self.nObj):
             source={'x':-1000, 'y':-1000, 'M':-1000., 'psfName':'', 'foc':0., 'radPSF':0., 'genPSF':0}
             self.objects.append(source)

         #initialize the image, the psf image and the background image
         self.image={}
         self.image['image']=num.zeros(self.realImData.shape,dtype=num.float32)
         self.image['psf']=num.zeros(self.realImData.shape,dtype=num.float32)
         self.image['bg']=num.zeros(self.realImData.shape,dtype=num.float32)

         #generation flags: 0 - not yet generated
         self.gen={}
         self.gen['psf']=0
         self.gen['image']=0
         self.gen['bg']=0
         self.gen['psfWrite']=0
         self.gen['bgWrite']=0
         self.gen['imageWrite']=0

         self.bgPars={}
         self.bgPars['Z']=0.0 

         self.PSFdir='../PSFS_'+self.filter+'+'
         commands.getoutput('mkdir '+self.PSFdir)

         #if useMemory has been set, load all psfs into memory to speed up the resampling
         self.usingMemory=useMemory
         #if useMemory<>'n':
         #    initMemory(self.PSFdir)

     #initialize the particular object, this does not generate the object
     # allows for modifiying the object by the same call, and the psf gen 
     # will recognize the adjustment and recreate the necessary psfs
     def initObject(self,N,xx,yy,M,radP,
                    star=False):
         x=xx
         y=yy
 #        x = round(xx*5.)/5.
 #        y = round(yy*5.)/5.
         oldSource=self.objects[N]
         if (oldSource['x']<>x) or (oldSource['y']<>y) or (oldSource['M']<>M) or (oldSource['star']<>star):
             source={'x':float(x), 'y':float(y), 'M':float(M), 'psfName':'','radPSF':radP, 'genPSF':0, 'star':star}
             self.objects[N]=source
             self.gen['psf']=0


     def initZpars(self,Z5,Z6,Z7,Z8,Z9,Z10,Z11):
         self.Zpars=[Z5,Z6,Z7,Z8,Z9,Z10,Z11]
         for ii in range(len(self.objects)):
             self.objects[ii]['genPSF']=0


     #produce the image which is PSFS+background
     # store this in self.image
     # write this to disk if requested, and will write the bg and psf images as well
     def produceImage(self,Z=0.0,dE=[0.0,0.0,0.0], despace=0.0, write='n'):
         if self.bgPars['Z']<>Z:
             self.gen['bg']=0

         if self.dE<>dE:
             for bum in range(self.nObj):
                 self.objects[bum]['genPSF']=0
             self.dE=dE

         if self.despace<>despace:
             for bum in range(self.nObj):
                 self.objects[bum]['genPSF']=0
             self.despace=despace

         self.produceImagePSF(dE=dE, despace=despace, write=write)


         #produce the image background
         if Z<>self.bgPars['Z']:
             self.image['bg']=num.zeros(self.realImData.shape,dtype=num.float32)+Z


         self.image['image']=self.image['psf']+self.image['bg']
         self.gen['image']=1

         if write<>'n':
             commands.getoutput('rm genImage.fits')
             hdu=pyfits.PrimaryHDU(self.image['image'])
             hdulist=pyfits.HDUList(hdu)
             hdulist.writeto('genImage.fits')

             self.gen['imageWrite']=1

     def genModelImage(self,imageName,Z=0.0,dE=[0.0,0.0,0.0], despace=0.0, addNoise=True,returnForXCorrelate=False):
         self.produceImage(Z=0.0,dE=dE, despace=despace, write='n')
         self.gen['image']=0
         image=self.image['image']
         noise=self.image['image']*0.0

         if addNoise:
             (Y,X)=noise.shape
             for i in range(Y):
                 for j in range(X):
                     noise[i,j]+=random.gauss(0.0,image[i,j]**0.5)
         else:
             print "     Not adding noise to the images."

         genWeightData=(self.realWeightData**2+noise**2)**0.5
         genImage=image+noise+self.realImData
         if returnForXCorrelate: return genImage

         self.realImageHan[1].data=genImage
         self.realImageHan[2].data=genWeightData

         if '.fits' in imageName:
             commands.getoutput('rm '+imageName)
             self.realImageHan.writeto(imageName)
         else:
             commands.getoutput('rm '+imageName+'.fits')
             self.realImageHan.writeto(imageName+'.fits')


     #return the generated image
     # return None if the image has not yet been created
     def returnImage(self):
         if self.gen['image']:
             return self.image['image']
         else:
             return None


     def imageDiff(self,name='',write='n'):
         if self.gen['image']:
             print
             print 'Creating the difference.'

             self.diffIm=self.image['image']*0.0
             self.diffIm+=(self.realImData-self.image['image'])*num.equal(0,self.killMask)
             if write<>'n':
                 print 'Writing the image'
                 os.system('rm '+name)
                 self.realImageHan[1].data=self.diffIm
                 if '.fits' in name:
                     self.realImageHan.writeto(name)
                 else:
                     self.realImageHan.writeto(name+'.fits')
                 #hdu=pyfits.PrimaryHDU(self.diffIm)
                 #hdulist=pyfits.HDUList(hdu)
                 #hdulist.writeto(name)
             print
         else:
             print
             print 'The image has not been generated yet!'
             print
         return


     def produceImagePSF(self,dE=[0.0,0.0,0.0],despace=0.0,write='n',templatePSF=False):        
         print 'Initializing PSF'

         self.image['psf']=num.zeros(self.realImData.shape,dtype=num.float32)        
         #generate the psfs and get their names
         for cci in range(self.nObj):
             if self.objects[cci]['genPSF']==0:
                 if not self.objects[cci]['star']:
                     #add the window center minus the image shape to get the correct psf
                     x=self.objects[cci]['x']+self.XCEN-self.XSIZE/2
                     y=self.objects[cci]['y']+self.YCEN-self.YSIZE/2
                     r=self.objects[cci]['radPSF']

                     print 'Creating psf '+str(cci)+': '+str(x)+' '+str(y)+' '+str(self.objects[cci]['M'])+' '+str(despace)+'...'

                     #plus one added to the coordinates here to get into the fits coordinates
                     self.objects[cci]['psfName']=TinyTim(self.filter,x+1,y+1,str(cci+self.ttRandom),self.PSFdir,self.detector,psfRad=r,despace=despace,Zpars=self.Zpars)

                     print '               done!'

                     self.objects[cci]['genPSF']=1
                 else:
                     psfNames=[]
                     print "Generating PSF for central position of Star"
                     l=len(self.objects[cci]['star'])/2
                     #add the window center minus the image shape to get the correct psf
                     x=self.objects[cci]['x']+self.XCEN-self.XSIZE/2+self.objects[cci]['star'][l][0]
                     y=self.objects[cci]['y']+self.YCEN-self.YSIZE/2+self.objects[cci]['star'][l][1]
                     r=self.objects[cci]['radPSF']

                     print 'Creating psf '+str(cci)+': '+str(x)+' '+str(y)+' '+str(self.objects[cci]['M'])+' '+str(despace)+'...'
                     #plus one added to the coordinates here to get into the fits coordinates
                     self.objects[cci]['psfName']=TinyTim(self.filter,x+1,y+1,str(cci+self.ttRandom),self.PSFdir,self.detector,psfRad=r,despace=despace,Zpars=self.Zpars)


                     print '               done!'

                     self.objects[cci]['genPSF']=1
             else:
                 print 'Psf '+str(cci)+': '+str(self.objects[cci]['x'])+' '+str(self.objects[cci]['y'])+' '+' '+str(self.objects[cci]['M'])+' '+str(despace)+' already created.'


         print 'Done Creating PSFs'

         print 'Combining TinyTim images.'


         print "Resampling with kernel parameters wx: "+str(dE[0])+" wy: "+str(dE[1])+" angle: "+str(dE[2])+"."
         for cci in range(self.nObj):
             imDat=superPSFResample(self.PSFdir+'/'+self.objects[cci]['psfName'],self.objects[cci]['x'],self.objects[cci]['y'],dE,useMemory=self.usingMemory,star=self.objects[cci]['star'],subSampFactor=self.subSampFactor)
             #multiply by the pixel area map
             pixMapMulti=self.PAMdata[int(self.objects[cci]['y']+self.YCEN-self.YSIZE/2),int(self.objects[cci]['x']+self.XCEN-self.XSIZE/2)]
             imDat*=pixMapMulti

             if templatePSF:
                 self.templatePSF=imDat
                 return

             (yLen,xLen)=imDat.shape


             x=int(self.objects[cci]['x'])
             y=int(self.objects[cci]['y'])

             startX=x-xLen/2
             endX=x+xLen/2
             startY=y-yLen/2
             endY=y+yLen/2

             #startX=self.objects[cci]['x']-xLen/2
             #endX=self.objects[cci]['x']+xLen/2
             #startY=self.objects[cci]['y']-yLen/2
             #endY=self.objects[cci]['y']+yLen/2


             rsX=max(0,startX)
             reX=min(len(self.image['psf'])-1,endX)
             rsY=max(0,startY)
             reY=min(len(self.image['psf'][0])-1,endY)

             """
             print '\n\n\n\n'
             print rsY,reY,rsX,reX,(yLen/2-(y-rsY)),(yLen/2+(reY-y)),(xLen/2-(x-rsX)),(xLen/2+(reX-x))
             print x,y,startX,endX,startY,endY
             print xLen,yLen
             print self.image['psf'][rsY:reY,rsX:reX].shape,(self.objects[cci]['M']*imDat[(yLen/2-(y-rsY)):(yLen/2+(reY-y)),(xLen/2-(x-rsX)):(xLen/2+(reX-x))]).shape
             print num.arange(rsY,reY),num.arange((yLen/2-(y-rsY)),(yLen/2+(reY-y)))
             print num.arange(rsX,reX),num.arange((xLen/2-(x-rsX)),(xLen/2+(reX-x)))
             print '\n\n\n\n'
             """
             self.image['psf'][rsY:reY,rsX:reX]+=self.objects[cci]['M']*imDat[(yLen/2-(y-rsY)):(yLen/2+(reY-y)),(xLen/2-(x-rsX)):(xLen/2+(reX-x))]


         if write<>'n':
             commands.getoutput('rm genImagePSF.fits')
             hdu=pyfits.PrimaryHDU(self.image['psf'])
             hdulist=pyfits.HDUList(hdu)
             hdulist.writeto('genImagePSF.fits')
             self.gen['psfWrite']=1
         else:
             self.gen['psfWrite']=0


         self.gen['psf']=1


     def getPhotometry(self):
         "not setup to handle photometry of stars"

         print 'Initializing PSF'

         self.image['psf']=num.zeros(self.realImData.shape,dtype=num.float32)        
         #generate the psfs and get their names
         for i in range(self.nObj):
             if self.objects[i]['genPSF']==0:

                 x=self.objects[i]['x']
                 y=self.objects[i]['y']
                 r=self.objects[i]['radPSF']

                 print 'Creating psf '+str(i)+': '+str(x)+' '+str(y)+' '+str(self.objects[i]['M'])+'...'

                 #plus one added to the coordinates here to get into the fits coordinates
                 self.objects[i]['psfName']=TinyTim(self.filter,x+1,y+1,str(i+self.ttRandom),self.PSFdir,self.detector,psfRad=r,despace=self.despace,Zpars=self.Zpars)

                 print '               done!'

                 self.objects[i]['genPSF']=1
             else:
                 print 'Psf '+str(i)+': '+str(self.objects[i]['x'])+' '+str(self.objects[i]['y'])+' '+str(self.objects[i]['M'])+' already created.'


         print 'Done Creating PSFs'

         print 'Getting Photometry'

         MAGS=[]
         FLUX=[]
         for i in range(self.nObj):
            imDat=superPSFResample(self.PSFdir+'/'+self.objects[i]['psfName'],self.objects[i]['x'],self.objects[i]['y'],self.dE,useMemory=self.usingMemory)
            #multiply by the pixel area map
            pixMapMulti=self.PAMdata[int(self.objects[i]['y']+self.YCEN-self.YSIZE/2),int(self.objects[i]['x']+self.XCEN-self.XSIZE/2)]
            imDat*=pixMapMulti


            (yLen,xLen)=imDat.shape
            

            startX=self.objects[i]['x']-xLen/2
            endX=self.objects[i]['x']+xLen/2
            startY=self.objects[i]['y']-yLen/2
            endY=self.objects[i]['y']+yLen/2

            x=int(self.objects[i]['x'])
            y=int(self.objects[i]['y'])

            rsX=max(0,startX)
            reX=min(len(self.image['psf']),endX)
            rsY=max(0,startY)
            reY=min(len(self.image['psf'][0]),endY)
            
            #print rsY,reY,rsX,reX,(yLen/2-(y-rsY)),(yLen/2+(reY-y)),(xLen/2-(x-rsX)),(xLen/2+(reX-x))
            flux=self.photflam*num.sum(self.objects[i]['M']*imDat[(yLen/2-(y-rsY)):(yLen/2+(reY-y)),(xLen/2-(x-rsX)):(xLen/2+(reX-x))])

            if self.detector=='uv1' or self.detector=='uv2':
                flux/=self.exptime
            FLUX.append(flux)
            mag=self.zpt-2.5*num.log(flux)/num.log(10)
            MAGS.append(mag)
         return (FLUX,MAGS)


#this is to do a grid-based binary minimization
class Binary_PSFFit:

    def __init__(self, realIm, primaryPars, filter,detector,dx,dy,useMemory='y'):
        self.targetP = PSFFit(realIm, 1, filter,despace=0.0,detector=detector,useMemory=useMemory)
        self.targetS = PSFFit(realIm, 1, filter,despace=0.0,detector=detector,useMemory=useMemory)
        bb=primaryPars[:]
        self.targetP.initObject(0,bb[1],bb[2],bb[3],bb[4])
        self.targetS.initObject(0,bb[1],bb[2],bb[3],bb[4])
        (Z,med)=getBG(self.targetP.realImData,bb[1],bb[2],30)
        self.z=Z
        (a,b)=self.targetP.realImData.shape

        self.M=bb[3]

        self.xcen=int(bb[1])
        self.ycen=int(bb[2])
        self.px=bb[1]
        self.py=bb[2]
        self.minx=min(self.xcen-dx,b)
        self.maxx=max(self.xcen+dx+1,0)
        self.miny=min(self.ycen-dy,a)
        self.maxy=max(self.ycen+dy+1,0)

    def doBinaryGrid(self,pDX,pDY,sDX,sDY,step=0.5):
        
        nxp=int(2*pDX/step+1)
        nyp=int(2*pDY/step+1)
        nxs=int(2*sDX/step+1)
        nys=int(2*sDY/step+1)

        self.primaryImages={}
        for xi in range(nxp):
            X=self.px-pDX+xi*step
            for yi in range(nyp):
                Y=self.py-pDY+yi*step

                key=str(xi)+'_'+str(yi)
                try:
                    print 'Primary ',X,Y
                    print
                    self.targetP.initObject(0,X,Y,1.0,self.targetP.objects[0]['radPSF'])
                    self.targetP.produceImage(Z=0.0,dE=[0.0,0.0,0.0])
                    self.primaryImages[key]=self.targetP.returnImage()[self.miny:self.maxy,self.minx:self.maxx]*1.0
                except: pass

        self.secondaryImages={}
        for xj in range(nxs):
            x=self.px-sDX+xj*step
            for yj in range(nys):
                y=self.py-sDY+yj*step

                key=str(xj)+'_'+str(yj)
                try:
                    print 'Secondary ',x,y
                    print
                    self.targetS.initObject(0,x,y,1.0,self.targetS.objects[0]['radPSF'])
                    self.targetS.produceImage(Z=0.0,dE=[0.0,0.0,0.0])
                    self.secondaryImages[key]=self.targetS.returnImage()[self.miny:self.maxy,self.minx:self.maxx]*1.0
                except: pass


        chiArr=[1000.]
        for xi in range(nxp):
            for yi in range(nyp):
                for xj in range(nxs):
                    for yj in range(nys):
                        for mi in range(1,41,2):
                            for Mi in range(11):
                                keyP=str(xi)+'_'+str(yi)
                                keyS=str(xj)+'_'+str(yj)
                                
                                if keyP not in self.primaryImages or keyS not in self.secondaryImages: continue
                                
                                x=self.px-sDX+xj*step
                                y=self.py-sDY+yj*step
                                X=self.px-pDX+xi*step
                                Y=self.py-pDY+yi*step
                                delta=((x-X)**2+(y-Y)**2)**0.5
                                
                                Mr=1.0
                                mr=float(mi)/50#20
                                
                                multi=1.+(Mi-5)*0.04
                                m=multi*self.M*mr/(Mr+mr)
                                M=multi*self.M*Mr/(Mr+mr)
                                
                                #m=self.M*mr/(Mr+mr)
                                #M=self.M*Mr/(Mr+mr)
 
                                tarIm=M*self.primaryImages[keyP]+m*self.secondaryImages[keyS]+self.z
                                #print num.sum(self.secondaryImages[keyS]),num.sum(self.primaryImages[keyP])
                                chi=self.fcn(tarIm)
                                if chi<=chiArr[0] and delta>1:
                                    chiArr=[chi,X,Y,M,x,y,m]
                                    print chiArr



        return chiArr
        
    def fcn(self,tarIm):
    

        #get chi squared
        chi=0.0
        realIm=self.targetP.realImData
        weightIm=self.targetP.realWeightData
        killMask=self.targetP.killMask
        eKillMask=num.equal(killMask[self.miny:self.maxy,self.minx:self.maxx],0.0)
        
        chi=num.sum(((tarIm-realIm[self.miny:self.maxy,self.minx:self.maxx])**2/(weightIm[self.miny:self.maxy,self.minx:self.maxx]**2))*eKillMask)
        npixel=num.sum(eKillMask)
        chi/=(npixel-6)#freepars is 6

        return chi



#this is to do the minimization
class Min4PSFFit:

    def __init__(self, realIm, numObj, objPars, filter='',detector='',despace=0.0,Zpars=False,useMemory='n',star=False,psfSubSampFactor=55):
        global mask
        mask=-1.0
        global target
        target =PSFFit(realIm, numObj, filter=filter,despace=despace,detector=detector,useMemory=useMemory,subSampFactor=psfSubSampFactor)

        for i in range(numObj):
            if not star:
                target.initObject(i,objPars[i][1],objPars[i][2],objPars[i][3],objPars[i][4])
            else:
                target.initObject(i,objPars[i][1],objPars[i][2],objPars[i][3],objPars[i][4],star=star)

        #flags for Minuit
        self.initMinVar=0
        self.runMinVar=0

        self.numPars=0
        self.Zpars=Zpars

        self.target=target
    """
    def crossCorrelate(self,objNum,despace,dE,subSampleFactor=15):
        genData=target.genModelImage('',Z=0.0,dE=dE, despace=despace, addNoise=False,returnForXCorrelate=True)
        imData=target.realImData*1.0

        x=target.objects[objNum]['x']
        y=target.objects[objNum]['y']
        M=target.objects[objNum]['M']
        r=target.objects[objNum]['radPSF']

        print 'Object parameters:',x,y,M,r
        xc=int((x-int(x))*subSampleFactor)+10*subSampleFactor
        yc=int((y-int(y))*subSampleFactor)+10*subSampleFactor
        c=4*subSampleFactor
        C=4*subSampleFactor+1

        genData=num.repeat(num.repeat(genData[int(y)-10:int(y)+11,int(x)-10:int(x)+11],subSampleFactor,axis=0),subSampleFactor,axis=1)[yc-c-2:yc+C+3,xc-c-2:xc+C+3]
        imData=num.repeat(num.repeat(imData[int(y)-10:int(y)+11,int(x)-10:int(x)+11],subSampleFactor,axis=0),subSampleFactor,axis=1)[yc-c:yc+C+1,xc-c:xc+C+1]

        genData/=num.sum(genData)
        imData/=num.sum(imData)
        

        (Y,X)=imData.shape
        #now do a correlation
        cross=[]
        for yi in range(0,3,1):
            cross.append([])
            for xi in range(0,3,1):
                genDataSub=genData[yi:Y+yi,xi:X+xi]
                c=num.sum(imData*num.conjugate(genDataSub))
                cross[len(cross)-1].append(c)
        cross=num.array(cross)

        if cross[1,0]<cross[1,2]: 
            print 'x too low'
        else:
            print 'x too high'

        if cross[0,1]<cross[2,0]: 
            print 'y too low'
        else:
            print 'y too high'
        print cross
        os.system('rm junk*fits')
        HDU=pyfits.PrimaryHDU(genData)
        List=pyfits.HDUList([HDU])
        HDU.writeto('junkGen.fits')
        HDU=pyfits.PrimaryHDU(imData)
        List=pyfits.HDUList([HDU])
        HDU.writeto('junk.fits')
        sys.exit()
        cross=signal.correlate(imData,genData,mode='full')
        print cross.shape
        mc=num.max(cross)
        w=num.where(cross==mc)
        print cross[58:62,58:63]
        print w
        


        return
    """


    def returnWesFit(self):
        return (target.bestChi,target.bestPars)

    #initialize the minuit object with the appropriate background and 
    #psf parameters.
    #NOTE: Must initialize the psf parameters first!!!
    def initMin(self,Z=0.0,wx=0.0,wy=0.0,rot=0.0,minx=0,maxx=256,miny=0,maxy=256,despace=0.0,steps=[],fitMask=[],fixFocus=False):

        #setup the initial step size for the minimization
        if steps==[]:
            steps=[100.,0.,0.,0.,0.,0.2,0.2,500.,0.,2.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0]  #angle step is large to get it to work!
            for i in range(target.nObj):
                steps.append(2.0)
                steps.append(2.0)
                steps.append(target.objects[i]['M']*0.15)
            self.steps=steps

        #setup the fit mask
        #values of fitmask==1 are taken as pixels to be ignored!
        if fitMask<>[]:
            print "Creating the fit mask."
            for i in range(len(fitMask)):
                for j in range(len(fitMask[i])):
                    if fitMask[i][j] or target.killMask[i][j]:
                        target.killMask[i][j]=1.0

        print 'Initializing the Minimizer'
        target.numPars=17+target.nObj*3  #5 parameters for the BG, 4 for image range, 
        #                            1 for the number of objects, and 3 for each source
        self.gMin=ROOT.TMinuit(target.numPars)

        self.gMin.SetFCN(fcn)

        self.arglist = array('d', 10*[0.])
        self.ierflg = ROOT.Long(1982)

        #self.arglist[0]=50000*25 #max calls
        #self.arglist[1]=1.0 #set error size
        #self.gMin.mnexcm("SET ERR", self.arglist, 1, self.ierflg)

        self.gMin.SetMaxIterations(50000)
        self.gMin.SetErrorDef(1.0)

        #setup a vector and args for Scipy
        self.ScipyArgs=[Z,minx,maxx,miny,maxy,target.nObj]
        self.ScipyVector=[wx,wy,rot]

        print 'Defining the initial paramters.'
        #( param num, param name, initial guess, approx error, range min, range max)
        self.gMin.DefineParameter(0, "Z", Z, steps[0], 0, 0)
        self.gMin.DefineParameter(1, "minX", minx, steps[1], 0, 0)
        self.gMin.DefineParameter(2, "maxX", maxx, steps[2], 0, 0)
        self.gMin.DefineParameter(3, "minY", miny, steps[3], 0, 0)
        self.gMin.DefineParameter(4, "maxY", maxy, steps[4], 0, 0)
        self.gMin.DefineParameter(5, "wx", wx, steps[5],0.2,1.)
        self.gMin.DefineParameter(6, "wy", wy, steps[6],0.2,1.)
        self.gMin.DefineParameter(7, "rotation", rot, steps[7], 0,90)
        self.gMin.DefineParameter(8, "nObj", target.nObj, 0, 0, 0)
        self.gMin.DefineParameter(9,"despace", despace, steps[9], -10,10)
        if not self.Zpars:
            self.gMin.DefineParameter(10,"Z5", 0.0, steps[10], -0.1,0.1)
            self.gMin.DefineParameter(11,"Z6", 0.0, steps[11], -0.1,0.1)
            self.gMin.DefineParameter(12,"Z7", 0.0, steps[12], -0.1,0.1)
            self.gMin.DefineParameter(13,"Z8", 0.0, steps[13], -0.1,0.1)
            self.gMin.DefineParameter(14,"Z9", 0.0, steps[14], -0.1,0.1)
            self.gMin.DefineParameter(15,"Z10",0.0, steps[15], -0.1,0.1)
            self.gMin.DefineParameter(16,"Z11",0.0, steps[16], -0.1,0.1)
        else:
            self.gMin.DefineParameter(10,"Z5", self.Zpars[0], steps[10], -0.1,0.1)
            self.gMin.DefineParameter(11,"Z6", self.Zpars[1], steps[11], -0.1,0.1)
            self.gMin.DefineParameter(12,"Z7", self.Zpars[2], steps[12], -0.1,0.1)
            self.gMin.DefineParameter(13,"Z8", self.Zpars[3], steps[13], -0.1,0.1)
            self.gMin.DefineParameter(14,"Z9", self.Zpars[4], steps[14], -0.1,0.1)
            self.gMin.DefineParameter(15,"Z10",self.Zpars[5], steps[15], -0.1,0.1)
            self.gMin.DefineParameter(16,"Z11",self.Zpars[6], steps[16], -0.1,0.1)

        self.gMin.FixParameter(0)
        self.gMin.FixParameter(1)
        self.gMin.FixParameter(2)
        self.gMin.FixParameter(3)
        self.gMin.FixParameter(4)
        self.gMin.FixParameter(5)
        self.gMin.FixParameter(6)
        self.gMin.FixParameter(7)
        self.gMin.FixParameter(8)
        if fixFocus:
            self.gMin.FixParameter(9)
        if not self.Zpars:
            self.gMin.FixParameter(10)
            self.gMin.FixParameter(11)
            self.gMin.FixParameter(12)
            self.gMin.FixParameter(13)
            self.gMin.FixParameter(14)
            self.gMin.FixParameter(15)
            self.gMin.FixParameter(16)
            

        self.numPars+=17

        for i in range(target.nObj):
            parName='objX'+str(i)
            xxx=target.objects[i]['x']
            yyy=target.objects[i]['y']
            MMM=target.objects[i]['M']
            self.gMin.DefineParameter(17+i*3,parName,xxx,steps[17+3*i], xxx-2, xxx+2)
            parName='objY'+str(i)
            self.gMin.DefineParameter(17+i*3+1,parName,yyy,steps[17+3*i+1],yyy-2,yyy+2)
            parName='objM'+str(i)
            self.gMin.DefineParameter(17+i*3+2,parName,MMM,steps[17+3*i+2],0,0)#0.0,max(MMM*3,20000))

            self.ScipyVector.append(xxx)
            self.ScipyVector.append(yyy)
            self.ScipyVector.append(MMM)

            self.numPars+=3

        self.initMinVar=1


    def fixPar(self,number):
        AA=self.gMin.GetNumFreePars()
        self.gMin.FixParameter(number)
        newAA=self.gMin.GetNumFreePars()
        target.freePars-=(AA-newAA)
        print "Parameter "+str(number)+"now held. "+str(newAA)+" free parameters."

    def freePar(self,number):
        if (number in [0,5,6,7,9,10,11,12,13,14,15,16]):
            AA=self.gMin.GetNumFreePars()
            self.gMin.Release(number)
            newAA=self.gMin.GetNumFreePars()
            target.freePars-=(AA-newAA)
            print "Parameter "+str(number)+"now free. "+str(newAA)+" free parameters."

    def runMin(self):
        if self.initMin:
            print "Running Minimization."
            self.runMinVar=1
            self.gMin.Migrad()#mxexcm("MIGRAD", self.arglist,2,self.ierflg)
            #self.gMin.mnsimp()
        else:
            print "Minuit not initialized."
        self.runMinVar=1

    def clearMem(self):
        print "Clearing Fit Memory"
        self.gMin.mncler()

    def runScipyMin(self):
        print optimize.fmin_powell(fcnScipy,self.ScipyVector,args=self.ScipyArgs)
        sys.exit()
        
    def runSimpMin(self):
        if self.initMin:
            print "Running Minimization."
            self.runMinVar=1
            self.gMin.mnsimp()#mxexcm("MIGRAD", self.arglist,2,self.ierflg)
            #self.gMin.mnsimp()
        else:
            print "Minuit not initialized."
        self.runMinVar=1

    def improveFit(self):
        if self.runMinVar:
            print "Trying to improve fit."
            self.runMinVar=1
            self.gMin.mnimpr()#mxexcm("MIGRAD", self.arglist,2,self.ierflg)
            #self.gMin.mnsimp()
        else:
            print "Fit not improved!"


    
        
    def CheckWesMin(self,betterPars):
        self.gMin.mnrset(0)
        self.gMin.DefineParameter(0, "Z", betterPars[0], self.steps[0], 0, 0)
        self.gMin.DefineParameter(5, "wx", betterPars[5], self.steps[5], 0.2, 0.8)
        self.gMin.DefineParameter(6, "wy", betterPars[6], self.steps[6], 0.2, 0.8)
        self.gMin.DefineParameter(7, "rotation", betterPars[7], self.steps[7], 0, 90)

        self.gMin.DefineParameter(9, "despace",betterPars[9], self.steps[9],0,0)
        self.gMin.DefineParameter(10, "despace",betterPars[10], self.steps[10],0,0)
        self.gMin.DefineParameter(11, "despace",betterPars[11], self.steps[11],0,0)
        self.gMin.DefineParameter(12, "despace",betterPars[12], self.steps[12],0,0)
        self.gMin.DefineParameter(13, "despace",betterPars[13], self.steps[13],0,0)
        self.gMin.DefineParameter(14, "despace",betterPars[14], self.steps[14],0,0)
        self.gMin.DefineParameter(15, "despace",betterPars[15], self.steps[15],0,0)
        self.gMin.DefineParameter(16, "despace",betterPars[16], self.steps[16],0,0)

        for i in range(target.nObj):
            parName='objX'+str(i)
            xxx=betterPars[17+3*i]
            self.gMin.DefineParameter(17+i*3,parName,xxx, self.steps[17+3*i],xxx-2, xxx+2)
            parName='objY'+str(i)
            yyy=betterPars[17+3*i+1]
            self.gMin.DefineParameter(17+i*3+1,parName,yyy, self.steps[17+3*i+1],yyy-2, yyy+2)
            parName='objM'+str(i)
            MMM=betterPars[17+3*i+2]
            print MMM,self.steps[17+3*i+2],max(MMM*3,1000)
            self.gMin.DefineParameter(17+i*3+2,parName,MMM, self.steps[17+3*i+2],0.0,max(MMM*3,1000))

        print "Re-Running Minimization!!!"
        self.runMinVar=1
        self.gMin.Migrad()


    def getMin(self,getWithoutMin=0):
        if self.runMinVar or getWithoutMin:
            good=0
            while not good:
                pBestFit=[]
                pErrors=[]
                for i in range(self.numPars):
                    pBestFit.append(ROOT.Double(1.))
                    pErrors.append(ROOT.Double(1.))
                    self.gMin.GetParameter(i,pBestFit[i],pErrors[i])
                    
                    print i,pBestFit[i],pErrors[i]
                (L,dof)=fcnEval(pBestFit)

                #make sure you didn't stumble on a better fit that Minuit didn't use
                if (L-target.bestChi)>0.1:
                    self.CheckWesMin(target.bestPars)
                else:
                    good=1

            print 'Best Fit Chi^2 is '+str(L)
            return (pBestFit,pErrors,L,dof)
        else:
            return (-1,-1)

    def getErrors(self,parsIn=[],whichPars=[9],step=0.01,minimum=False,maximum=False):
        go=0
        if len(parsIn)>0:   ###if the parameters are provided, take them. Otherwise run the fitting mechanism.
            parsG=parsIn[:]
            (lG,dof)=fcnEval(parsG)
            go=1
        elif self.runMinVar:
            (parsG,errsG,lG)=getMin()
            go=1

        print "Current Minimum is: ",lG

        if go:  #do for the brightnesses
            errs=[]
            for ii in whichPars:   #loop over each M provided
                origM=parsG[ii]

                done=0
                diff=0.0
                LarrUp=[]
                diffarrUp=[]
                while not done:
                    diff+=step
                    #plus 
                    parsG[ii]=origM+diff
                    (pL,dof)=fcnEval(parsG)

                    LarrUp.append(pL)
                    diffarrUp.append(diff)
                    if LarrUp[len(LarrUp)-1]>=lG+1 or (maximum<>False and diff>maximum): #cut when at the right chi-squared, or beyond limits
                        done=1

                done=0
                diff=0.0
                LarrDown=[]
                diffarrDown=[]
                while not done:
                    diff+=step
                    #minus
                    parsG[ii]=origM-diff
                    (mL,dof)=fcnEval(parsG)
                    LarrDown.append(mL)
                    diffarrDown.append(diff)
                    if LarrDown[len(LarrDown)-1]>=lG+1 or (minimum<>False and diff>minimum):#cut when at the right chi-squared, or beyond limits
                        done=1
                print lG
                bbb=len(LarrUp)-1
                print LarrUp[bbb],diffarrUp[bbb]
                bbbb=len(LarrDown)-1
                print LarrDown[bbbb],diffarrDown[bbbb]
                
                parsG[ii]=origM

                errs.append([diffarrUp[bbb],-diffarrDown[bbbb]])
            return errs


        else:
            return

    def getPhotometry(self,getWithoutMin=False):
        if not (self.runMinVar or getWithoutMin):
            print "Haven't produced a minimum yet, nor have you stated that this is OK!"
            return None
        else:
            return self.target.getPhotometry()
            
    def createImages(self,doWithoutMin=False):
        if not (self.runMinVar or doWithoutMin):
            print "Haven't produced a minimum yet, nor have you stated that this is OK!"
            return None
        else:
            self.target.produceImage(write='y')
            modelName='model_%s'%(self.target.filename)
            psfName='modelpsf_%s'%(self.target.filename)
            diffName='diff_%s'%(self.target.filename)
            os.system('mv genImage.fits %s'%(modelName))
            os.system('mv genImagePSF.fits %s'%(psfName))
            self.target.imageDiff(name=diffName,write='y')

            return (modelName, psfName, diffName)
