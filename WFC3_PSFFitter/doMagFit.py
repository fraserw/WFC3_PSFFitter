#! /usr/bin/env python

import WFC3_PSFFits as P
import WFC3_phot as phot
import sys
import MySQLdb
import os
import sourceFind
import numpy as num

def getR(image,Z):
    w=num.where(image==num.max(image))

    peakY=w[0][0]
    peakX=w[1][0]

    w=num.where(image[peakY-6:peakY+7,peakX-6:peakX+7]>Z+num.sqrt(Z)*2)

    dx=(num.max(w[0])-num.min(w[0]))/2.
    dy=(num.max(w[1])-num.min(w[1]))/2.
    return (dx+dy)/2.



conn = MySQLdb.connect(host = 'localhost',
                         user = 'fraserw',
                         passwd = 'ugly_411',
                         db = 'HST')
curs = conn.cursor()


han=open('centroids_17.out')
centData=han.readlines()
han.close()
for i in centData:
    s=i.split()
    obj=s[0]
    I=int(float(s[1]))
    fn=s[2].split('/')[1]
    x=float(s[3])
    y=float(s[4])


    #make sure the table entry exists
    curs.execute("""SELECT dir from fits17_CalWF31 WHERE dir=%s""",(obj,))
    result=curs.fetchone()
    if result==None:
        curs.execute("""INSERT INTO fits17_CalWF31 (dir) VALUES (%s)""",(obj,))
        conn.commit()

    os.chdir(obj)


    #get the image background
    outObj=P.PSFFit(fn,1,useMemory='n')
    (z,med)=P.getBG(outObj.realImData,int(x),int(y),60)
    print 'Background:',z,med
    print

    

    #choose the box size for chi fitting
    if I<4: 
        r=4
        f=2000.
    else: 
        r=2
        f=30.




    #initialize the fit and run it
    minuit=P.Min4PSFFit(fn,1,[[0,x,y,f,3.5]],useMemory='y',psfSubSampFactor=105)
    minuit.initMin(minx=x-r,maxx=x+r+1,miny=y-r,maxy=y+r+1,Z=z,wx=0.0,wy=0.0,rot=0.0,fitMask=[],fixFocus=True)
    minuit.fixPar(17) #fix x centroid
    minuit.fixPar(18) #fix y centroid
    minuit.runMin()
                
    #get the best-fit parameters. Since we did the centroid elsewhere, all we care about here is f.
    (pars,errs,chi,dof)=minuit.getMin()

    f=pars[len(pars)-1]
    
    #get the psf matching photometry
    (flux,m_fit)=minuit.getPhotometry()
    m_fit=m_fit[0]
    
    #create the PSF image, the model image (PSF+background), and the difference image
    (model,psf,diff)=minuit.createImages()
    minuit.gMin.mncler()
    os.chdir('..')


    #now get the image data and do the psf photo centroid

    psfImage=minuit.target.image['image']
    #image=minuit.target.realImageHan[1].data

    #dr=getR(psfImage,z)
    #if I<4:
    #    dr=min(dr,4.)
    #else:
    #    dr=min(dr,2.)   
    #(cent_x,cent_y)=sourceFind.centroid(image,x,y,dr=dr,iterations=100,tol=0.01,repeatFactor=50)
    #print 'Fraser centroid using aperture radius %s: %s %s'%(dr,cent_x,cent_y)
    
    #finally get the aperture photometry
    #use the model psf image to provide the best aperture centroid
    rad=3.
    (a,b,c,d,e)=phot.phot(obj+'/'+fn,[[x,y]],useOtherData=psfImage,radii=[rad],cent='centroid',psfInterp=False)
    cent_x=d[0]
    cent_y=e[0]
    print 'Iraf centroid:',cent_x,cent_y

    (a,b,c)=phot.phot(obj+'/'+fn,[[cent_x,cent_y]],radii=[rad],cent='none',psfInterp=True)
    m_aper=a[0][0]

    redo=False
    if m_aper>24 and (i==4 or i==6):
        rad=3.5
        redo=True
    elif m_aper>24 and (i==5 or i==7):
        rad=3.5
        redo=True
           
    if redo:
        print 'Redoing photometry with a smaller magnitude.'
        (a,b,c)=phot.phot(obj+'/'+fn,[[x,y]],radii=[rad],cent='none',psfInterp=True)
        m_aper=a[0][0]

    print x,y
    print cent_x,cent_y
    print m_aper,m_fit

    tuple=(I, fn, I, x, I, y, I, f, I, m_fit, I, chi, I, cent_x, I, cent_y, I, m_aper, I, rad, obj)
    command="""UPDATE fits17_CalWF31 SET fn%s=%s, x%s=%s, y%s=%s, m%s=%s, fitmag%s=%s, chi%s=%s, cent_x%s=%s, cent_y%s=%s, aper_mag%s=%s, rad%s=%s where dir=%s"""
    print
    print "Updating database:"
    print command%tuple
    print

    curs.execute(command,tuple)
    conn.commit()
    

curs.close()
conn.close()


"""
#create table command
#CREATE TABLE fits17_CalWF31 (i int NOT NULL AUTO_INCREMENT, dir VARCHAR(15), fn0 VARCHAR(15), x0 float DEFAULT -32768., y0 float DEFAULT -32768, m0 float DEFAULT -32768., fitmag0 float DEFAULT -32768., chi0 float DEFAULT -32768., cent_x0 float DEFAULT -32768.,cent_y0 float DEFAULT -32768., aper_mag0 float DEFAULT -32768., rad0 float DEFAULT -32768., aCorr0 float DEFAULT -32768., fn1 VARCHAR(15), x1 float DEFAULT -32768., y1 float DEFAULT -32768, m1 float DEFAULT -32768., fitmag1 float DEFAULT -32768., chi1 float DEFAULT -32768., cent_x1 float DEFAULT -32768.,cent_y1 float DEFAULT -32768., aper_mag1 float DEFAULT -32768., rad1 float DEFAULT -32768., aCorr1 float DEFAULT -32768., fn2 VARCHAR(15), x2 float DEFAULT -32768., y2 float DEFAULT -32768, m2 float DEFAULT -32768., fitmag2 float DEFAULT -32768., chi2 float DEFAULT -32768., cent_x2 float DEFAULT -32768.,cent_y2 float DEFAULT -32768., aper_mag2 float DEFAULT -32768., rad2 float DEFAULT -32768., aCorr2 float DEFAULT -32768., fn3 VARCHAR(15), x3 float DEFAULT -32768., y3 float DEFAULT -32768, m3 float DEFAULT -32768., fitmag3 float DEFAULT -32768., chi3 float DEFAULT -32768., cent_x3 float DEFAULT -32768.,cent_y3 float DEFAULT -32768., aper_mag3 float DEFAULT -32768., rad3 float DEFAULT -32768., aCorr3 float DEFAULT -32768., fn4 VARCHAR(15), x4 float DEFAULT -32768., y4 float DEFAULT -32768, m4 float DEFAULT -32768., fitmag4 float DEFAULT -32768., chi4 float DEFAULT -32768., cent_x4 float DEFAULT -32768.,cent_y4 float DEFAULT -32768., aper_mag4 float DEFAULT -32768., rad4 float DEFAULT -32768., aCorr4 float DEFAULT -32768., fn5 VARCHAR(15), x5 float DEFAULT -32768., y5 float DEFAULT -32768, m5 float DEFAULT -32768., fitmag5 float DEFAULT -32768., chi5 float DEFAULT -32768., cent_x5 float DEFAULT -32768.,cent_y5 float DEFAULT -32768., aper_mag5 float DEFAULT -32768., rad5 float DEFAULT -32768., aCorr5 float DEFAULT -32768., fn6 VARCHAR(15), x6 float DEFAULT -32768., y6 float DEFAULT -32768, m6 float DEFAULT -32768., fitmag6 float DEFAULT -32768., chi6 float DEFAULT -32768., cent_x6 float DEFAULT -32768.,cent_y6 float DEFAULT -32768., aper_mag6 float DEFAULT -32768., rad6 float DEFAULT -32768., aCorr6 float DEFAULT -32768., fn7 VARCHAR(15), x7 float DEFAULT -32768., y7 float DEFAULT -32768, m7 float DEFAULT -32768., fitmag7 float DEFAULT -32768., chi7 float DEFAULT -32768., cent_x7 float DEFAULT -32768.,cent_y7 float DEFAULT -32768., aper_mag7 float DEFAULT -32768., rad7 float DEFAULT -32768., aCorr7 float DEFAULT -32768., H float DEFAULT -32768., dH float DEFAULT -32768., primary key(i)) AUTO_INCREMENT=0;
"""

