#!/usr/bin/python

from __future__ import print_function

import matplotlib
matplotlib.use('Agg')
from astropy.table import Table,vstack
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import sys
import os
import glob
from subim import extract_subim
from separation import separation
from image_utils import find_bbox
from overlay import show_overlay
import subprocess
import montage_wrapper
from shapely.geometry import Polygon

scale=1.0 # scaling factor for region sizes

def ellipse_poly(x0,y0,a,b,pa,n=200):
    theta=np.linspace(0,2*np.pi,n,endpoint=False)
    st = np.sin(theta)
    ct = np.cos(theta)
    pa = np.deg2rad(pa+90)
    sa = np.sin(pa)
    ca = np.cos(pa)
    p = np.empty((n, 2))
    p[:, 0] = x0 + a * ca * ct - b * sa * st
    p[:, 1] = y0 + a * sa * ct + b * ca * st
    return Polygon(p)

def ellipse(r,ra,dec):
    return ellipse_poly((ra-r['RA'])*np.cos(np.pi*dec/180.0)*3600,(r['DEC']-dec)*3600,r['Maj']*3600,r['Min']*3600,r['PA'])

def intersects(r1,r2):
    return r1['ellipse'].intersects(r2['ellipse'])

def mogrify(filename):
    command='mogrify -quality 90 -trim '+filename
    p = subprocess.Popen(command,shell=True)
    return p


if __name__=='__main__':

    tname=sys.argv[1]
    t=Table.read(tname) # filtered table, i.e. t=ot[ot['Flag_final']==1]
    print('Processing',len(t),'images')

    ot=Table.read('Flowchart_sources_XMMLSS-HighRes-DoubleDetect-2023-05-04_fixed.fits')
    # large source table for neighbours
    lt=ot[(ot['DC_Maj']>8/3600.0)]
    lt=lt[lt['Total_flux']>0.004]
    cname='Source_Name'
    # read in the big files that have all the data
    print('Reading data...')
    gals=Table.read('/beegfs/general/halec001/MIGHTEE/Catalogues/Ks-band/XMM-LSS/XMMFULL_DR3_UNMASKED_Ks_2022_11_11_2.0as_IRAC_2.8as_ALL_CH_prefilter.fits')
    gals['RA'].name='ra'
    gals['DEC'].name='dec'
    lofarfile=fits.open('/beegfs/general/halec001/MIGHTEE/Images/MeerKAT/XMM-LSS/MIGHTEE_Continuum_DR1_XMMLSS_5p0arcsec_I_v1.fits')
    spitzerfiles=[]
    spitzerwcs=[]
    ibandfiles=[]
    ibandwcs=[]
    for i in range(1,4):
        print(i)
        spitzerfiles.append(fits.open('/beegfs/general/halec001/MIGHTEE/Images/3p6um-band/XMM-LSS/xmm%i_SERVS_ch1.fits' % i))
        spitzerwcs.append(WCS(spitzerfiles[-1][0].header))
        ibandfiles.append(fits.open('/beegfs/general/halec001/MIGHTEE/Images/i-band/XMM-LSS/XMM%i_HSC-I_DR3.fits' % i))
        ibandwcs.append(WCS(ibandfiles[-1][0].header))
    marker_ra=None
    marker_dec=None
    
    title=None

    start=int(sys.argv[2])
    try:
        end=int(sys.argv[3])+1
    except:
        end=start+1
    if end>len(t):
        end=len(t)
        
    for i in range(start,end):

        r=t[i]
        print(r)
        sourcename=r['Source_Name']

        iimage=sourcename+'_I.png'
        ipimage=sourcename+'_Ip.png'
        simage=sourcename+'_S.png'
        spimage=sourcename+'_Sp.png'
        manifestname=sourcename+'-manifest.txt'

        if os.path.isfile(manifestname):
            print('Selected output file exists already')
            continue

        ra,dec=r['RA'],r['DEC']

         # new part of the algorithm: first look for any ellipses that
        # intersect the current one(s) (using ot)

        iter=0
        checktable=Table(r)
        checktable['ellipse']=None
        checktable['ellipse'][0]=ellipse(r,ra,dec)
        ctlen=1
        while True:
            tcopy=ot
            tcopy['dist']=np.sqrt((np.cos(dec*np.pi/180.0)*(tcopy['RA']-ra))**2.0+(tcopy['DEC']-dec)**2.0)*3600.0
            filt=tcopy['dist']<180
            filt&=tcopy['dist']>0
            tcopy=tcopy[filt]
            if len(tcopy)==0:
                print('No neighbours found!')
                break
            tcopy['ellipse']=None
            for j in range(len(tcopy)):
                tcopy[j]['ellipse']=ellipse(tcopy[j],ra,dec)
            intersect=[]
            for row in checktable:
                filt=[]
                for row2 in tcopy:
                    if row[cname]==row2[cname]:
                        result=False
                    else:
                        result=intersects(row,row2)
                    print('Check intercept between',row[cname],'and',row2[cname],'result is',result)
                    
                    filt.append(result)
            print('Filt is',filt)
            itable=tcopy[filt]
            #print 'Length of checktable is',len(checktable)
            #print 'Length of itable is',len(itable)
            interim=vstack((checktable,itable))
            #print 'Length of interim table is',len(interim)
            names=list(set(interim[cname]))
            #print 'Names are',names
            filt=[False]*len(interim)
            for n in names:
                pos=np.argmax(interim[cname]==n)
                filt[pos]=True
            checktable=interim[filt]
            print('New checktable length is',len(checktable))
            if len(checktable)==ctlen:
                print('Checktable length converged!')
                break
            ctlen=len(checktable)
            iter+=1
            if iter==10:
                print('Not converged!')
                break

        #print checktable
        '''
        import matplotlib.pyplot as plt
        for r in checktable:
            x,y=r['ellipse'].exterior.xy
            plt.plot(x,y)
        plt.savefig('poly.pdf')
        '''
        #ra,dec,size=find_bbox(checktable)
        
        # resize the image to look for interesting neighbours (using lt)
        iter=0
        flux=r['Total_flux']
        while True:
            startra,startdec=ra,dec
            tcopy=lt
            tcopy['dist']=np.sqrt((np.cos(dec*np.pi/180.0)*(tcopy['RA']-ra))**2.0+(tcopy['DEC']-dec)**2.0)*3600.0
            filt=tcopy['dist']<60
            filt&=tcopy['Total_flux']>flux/3.0 # look for similar flux
            tcopy=tcopy[filt]
            print('Iter',iter,'found',len(tcopy),'neighbours')

            # include checktable
            interim=vstack((tcopy,checktable))
            # de-dupe again
            names=list(set(interim[cname]))
            print('Names are',names)
            filt=[False]*len(interim)
            for n in names:
                pos=np.argmax(interim[cname]==n)
                filt[pos]=True
            tcopy=interim[filt]
            
            ra=np.mean(tcopy['RA'])
            dec=np.mean(tcopy['DEC'])
            
            newra,newdec,size=find_bbox(tcopy,scale=1)
            print('     size is',size)

            if startra==ra and startdec==dec:
                break
            iter+=1
            if iter==10:
                break

        # now find the bounding box of the resulting collection
        ra,dec,size=find_bbox(tcopy,scale=1)

        if np.isnan(size):
            ra=r['RA']
            dec=r['DEC']
            size=60

        if size>180.0:
            # revert just to original
            ra,dec=r['RA'],r['DEC']
            tcopy=checktable
            ra,dec,size=find_bbox(tcopy,scale=1)

        if size>180:
            size=180.0
        if size<60:
            size=60.0
        size=(int(0.5+size/10))*10
        print('final size is',size)
        print('final RA DEC is',ra,dec)
 
        size/=3600.0

        seps=separation(ra,dec,ot['RA'],ot['DEC'])
        ots=ot[seps<size*2]


        ls=[]
        for nr in ots:
            if nr['Source_Name']==r['Source_Name']:
                ls.append('solid')
            else:
                ls.append('dashed')

        rms=r['Isl_rms']
       
        lhdu=extract_subim(lofarfile,ra,dec,size)

        # match to correct files
        for i in range(3):
            ysize,xsize=np.shape(spitzerfiles[i][0].data)
            x,y=spitzerwcs[i].wcs_world2pix(ra,dec,0)
            if x>0 and y>0 and x<xsize and y<ysize:
                spitzerfile=spitzerfiles[i]
                break
        else:
            raise RuntimeError('Cannot find Spitzer file')

        for i in range(3):
            ysize,xsize=np.shape(ibandfiles[i][0].data)
            x,y=ibandwcs[i].wcs_world2pix(ra,dec,0)
            if x>0 and y>0 and x<xsize and y<ysize:
                ibandfile=ibandfiles[i]
                break
        else:
            raise RuntimeError('Cannot find I-band file')
        
        try:
            shdu=extract_subim(spitzerfile,ra,dec,size)
        except RuntimeError:
            continue
        ihdu=extract_subim(ibandfile,ra,dec,size)
        lhdu.writeto(sourcename+'_m.fits',overwrite=True)
        ihdu.writeto(sourcename+'_i.fits',overwrite=True)
        shdu.writeto(sourcename+'_s.fits',overwrite=True)
        #montage_wrapper.mGetHdr(sourcename+'_i.fits',sourcename+'_i.hdr') 
        #montage_wrapper.mProject(sourcename+'_s.fits',sourcename+'_so.fits',sourcename+'_i.hdr')
        #ihdu[0].data=np.where(ihdu[0].data>49999,np.nan,ihdu[0].data)
        #shdu=fits.open(sourcename+'_so.fits')

        pg=gals[(np.abs(gals['ra']-ra)<(size/np.cos(dec*np.pi/180.0))) & (np.abs(gals['dec']-dec)<size)]

        try:
            peak=r['Peak_flux']
        except:
            peak=None
        print('Peak flux of the object is',peak)
        plist=[]
        try:
            if not os.path.isfile(iimage):
                show_overlay(lhdu,ihdu,ra,dec,size,firsthdu=None,overlay_cat=ots,overlay_scale=scale,coords_color='limegreen',contour_color='white',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=5,lw=2,save_name=iimage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,ellipse_color='cyan',marker_color='cyan',title=title,lw_ellipse=3,ellipse_style=ls,minimum=peak/2,noisethresh=1)
                plist.append(mogrify(iimage))
            if not os.path.isfile(ipimage):
                show_overlay(lhdu,ihdu,ra,dec,size,firsthdu=None,overlay_cat=ots,overlay_scale=scale,coords_color='limegreen',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=5,lw=2,save_name=ipimage,show_lofar=False,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,lw_ellipse=3,ellipse_style=ls,ellipse_color='cyan',minimum=peak/2,noisethresh=1,plotpos=pg,ppsize=350)
                plist.append(mogrify(ipimage))
            if not os.path.isfile(simage):
                show_overlay(lhdu,shdu,ra,dec,size,firsthdu=None,overlay_cat=ots,overlay_scale=scale,coords_color='limegreen',contour_color='white',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=5,lw=2,save_name=simage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,lw_ellipse=3,ellipse_style=ls,ellipse_color='cyan',minimum=peak/2)
                plist.append(mogrify(simage))
            if not os.path.isfile(spimage):
                show_overlay(lhdu,shdu,ra,dec,size,firsthdu=None,overlay_cat=ots,overlay_scale=scale,coords_color='limegreen',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=5,lw=2,save_name=spimage,show_lofar=False,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',ellipse_color='cyan',title=title,lw_ellipse=3,ellipse_style=ls,minimum=peak/2)
                plist.append(mogrify(spimage))

            with open(manifestname,'w') as manifest:
                manifest.write('%i,%s,%s,%s,%s,%s,%f,%f,%f\n' % (i,iimage,ipimage,simage,spimage,sourcename,ra,dec,size*3600.0))

            for p in plist:
                p.wait()
        except Exception as e:
            with open(manifestname.replace('manifest','error'),'w') as outfile:
                outfile.write(str(e))

            
                         
