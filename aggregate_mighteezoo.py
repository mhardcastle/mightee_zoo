# Code to analyse LGZ click data to generate catalogues.
# Input is four csv files generated by galzoo_export.py
# Output is a source catalogue containing association and
# optical ID information, and a component catalogue listing
# basic info for all source components.
#
# Written by J. Croston and modified by M.J. Hardcastle
#
# this version for MIGHTEE Zoo

from __future__ import print_function

import sys
import pandas as pd
import numpy as np
import time
import os

from astropy.table import Table
from itertools import combinations
from collections import defaultdict
from astropy.coordinates import SkyCoord
import astropy.units as u

from scipy.interpolate import interp1d

pixim=866.0
rpix=433.0
dpix=433.0

def getname(srcra,srcdec):

    # Determine final source name from coordinates
    
    sc=SkyCoord(srcra,srcdec,unit='deg',frame='icrs')
    s=sc.to_string(style='hmsdms',sep='',precision=2)
    ilt=str('MGTC'+s).replace(' ','')[:-1]
    return ilt
        

def pointInEllipse(x,y,xp,yp,d,D,angle):
    
    #tests if a point[xp,yp] is within
    #boundaries defined by the ellipse
    #of center[x,y], diameter d D, and tilted at angle

    cosa=np.cos(angle*np.pi/180.0)
    sina=np.sin(angle*np.pi/180.0)
    dd=(d/2.0)*(d/2.0)
    DD=(D/2.0)*(D/2.0)

    a = (cosa*(xp-x)+sina*(yp-y))**2.0
    b = (sina*(xp-x)-cosa*(yp-y))**2.0
    ellipse=(a/dd)+(b/DD)

    if ellipse <= 1:
        return True
    else:
        return False

def match_click_ellipse(cat,click,nsource):

    # Matches click positions to component ellipses from catalogue cat
    
    x=click['x']
    y=click['y']
    user=click['user_name']
    classification_id=click['classification_id']
    # Convert coords to ra, dec
    rx = ra+(x-rpix)*(-ascale)/(3600.0*np.cos(dec*np.pi/180.0))
    dy = dec+(dpix-y)*ascale/3600.0

    matches=[]
    smatches=[]
    lista=[]
    classesa=[]
    assocpos=[]

    for s in cat:
        # Ignore target source and loop through others with centres in reg
        if s['Source_Name'] != nsource:
            # These are the possible matches
            xcheck=rpix+(s['RA']-ra)*np.cos(dec*np.pi/180.0)*3600.0*(-pscale)
            ycheck=dpix-(s['DEC']-dec)*3600.0*pscale
            realy=pixim-ycheck

            if(pointInEllipse(xcheck,realy,x,pixim-y,2.0*pscale*s['Maj'],2.0*pscale*s['Min'],s['PA']+90.0) is True):
                matches.append(s)
                        
    if len(matches)>0:
        # At least one match, so pick one where closest to centre
        mdist=10000.0
        for m in matches:
            val = np.sqrt((np.abs(rx-m['RA'])*np.cos(dec*np.pi/180.0))**2.0+(np.abs(dy-m['DEC'])**2.0))
            if val < mdist:
                mdist=val
                goodmatch=m
                # Pick best match
                
                #smatches.append(goodmatch['Source_Name'])
        asrc=goodmatch['Source_Name']
        newclist=assocdict.get(asrc,[])
        print("In assoc code, classification ID is "+str(classification_id))
        if classification_id not in newclist:
            newclist.append(classification_id)
            assocdict.update({asrc:newclist})
            lista.append({'Assoc_Name':goodmatch['Source_Name'],'Source_Name':nsource,'User':user, 'Class':classification_id})
            assocpos.append((goodmatch['Source_Name'],goodmatch['RA'],goodmatch['DEC'],goodmatch['Total_flux'],goodmatch['Maj'],goodmatch['Min'],goodmatch['PA']))
            classesa.append(classification_id)
    else:
        matched_orig=False
        # did the user click in the solid ellipse?
        for s in cat:
            if s['Source_Name'] == nsource:
                # These are the possible matches
                xcheck=rpix+(s['RA']-ra)*np.cos(dec*np.pi/180.0)*3600.0*(-pscale)
                ycheck=dpix-(s['DEC']-dec)*3600.0*pscale
                realy=pixim-ycheck

                if(pointInEllipse(xcheck,realy,x,pixim-y,2.0*pscale*s['Maj'],2.0*pscale*s['Min'],s['PA']+90.0) is True):
                    matched_orig=True
        if not matched_orig:
            print('Warning: click not matched to any ellipse',x,y,rx,dy)
    return classesa,lista,assocpos

def match_source_cat(cat,mra,mdec):

    # Matches components from different source catalogues because argh.
    
    rx = mra
    dy = mdec

    match=None
    # Use 3 arcsec as minimum matching dist

    nflux=None
    nmaj=None
    nmin=None
    nposang=None
    
    if len(cat)>0:
        mdist=0.00083
        for s in cat:
            val = np.sqrt(((rx-s['RA'])*np.cos(dy*np.pi/180.0))**2.0+(dy-s['DEC'])**2.0)
            if val < mdist:
                mdist=val
                match=s['Source_Name']
                nflux=s['Total_flux']
                nmaj=s['Maj']
                nmin=s['Min']
                nposang=s['PA']
        
    nsname=match
        
    return nsname,nflux,nmaj,nmin,nposang


def match_opt(click,cat):
    # Matches optical click positions to merged optical catalogue
    
    x=click['x']
    y=click['y']
    user=click['user_name']

    rx = ra+(x-rpix)*(-ascale)/(3600.0*np.cos(dec*np.pi/180.0))
    dy = dec+(dpix-y)*ascale/3600.0
  
    poss="None"
    pra=np.nan
    pdec=np.nan
    dist=3600.0*np.sqrt((np.abs(rx-cat[rast])*np.cos(dec*np.pi/180.0))**2.0+(np.abs(dy-cat[decst])**2.0))
    if np.any(dist<1.5): ## fixed association radius
        i=np.argmin(dist)
        poss=str(cat[i][idname])
        pra=cat[i][rast]
        pdec=cat[i][decst]
        thres=dist[i]

    return poss,pra,pdec,user
                
def write_csv_out(col_order,data_list,outfile,label):

    outpd=pd.DataFrame(data_list)[col_order]
    outpd.to_csv(outfile,index_label=label)

def write_fits_out(col_order,data_list,outfile):

    outpd=pd.DataFrame(data_list)[col_order]
    fitst=Table.from_pandas(outpd)
    fitst.write(outfile,overwrite=True)
    
def tally_votes(matches,numseen):

     # For each unique ID match, sum the number of votes and compare to the number of times the
     # source was seen to determine the quality flag. 

    
    uniquematches=set(matches)
    goods=[]
    for m in uniquematches:
            numioccur=matches.count(m)
     
            frac=float(numioccur)/float(numseen)

            if frac > 0.99:
                aflag=1                
            elif frac>=2.0/3.0:
                aflag=2
            else:
                aflag=3
                
            if aflag<=2:
                # Good match
                goods.append((m,frac))


    return goods

# Main code

# We assume this is running in the directory where galzoo_export has run

#field_l=sys.argv[1]
#field_h=sys.argv[2]

assocfile='lofar-associations.csv'
idfile='lofar-opticalids.csv'
subjfile='lofar-allclassifications.csv'
probfile='lofar-classification-problems.csv'

# Required input is the four output files produced by galzoo_export.py
    
assocs = pd.read_csv(assocfile)
ids = pd.read_csv(idfile)
subjects = pd.read_csv(subjfile)
problems = pd.read_csv(probfile)

#sourcecat = "/data/lofar/mjh/dr2/0h/lotss_dr2_ra0inner_masked_v0.9.srl.fits"
#sourcecat='/data/lofar/DR2/catalogues/LoTSS_DR2_v100.srl.fits'
#optcat='/beegfs/lofar/mjh/dr2/dr2_combined.fits' #  % (field_l,field_h)

wd=os.getcwd()

if 'cosmos' in wd:
    sourcecat='../cosmos_new/Flowchart_sources_COSMOS-HighRes-DoubleDetect-2023-05-04_fixed.fits'
    optcat='/beegfs/general/halec001/MIGHTEE/Catalogues/Ks-band/COSMOS/COSMOSFULL_DR3_UNMASKED_Ks_2023_03_11_2.0as_IRAC_2.8as_ALL_CH_prefilter.fits'
elif 'cdfs' in wd:
    sourcecat='../cdfs/Flowchart_sources_CDFS-DEEP-HighRes-DoubleDetect-2023-05-04_fixed.fits'
    optcat='/beegfs/general/halec001/MIGHTEE/Catalogues/Ks-band/CDFS/CDFSFULL_DR3_UNMASKED_Ks_2022_11_30_2.0as_IRAC_2.8as_ALL_CH_prefilter.fits'
elif 'xmm' in wd:
    sourcecat='../xmm-lss/Flowchart_sources_XMMLSS-HighRes-DoubleDetect-2023-05-04_fixed.fits'
    optcat='/beegfs/general/halec001/MIGHTEE/Catalogues/Ks-band/XMM-LSS/XMMFULL_DR3_UNMASKED_Ks_2022_11_11_2.0as_IRAC_2.8as_ALL_CH_prefilter.fits'
opt=Table.read(optcat)
idname='ID'
#opt['ID']=np.where(opt['UNWISE_OBJID']!="N/A",opt['UNWISE_OBJID'],opt['UID_L'])
rast='RA'
decst='DEC'


# Generate unique list of subjects with either associations, IDs or both
# This is done by source_name, because the subject IDs are not unique between batches

newsubj=subjects.drop_duplicates(subset='source_name')
subs_unique=set(newsubj.source_name)

# Read in catalogues

t=Table.read(sourcecat)
srct=t
# Convert axes to arcsec, as the code expects
t['Maj']*=3600.0
t['Min']*=3600.0

# Set up a bunch of lists

badassocs=[]
subjs_noassoc=[]
subjs_noids=[]
alist=[]
tlist=[]
ilist=[]
multlist=[]
multsets=[]
subjs_mult=[]
artlist=[]
blendlist=[]
zoomlist=[]
brokenlist=[]
missinglist=[]
finalcat=[]
assoccat=[]
compcat=[]
classa=[]
sources_unique=[]
classid=[]
catalogued=[]
overlap=[]
seendict={}
fluxdict={}
sizedict={}
mindict={}
padict={}
coordsdict=defaultdict(list)
ocoordsdict=defaultdict(list)
probdict=defaultdict(list)
optiddict=defaultdict(list)
assocdict=defaultdict(list)
badclickdict=defaultdict(int)
debug=True

i=0

# Loop through all sources by source_name to match clicks with associated regions and optical IDs

new_subs_unique=[]
nonexist=[]
nct=0
for scount,subj in enumerate(subs_unique):

    if debug:
        print('subj is %s (%i of %i)' % (subj,scount,len(subs_unique)))
    probflag=0
    assocs_s = assocs.query('source_name == @subj')
    ids_s = ids.query('source_name == @subj')
    probs_s = problems.query('source_name == @subj')
    subj_details=newsubj.query('source_name == @subj')
    subj_class=subjects.query('source_name == @subj')
    num_assoc = len(assocs_s)
    num_ids = len(ids_s)
    num_subj = len(subj_details)
    num_probs = len(probs_s)
    source_name = subj

    ra = subj_details.iloc[0]['ra']     # ra at pixel 433
    dec = subj_details.iloc[0]['dec']   # dec at pixel 433
    size = subj_details.iloc[0]['size']   # image size in arcsec
    pscale=pixim/size
    ascale=size/pixim    
    
    sourcematches=[]
    imatches=[]

    tabval=1
    rowt=srct[srct['Source_Name']==subj]
    if len(rowt)>0:
        flux=rowt[0]['Total_flux']
        majsize=rowt[0]['Maj']
        minsize=rowt[0]['Min']
        posang=rowt[0]['PA']
        sourcera=rowt[0]['RA']
        sourcedec=rowt[0]['DEC']
    
    new_subs_unique.append(source_name)

    coordsdict.update({source_name:(sourcera,sourcedec)})
    fluxdict.update({source_name:flux})
    sizedict.update({source_name:majsize})
    mindict.update({source_name:minsize})
    padict.update({source_name:posang})
    
    # Now set up dict of how many classifiers have seen this source as main subject
    # This does not include sights as a potential association, which are dealt with later

    users=subj_class.user_name
    numseen=len(subj_class)
    seendict.update({source_name:numseen})
    
    # Also set up a dict of which catalogue each source is in for future ref

    # Identify problems and tally votes for problems
    
    if num_probs>0:
        probflag=1
        acount=0
        bcount=0
        zcount=0
        brcount=0
        mcount=0
        for ind,p in probs_s.iterrows():
            if p['problem']=="Artefact":
                acount=acount+1
            if p['problem']=="Blend":
                bcount=bcount+1
            if p['problem']=="Too zoomed in" or p['problem']=="Too zoomed in ":
                zcount=zcount+1
            if p['problem']=="Host galaxy broken up":
                brcount=brcount+1
            if p['problem']=="Image missing":
                mcount=mcount+1
    
        artfrac=float(acount)/float(numseen)
        blendfrac=float(bcount)/float(numseen)
        zoomfrac=float(zcount)/float(numseen)
        brokefrac=float(brcount)/float(numseen)
        missfrac=float(mcount)/float(numseen)

        probdict.update({source_name:(artfrac,blendfrac,zoomfrac,brokefrac,missfrac)})    
            
    else:
        probdict.update({source_name:(0.0,0.0,0.0,0.0,0.0)}) 
    numgood=0
    numigood=0

    # Loop through association clicks for this source to match each click to catalogued radio source
    
    if num_assoc>0:

        rsize=2.0*size/3600.0
        dsize=2.0*size/3600.0
        subcat=srct[(np.abs(srct['RA']-ra)*np.cos(dec*np.pi/180.0)<rsize) & (np.abs(srct['DEC']-dec)<dsize)]
        for index,news in assocs_s.iterrows():
            classesa,lista,assocpos=match_click_ellipse(subcat,news,source_name)
            for thing in classesa:
                classa.append(thing)
            for thing in lista:
                alist.append(thing)
            for thing in assocpos:
                nam=thing[0]
                nra=thing[1]
                ndec=thing[2]
                nflux=thing[3]
                nmajsiz=thing[4]
                nminsiz=thing[5]
                nposang=thing[6]
                coordsdict.update({nam:(nra,ndec)})
                fluxdict.update({nam:nflux})
                sizedict.update({nam:nmajsiz})
                mindict.update({nam:nminsiz})
                padict.update({nam:nposang})
        del(subcat)
    else:
        # No associations for this subject
        
        subjs_noassoc.append(source_name)


    # Loop through ID clicks for this source to match each click to a catalogued optical source

    if num_ids>0:

        rsize=0.5*size/3600.0
        dsize=0.5*size/3600.0

        # Define sub-catalogues for the relevant sky regions
        
        osubcat=opt[(np.abs(opt[rast]-ra)*np.cos(dec*np.pi/180.0)<rsize) & (np.abs(opt[decst]-dec)<dsize)]
        
         
        for index,news in ids_s.iterrows():

            # Optical match
            poss,ora,odec,user=match_opt(news,osubcat)
            classification_id=news['classification_id']
            ocoordsdict.update({poss:(ora,odec)})
            if poss != "None":
                if debug: print('Matched with',poss)
                # check if ID/classification pair already exists
                classlist=optiddict.get(poss,[])
                if classification_id not in classlist:
                    imatches.append(poss)
                    classlist.append(classification_id)
                    optiddict.update({poss:classlist})
                else:
                    print("Found duplicate PS ID/class pair: "+str(poss)+": "+str(classification_id))
                
                cflag=0
            else:
                badclickdict[subj]+=1
                if debug: print('A click for this source was not matched to the optical catalogue',badclickdict[subj])

        # Now tally up votes to generate a list of source_name, optical_ID name
        # This will later be matched to the associated sources where appropriate
      
        goods=tally_votes(imatches,numseen)
        numigood=len(goods)
        
        if numigood==1:
                       
            ilist.append({'Source_Name':source_name,'Opt_id':goods[0][0],'Numclass':numseen,'Quality':goods[0][1]})

        elif numigood>1:
            for thing in goods:
                multlist.append({'Source_Name':source_name,'Opt_id':thing[0],'Numclass':numseen,'Quality':thing[1]})
            subjs_mult.append(source_name)
            
        else:
            subjs_noids.append(source_name)
        
    else:
# No IDs
        subjs_noids.append(source_name)
            
    i=i+1

# Sort associated source_names into sets an individual clicker has associated together

print('Creating sets')

classinassocs=set(classa)

listofsets=[]
for classed in classinassocs:
    
    # Define a set of associated sources and add to list of sets

    newset=[]
    allassocs=filter(lambda x: x['Class'] == classed, alist)


    for thing in allassocs:
        if thing['Source_Name'] in newset:
            newset.append(thing['Assoc_Name'])
        else:
            newset.append(thing['Assoc_Name'])
            newset.append(thing['Source_Name'])

    listofsets.append(newset)
    

# Create a dictionary that lists which additional subjects an associated subject is seen with.
# The dictionary key is the subject, the value lists all other subjects whose "numseen" value should be
# included in the denominator for vote tallying for that key. 

print('Creating associated sources')

seenwithdict=defaultdict(list)

for assocsrc in alist:
    if assocsrc['Assoc_Name'] in new_subs_unique:
        sourcesubj=assocsrc['Source_Name']
        asubj=assocsrc['Assoc_Name']
        if sourcesubj not in seenwithdict[asubj]:
            seenwithdict[asubj].append(sourcesubj)

# For each main-subject source_name, consider all sets that include it and loop round
# to tally the votes for that set

allgoodsets=[]
setscore=[]

for subj in new_subs_unique:

    goodsets=[]
    subjseen=seendict.get(subj)
    addseen=0
    addassoc=seenwithdict.get(subj)

    if addassoc is not None:
        for asrc in addassoc:
            addseen=addseen+seendict.get(asrc)

    # subjseen is the total tally of how many times source_name was seen EITHER as subject or
    # in the field of another subject
    
    subjseen=subjseen+addseen
        
    subjset=[x for x in listofsets if subj in x]
 
    flat_list = [item for sublist in subjset for item in sublist]
    uniqelems=set(flat_list)
    nelem=len(uniqelems)
    
    # Consider all combinations of the unique sources anyone has suggested may be associated

    c = [comb for i in range(1,nelem) for comb in combinations(uniqelems, i + 1)]

    newc=set(c)

    # Loop through all unique combinations to tally the votes for that combination
    
    for tesc in newc:
        ts=set(tesc)
        vote=0
        for slist in subjset:

            # This combination gets a vote for every listed set of which it is a subset
            # because we want to catch cases where clickers agree a common subset even
            # if full set is not in agreement

            if ts.issubset(slist):
                
                vote=vote+1

        # Calculate the score for this possible combination
                
        score=float(vote)/float(subjseen)
       
        if score>0.667:
            goodsets.append(ts)
            setscore.append((ts,score))


    left=[]

    # Remove any good sets that are either a subset of a larger set, or that don't include the
    # subject -- the latter should be tabulated elsewhere to avoid duplication

    if len(goodsets)>1:
 
        rem=[]
        for item in goodsets:
            for otheritem in goodsets:
                if item.issubset(otheritem) and len(item) < len(otheritem):
                    rem.append(item)
            if subj not in item:
                rem.append(item)
        for item in goodsets:
            if item not in rem:
                left.append(item)
    else:
        rem=[]
        for item in goodsets:
            if subj not in item:
                rem.append(item)
            else:
                left.append(item)
            
        
    numsets=len(left)
    for item in left:
        if item not in allgoodsets:
            allgoodsets.append(item)

# The list "allgoodsets" should now contain the largest set with a score>2.0/3.0 containing each target source, where such a set exists. 
        
ngsets=len(allgoodsets)

print("Total number of associated sources is ",ngsets)

# Now we need to generate the final list of sources, including optical IDs and replacing individual components with associated sets where appropriate 

print('Generating final list')

finalcat=[]

for finalset in allgoodsets:
    oflag=0
    newscore=[x for x in setscore if x[0]==finalset]
    # Identify members that are subjects and find their IDs
    setsubjs=[x for x in finalset if x in new_subs_unique]
    nflux=0
    apb=0
    bpb=0
    zpb=0
    hpb=0
    imb=0
    if len(setsubjs)>0:
        
        for fsubj in setsubjs:
                            
            optids=[]
            findilist=[x for x in ilist if x['Source_Name']==fsubj]
            if len(findilist)>0:
                if len(findilist)>1:
                    print("Subject has multiple entries in ID table -- this shouldn't happen!")
                else:
                    nam=findilist[0]['Opt_id']
                    numseen=findilist[0]['Numclass']
                    qual=findilist[0]['Quality']
                    optids.append((nam,numseen,qual))
            fapb,fbpb,fzpb,fhpb,fimb=probdict.get(fsubj)
            if fapb>apb:
                apb=fapb
            if fbpb>bpb:
                bpb=fbpb
            if fzpb>zpb:
                zpb=fzpb
            if fhpb>hpb:
                hpb=fhpb
            if fimb>imb:
                imb=fimb
                
        if len(optids)==0:
            optid="None"
            optra=np.nan
            optdec=np.nan
            qual=np.nan
        elif len(optids)>1:
            optid="Mult"
            optra=np.nan
            optdec=np.nan
            qual=np.nan
            sets_mult.append(fsubj)
            for thing in optids:
                nam=thing[0]
                numseen=thing[1]
                tqual=thing[2]
                multlist.append({'Source_Name':subj,'Opt_id':nam,'Numclass':numseen,'Quality':tqual})
        else:
            optid=optids[0][0]
            qual=optids[0][2]
            optra,optdec=ocoordsdict.get(nam)
            
        ralist=[]
        declist=[]
        coordslist=[]
        badclicklist=[]
        msiz=0
        for elem in finalset:
            if elem in catalogued:
                overlap.append(elem)
            else:
                catalogued.append(elem)
                
            aflux=fluxdict.get(elem)
            nflux=nflux+aflux
            acoords=coordsdict.get(elem)
            osiz=sizedict.get(elem)
            if osiz>msiz:
                msiz=osiz
            ralist.append(acoords[0])
            declist.append(acoords[1])
            coordslist.append((acoords[0],acoords[1]))
            badclicklist.append(badclickdict[elem])
                            
        ramean=sum(ralist)/float(len(ralist))
        decmean=sum(declist)/float(len(declist))
        print('badclicklist is',badclicklist)
        badclick=sum(badclicklist)
        ndists=[]
        for elem in coordslist:
            for otherelem in coordslist:
                ndists.append(np.sqrt((np.abs(elem[0]-otherelem[0])*np.cos(elem[1]*np.pi/180.0))**2.0+(np.abs(elem[1]-otherelem[1])**2.0)))

        estsize=max(ndists)*3600.0
        # Use largest component size if estsize is smaller (e.g. bc of very overlapping regs in association)
        if estsize<msiz:
            estsize=msiz

        # Generate name for new source
    
        setname=getname(ramean,decmean)
            
        # put elements in associations table

        for elem in finalset:
            aflux=fluxdict.get(elem)
            amaj=sizedict.get(elem)
            amin=mindict.get(elem)
            apa=padict.get(elem)
            ara,adec=coordsdict.get(elem)
            if elem in overlap:
                oflag=1
            assoccat.append({'Assoc_Name':elem,'Source_Name':setname})
            compcat.append({'Comp_Name':elem,'Assoc':1,'Source_Name':setname,'Comp_RA':ara,'Comp_DEC':adec,'Comp_Flux':aflux,'Comp_Maj':amaj,'Comp_Min':amin,'Comp_PA':apa})

        ncomps=len(finalset)
        finalcat.append({'Source_Name':setname,'RA':ramean,'Dec':decmean,'OptID_Name':optid, 'optRA':optra, 'optDec':optdec,'Flux':nflux,'Size':estsize,'Assoc':ncomps, 'Assoc_Qual':newscore[0][1],'ID_Qual':qual,'Compoverlap':oflag,'Art_prob':apb,'Blend_prob':bpb,'Zoom_prob':zpb,'Hostbroken_prob':hpb,'Imagemissing_prob':imb,'Badclick':badclick})

            
    else:
        print("Set does not contain any LGZ subjects --- this shouldn't happen!")
            
# Now go through single sources, excluding members of sets

for subj in new_subs_unique:
    if subj in catalogued:
        overlap.append(subj)
        continue
    else:
        catalogued.append(subj)
    # Get source RA and Dec:
    oflag=0
    ncoords=coordsdict.get(subj)
    nra=ncoords[0]
    ndec=ncoords[1]
    nflux=fluxdict.get(subj)
    nmajsiz=sizedict.get(subj)
    nmin=mindict.get(subj)
    npa=padict.get(subj)
    
    catalogued.append(subj)
    apb,bpb,zpb,hpb,imb=probdict.get(subj)
        
    # Source is not in an associated source, so tabulate with its ID if it exists or none if it doesn't

    if subj in subjs_mult:
        optid="Mult"
        optra=np.nan
        optdec=np.nan
        qual=np.nan
    elif subj in subjs_noids:
        optid="None"
        optra=np.nan
        optdec=np.nan
        qual=np.nan
    else:
        findilist=[x for x in ilist if x['Source_Name']==subj]
                   #filter(lambda x: x['Source_Name'] == subj, ilist)

        if len(findilist)==0:
            optid="None"
            optra=np.nan
            optdec=np.nan
            qual=np.nan

        else:        
            optid=findilist[0]['Opt_id']
            qual=findilist[0]['Quality']
            optra,optdec=ocoordsdict.get(optid)

    if subj in overlap:
        oflag=1
    finalcat.append({'Source_Name':subj,'RA':nra,'Dec':ndec,'OptID_Name':optid, 'optRA':optra, 'optDec':optdec,'Flux':nflux,'Size':nmajsiz,'Assoc':0, 'Assoc_Qual':np.nan,'ID_Qual':qual,'Compoverlap':oflag,'Art_prob':apb,'Blend_prob':bpb,'Zoom_prob':zpb,'Hostbroken_prob':hpb,'Imagemissing_prob':imb,'Badclick':badclickdict[subj]})
    compcat.append({'Comp_Name':subj,'Assoc':0,'Source_Name':subj,'Comp_RA':nra,'Comp_DEC':ndec,'Comp_Flux':nflux,'Comp_Maj':nmajsiz,'Comp_Min':nmin,'Comp_PA':npa})
            
    

print("Num of sources not in final catalogue: ",len(nonexist))

for thing in nonexist:
    print(thing)
                    
# Write out final catalogue, associated source catalogue and artefacts catalogue

print('Writing final catalogue')

write_fits_out(['Source_Name','RA','Dec','OptID_Name','optRA','optDec','Flux','Size','Assoc','Assoc_Qual','ID_Qual','Compoverlap','Art_prob','Blend_prob','Zoom_prob','Hostbroken_prob','Imagemissing_prob','Badclick'],finalcat,"MGZ-cat.fits")
write_fits_out(['Comp_Name','Assoc','Source_Name','Comp_RA','Comp_DEC','Comp_Flux','Comp_Maj','Comp_Min','Comp_PA'],compcat,'MGZ-comps.fits')

print('Multlist contains',len(multlist),'sources')
if len(multlist)>0:
    #'Source_Name':subj,'Opt_id':nam,'Numclass':numseen,'Quality':tqual
    write_fits_out(['Source_Name','Opt_id','Numclass','Quality'],multlist,'MGZ-multiple.fits')
