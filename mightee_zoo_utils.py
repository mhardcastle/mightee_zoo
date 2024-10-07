#!/usr/bin/env python
# MIGHTEE zoo upload (and eventually other) tools

from __future__ import print_function
from astropy.table import Table
from panoptes_client import Panoptes, Project, SubjectSet, Subject, Workflow
import os
import glob
from time import sleep
import datetime
import threading

target='/beegfs/car/mjh/mightee_zoo/'
for k in ['PANOPTES_PASSWORD']:
    if k not in os.environ:
        raise RuntimeError('Variable %s not set' % k)

def upload_images(name,maxcount=None):
    print('Create subject set and upload images for',name)
    wd=os.getcwd()
    Panoptes.connect(username='mjh22',password=os.environ['PANOPTES_PASSWORD'])
    os.chdir(target+name)
    project = Project.find(slug='mjh22/mightee-galaxy-zoo')

    # check for subject sets already existing
    ssets=project.links.subject_sets
    for s in ssets:
        if s.display_name==name:
            # append to this subject set
            print('Found subject set',s,'already exists')
            subject_set=s
            uploaded=[]
            for subject in subject_set.subjects:
                uploaded.append(subject.metadata['source_name'])
            new_set=False
            break
    else:
        new_set=True
        uploaded=[]
        subject_set = SubjectSet()
        subject_set.display_name=name
        subject_set.links.project=project
        subject_set.save()
        print('Made subject set')

    new_subjects = []
    count=0
    g=glob.glob('*-manifest.txt')
    for i,f in enumerate(g):
        bits=open(f).readlines()[0].split(',')
        sourcename=bits[5]
        if sourcename in uploaded:
            print('Skipping',sourcename,'as already uploaded')
            continue
        metadata={'subject_id':int(bits[0]),'ra':float(bits[6]),'dec':float(bits[7]),'#size':float(bits[8]),'source_name':sourcename}
        print('Upload doing',bits[5],'%i/%i' % (i,len(g)))
        subject = Subject()
        subject.links.project = project
        subject.metadata.update(metadata)
        for location in bits[1:5]:
            subject.add_location(location)
        subject.save()
        new_subjects.append(subject)
        count+=1
        if maxcount and count>maxcount:
            print('Max batch number exceeded, stopping')
            break

    subject_set.add(new_subjects)

    if new_set:
        print('Linking new subject set')
        workflow=Workflow(24522)
        workflow.links.subject_sets.add(subject_set)
    print('Done!')
    
if __name__=='__main__':    
    import sys
    try:
        maxc=int(sys.argv[2])
    except:
        maxc=None
    upload_images(sys.argv[1],maxcount=maxc)
    
