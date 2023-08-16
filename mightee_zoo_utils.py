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

def upload_images(id):
    print('Create subject set and upload images for',id)
    wd=os.getcwd()
    Panoptes.connect(username='mjh22',password=os.environ['PANOPTES_PASSWORD'])
    os.chdir(target+id)
    project = Project.find(slug='mjh22/mightee-galaxy-zoo')
    subject_set = SubjectSet()

    subject_set.display_name=id
    subject_set.links.project=project
    subject_set.save()
    print('Made subject set')
    new_subjects = []
    g=glob.glob('*-manifest.txt')
    for i,f in enumerate(g):
        bits=open(f).readlines()[0].split(',')
        metadata={'subject_id':int(bits[0]),'ra':float(bits[6]),'dec':float(bits[7]),'#size':float(bits[8]),'source_name':bits[5]}
        print('Upload doing',bits[5],'%i/%i' % (i,len(g)))
        subject = Subject()
        subject.links.project = project
        subject.metadata.update(metadata)
        for location in bits[1:5]:
            subject.add_location(location)
        subject.save()
        new_subjects.append(subject)

    subject_set.add(new_subjects)

    workflow=Workflow(24522)
    workflow.links.subject_sets.add(subject_set)
    print('Done!')
    
if __name__=='__main__':    
    import sys
    upload_images(sys.argv[1])
    
