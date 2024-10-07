# take the output CSV file and split it per field

import csv
import json

infile=open('mightee-galaxy-zoo-classifications.csv')
cosmos=open('cosmos-classifications.csv','w')
cdfs=open('cdfs-classifications.csv','w')
xmm=open('xmmlss-classfications.csv','w')

line=infile.readline()
cosmos.write(line)
cdfs.write(line)
xmm.write(line)
reader=csv.reader(infile)
cosmos_writer=csv.writer(cosmos)
cdfs_writer=csv.writer(cdfs)
xmm_writer=csv.writer(xmm)
for row in reader:
    print(row[0],row[1])
    sd=json.loads(row[12])
    subj=list(sd.keys())[0]
    kd=sd[subj]
    ra=kd['ra']
    if 52<ra<54:
        writer=cdfs_writer
    elif 33<ra<38:
        writer=xmm_writer
    elif 149<ra<151:
        writer=cosmos_writer
    else:
        writer=None
    writer.writerow(row)
    
    
    
