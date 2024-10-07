from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys

t=Table.read(sys.argv[1])
ilt=[]
if t['RA'].unit==u.deg:
    sc=SkyCoord(t['RA'],t['DEC'],frame='icrs')
else:
    sc=SkyCoord(t['RA']*u.deg,t['DEC']*u.deg,frame='icrs')
strings=sc.to_string(style='hmsdms',sep='',precision=2)
for s in strings:
    ilt.append(str('MGTC'+s).replace(' ','')[:-1])
t['Source_Name']=ilt
if 'Maj' not in t.columns:
    # insert dummy
    t['Maj']=60
    t['Min']=60
    t['PA']=0

t.write(sys.argv[2])
