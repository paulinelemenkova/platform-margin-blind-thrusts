#!/usr/bin/env python3
import numpy as np, math
from scipy.spatial import cKDTree
LAT0=36.65; KX=111.32*math.cos(math.radians(LAT0)); KY=111.195
EQ=(67.3886,36.5461)
def segs(fn):
    S=[];cur=[];z=None
    for l in open(fn):
        if l.startswith('>'):
            if cur:S.append((z,np.array(cur)))
            cur=[]; z=float(l.split('-Z')[1].split()[0]) if '-Z' in l else None
        elif l.strip() and l[0]!='#':
            t=l.split();cur.append((float(t[0]),float(t[1])))
    if cur:S.append((z,np.array(cur)))
    return S
xy=lambda a:np.column_stack([(a[:,0]-66)*KX,(a[:,1]-LAT0)*KY])
def dens(S):
    pts=[]
    for _,s in S:
        p=xy(s)
        for a,b in zip(p[:-1],p[1:]):
            n=max(1,int(np.hypot(*(b-a))/0.2)); pts.append(a+(b-a)*np.linspace(0,1,n+1)[:,None])
    return np.vstack(pts)
T=cKDTree(dens(segs('qf.gmt')))
towns={"C01":(67.110,36.709),"C02":(67.686,36.698),"C03":(66.897,36.758),"C04":(68.018,36.262),"C05":(65.752,36.667),"C06":(68.860,36.728)}
rows=[]
for z,s in segs('struct_cand.gmt'):
    p=xy(s); d=np.hypot(*(p[:,None,:]-p[None,:,:]).transpose(2,0,1)); i,j=np.unravel_index(d.argmax(),d.shape)
    L=d[i,j]; dx,dy=p[j]-p[i]; az=math.degrees(math.atan2(dx,dy))%180
    c=s.mean(0); deq=math.hypot((c[0]-EQ[0])*KX,(c[1]-EQ[1])*KY)
    if L<15 or deq>75: continue
    dq=T.query(dens([(z,s)]))[0].min()
    nt,dt=min(((n,math.hypot((c[0]-a)*KX,(c[1]-b)*KY)) for n,(a,b) in towns.items()),key=lambda t:t[1])
    rows.append((c[0],c[1],L,az,z,dq,deq,nt,dt,4.49+1.49*math.log10(L)))
rows.sort()
for k,r in enumerate(rows,1): print("S%d %.3f %.3f L%.1f st%.0f p%.2f dQF%.1f dEQ%.0f %s %.0f Mw%.1f"%((k,)+r))
