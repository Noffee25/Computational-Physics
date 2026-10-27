import numpy as np
from scipy import stats
names={1:'RANDU',2:'RANDM',3:'RAND',4:'MINSTD',5:'R250',6:'RANNUM',7:'MT19937',8:'MIDSQ',9:'LOGISTIC'}
def load(fn):
    hdr={}
    for l in open(fn):
        if l.startswith('# N '): hdr['N']=int(l.split()[2])
        if l.startswith('# sumx'): hdr['sx'],hdr['sxx']=map(float,l.split()[3:5])
    d=np.loadtxt(fn,comments='#')
    return hdr,d
for g in range(1,10):
    h,d=load(f'out/lag/lag_{g}.txt'); N=h['N']; k=d[:,0]
    cnt=N-k
    r=d[:,1:]/cnt[:,None]; z=r*np.sqrt(cnt)[:,None]
    Q=(z[:,0]**2).sum(); pQ=stats.chi2.logsf(Q,300)/np.log(10)
    Q2=[(z[:64,c]**2).sum() for c in (1,2,3)]
    imax=np.argmax(np.abs(z[:,0]))
    print('%-9s N=%d r11(1)=%.3e z=%.2f | max|z11| %.2f at k=%d | Q300=%.0f log10p=%.1f | Q64(12,21,22)=%s | z21(1)=%.1f z22(1)=%.1f'%(
        names[g],N,r[0,0],z[0,0],abs(z[imax,0]),k[imax],Q,pQ,['%.0f'%q for q in Q2],z[0,2],z[0,3]))
