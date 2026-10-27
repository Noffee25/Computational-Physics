import numpy as np
from scipy import stats
names={1:'RANDU',2:'RANDM',3:'RAND',4:'MINSTD',5:'R250',6:'RANNUM',7:'MT19937',8:'MIDSQ',9:'LOGISTIC'}
P={1:2**29,2:2**29,3:2**31-2,4:2**31-2}
K=2**20; L=2**23; KC=1024
for g in range(1,10):
    d=np.loadtxt(f'out/sl/sl_{g}.txt',comments='#')
    X=d[:,1]; Xc=d[:,2]
    pf=stats.chi2.sf(X,K-1); pc=stats.chi2.sf(Xc,KC-1)
    zf=(X-(K-1))/np.sqrt(2*(K-1)); zc=(Xc-(KC-1))/np.sqrt(2*(KC-1))
    ksf=stats.kstest(pf,'uniform'); ksc=stats.kstest(pc,'uniform')
    pred=-(K-1)*(L/P[g])/np.sqrt(2*(K-1)) if g in P else 0
    cum=d[:,3]/(K-1)
    print('%-8s fine: mean z %.2f (pred %.2f) sd z %.2f  KS D=%.3f p=%.2e | coarse: mean z %.2f KS D=%.3f p=%.2e | cum X2/df at 1/4,1/2,1: %.3f %.3f %.4f'%(
        names[g],zf.mean(),pred,zf.std(ddof=1),ksf.statistic,ksf.pvalue,zc.mean(),ksc.statistic,ksc.pvalue,cum[15],cum[31],cum[63]))
