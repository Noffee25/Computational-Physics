import numpy as np
from scipy import stats
names={1:'RANDU',2:'RANDM',4:'MINSTD',7:'MT19937',9:'LOGISTIC'}
P={1:2**29,2:2**29,4:2**31-2,7:np.inf,9:np.inf}
for g in (1,2,4,7,9):
    for lg in range(20,31):
        fn=f'out/lagN/lag_{g}_{lg}.txt' if lg!=29 else f'out/lag/lag_{g}.txt'
        try: d=np.loadtxt(fn,comments='#')
        except Exception: continue
        N=2**lg; k=d[:,0]; z=d[:,1]/np.sqrt(N-k)
        Q=(z**2).sum(); Q2=(d[:64,4]**2/(N-k[:64])).sum()
        print('%-8s N=2^%d N/P=%.4f Q300=%.1f (Q/300=%.3f, 1-N/P=%.3f) z(Q)=%.1f Q64_22=%.1f'%(names[g],lg,N/P[g],Q,Q/300,1-N/P[g],(Q-300)/np.sqrt(600),Q2))
