import numpy as np, json, math
from scipy import stats
from fractions import Fraction
names={1:'RANDU',2:'RANDM',3:'RAND',4:'MINSTD',5:'R250',6:'RANNUM',7:'MT19937',8:'MIDSQ',9:'LOGISTIC'}
TESTS=['RUNS','GAP','PERM','COUP','MAXT','BDAY','LAG1']
def stirling2(n,k):
    S=[[0]*(k+1) for _ in range(n+1)]; S[0][0]=1
    for i in range(1,n+1):
        for j in range(1,min(i,k)+1):
            S[i][j]=j*S[i-1][j]+S[i-1][j-1]
    return S[n][k]
def probs():
    P={}
    r=np.arange(1,6); P['RUNS']=np.append(r/np.array([math.factorial(x+1) for x in r]),1/720)
    p=0.1; P['GAP']=np.append(p*(1-p)**np.arange(40),(1-p)**40)
    P['PERM']=np.full(24,1/24)
    d=8; t=40
    pc=[Fraction(math.factorial(d),d**rr)*stirling2(rr-1,d-1) for rr in range(d,t)]
    pc.append(1-Fraction(math.factorial(d),d**(t-1))*stirling2(t-1,d))
    P['COUP']=np.array([float(x) for x in pc])
    P['MAXT']=np.full(100,0.01)
    P['BDAY']=np.array(json.load(open('out/bday_null.json'))['probs'])
    return P
PR=probs()
for k,v in PR.items(): assert abs(v.sum()-1)<1e-9,(k,v.sum())
def parse(fn):
    L=open(fn).read().split('\n'); d={}
    for i,l in enumerate(L):
        for t in TESTS:
            if l.startswith(t):
                d[t]=np.array(L[i+1].split(),float)
                d[t+'_hdr']=l
    return d
def chi2p(obs,p):
    n=obs.sum()
    if n==0: return np.nan, -np.inf
    e=n*p; X=((obs-e)**2/e).sum(); df=len(p)-1
    return X, stats.chi2.logsf(X,df)/np.log(10)
def analyse(fn):
    d=parse(fn); res={}
    for t in TESTS[:-1]:
        X,lp=chi2p(d[t],PR[t]); res[t]=(X,lp,d[t].sum())
    s,cnt=d['LAG1']; z=s/np.sqrt(cnt)
    res['LAG1']=(z, (stats.norm.logsf(abs(z))+np.log(2))/np.log(10), cnt)
    res['hdr']={t:d.get(t+'_hdr','') for t in TESTS}
    return res
if __name__=='__main__':
    import sys
    files=sys.argv[1:] or [f'out/kn/kn_{g}_v{v}.txt' for v in (0,1) for g in range(1,10)]
    for fn in files:
        r=analyse(fn)
        print('%-26s'%fn.split('/')[-1], '  '.join('%s %6.2f'%(t, r[t][1]) for t in TESTS), r['hdr']['RUNS'].strip(), r['hdr']['COUP'].strip())
