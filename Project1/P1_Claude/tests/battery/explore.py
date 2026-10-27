import numpy as np
from fractions import Fraction
# MIDSQ trajectories
def midsq(s):
    seen={}
    n=0
    while s not in seen:
        seen[s]=n; s=(s*s//10000)%100000000; n+=1
    return n, n-seen[s], s
for s in [777, 12345678, 31415926, 27182818, 77777777, 10101011]:
    print('MIDSQ seed',s, 'steps to cycle, cycle len, entry', midsq(s))
best=max(range(1,200000,7), key=lambda s: midsq(s)[0])
print('best small seed', best, midsq(best))
# modular powers
m=2**31
for a in (65539,6539):
    print(a, 'a^(2^28) mod 2^31 =', pow(a,2**28,m), '2^30+1=',2**30+1, ' a^(2^27)=',pow(a,2**27,m))
# LLL spectral test
def lll(B, delta=Fraction(99,100)):
    import copy
    B=[list(map(int,b)) for b in B]
    n=len(B)
    def dot(u,v): return sum(x*y for x,y in zip(u,v))
    def gs(B):
        Bs=[]; mu=[[Fraction(0)]*n for _ in range(n)]
        for i in range(n):
            v=[Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j]=Fraction(dot(B[i],Bs[j]))/dot(Bs[j],Bs[j]) if dot(Bs[j],Bs[j]) else 0
                v=[a-mu[i][j]*b for a,b in zip(v,Bs[j])]
            Bs.append(v)
        return Bs,mu
    k=1
    Bs,mu=gs(B)
    while k<n:
        for j in range(k-1,-1,-1):
            q=round(mu[k][j])
            if q:
                B[k]=[a-q*b for a,b in zip(B[k],B[j])]
                Bs,mu=gs(B)
        if dot(Bs[k],Bs[k]) >= (delta-mu[k][k-1]**2)*dot(Bs[k-1],Bs[k-1]):
            k+=1
        else:
            B[k],B[k-1]=B[k-1],B[k]; Bs,mu=gs(B); k=max(k-1,1)
    return B
def dual_short(a,m,d):
    # dual lattice: h with h0 + h1 a + ... = 0 mod m ; basis: (m,0..), (-a,1,0..), (-a^2,0,1..)...
    B=[]
    B.append([m]+[0]*(d-1))
    for i in range(1,d):
        v=[0]*d; v[0]=-pow(a,i,m); v[i]=1; B.append(v)
    R=lll(B)
    R.sort(key=lambda v: sum(x*x for x in v))
    return R[0], sum(x*x for x in R[0])**0.5
for name,a,mm in [('RANDU',65539,2**31),('RANDM',6539,2**31),('MINSTD',16807,2**31-1)]:
    for d in range(2,7):
        v,L=dual_short(a,mm,d)
        print(name,d,'shortest dual vector',v,'length %.1f'%L,'plane spacing %.4g'%(1/L))
