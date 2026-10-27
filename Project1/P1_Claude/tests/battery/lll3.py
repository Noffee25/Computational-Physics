import sys; sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/battery")
from fractions import Fraction
import itertools, numpy as np
def lll(B, delta=Fraction(99,100)):
    B=[list(map(int,b)) for b in B]; n=len(B)
    def dot(u,v): return sum(x*y for x,y in zip(u,v))
    def gs(B):
        Bs=[]; mu=[[Fraction(0)]*n for _ in range(n)]
        for i in range(n):
            v=[Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j]=Fraction(dot(B[i],Bs[j]))/dot(Bs[j],Bs[j])
                v=[a-mu[i][j]*b for a,b in zip(v,Bs[j])]
            Bs.append(v)
        return Bs,mu
    k=1; Bs,mu=gs(B)
    while k<n:
        for j in range(k-1,-1,-1):
            q=round(mu[k][j])
            if q: B[k]=[a-q*b for a,b in zip(B[k],B[j])]; Bs,mu=gs(B)
        if dot(Bs[k],Bs[k])>=(delta-mu[k][k-1]**2)*dot(Bs[k-1],Bs[k-1]): k+=1
        else: B[k],B[k-1]=B[k-1],B[k]; Bs,mu=gs(B); k=max(k-1,1)
    return B
def short(a,j,k,M):
    b1=pow(a,k-j,2**31); b2=pow(a,k,2**31)   # x_{n-j}=b1 x_{n-k}, x_n=b2 x_{n-k}
    B=[[M,0,0],[-b1,1,0],[-b2,0,1]]
    R=lll(B)
    vs=set()
    for c in itertools.product(range(-2,3),repeat=3):
        v=tuple(sum(c[i]*R[i][t] for i in range(3)) for t in range(3))
        if any(v): vs.add(v)
    return sorted(vs,key=lambda v:sum(x*x for x in v))[:6]
for a,j,k in [(65539,69,213),(6539,5,188),(65539,1,2)]:
    for M in (2**31,2**30,2**29,2**28):
        print(a,(j,k),'mod 2^%d'%int(np.log2(M)),short(a,j,k,M)[:3])
