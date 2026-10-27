import math, numpy as np, sys
sys.path.insert(0,'.')
from spec_verify import shortest
from coset_bruteforce import h_e
for (a,e,t) in [(27,11,3),(2443,12,2)]:
    M=2**e
    for k in range(1,e):
        n2,sv,_=shortest(a,2**k,t)
        s=sv[0]
        g=sum(c*pow(a,j,M) for j,c in enumerate(s))%M
        v2 = (g & -g).bit_length()-1 if g else e
        # empirical
        xs=[1]
        while True:
            y=a*xs[-1]%M
            if y==1: break
            xs.append(y)
        r=np.unique([(x*g)%M for x in xs])
        gaps=np.diff(np.concatenate([r,[r[0]+M]]))
        print(a,e,t,'k',k,'s',s,'|s|^2',n2,'v2(g)',v2,'h(k)',h_e(k,e),'h(v2)',h_e(v2,e),'emp gap',gaps.max()/M, 'formula',h_e(k,e)/math.sqrt(n2))
