from fractions import Fraction as F
import math
def dedekind(h,k):
    # s(h,k) via reciprocity, gcd(h,k)=1
    s=F(0); sign=1
    h%=k
    while h!=0:
        s+= sign*(F(h,k)+F(k,h)+F(1,h*k))/12 - sign*F(1,4)
        # s(h,k) = -s(k,h) + (...) ; s(k,h)=s(k mod h, h)
        h,k = k%h, h
        sign=-sign
    return s
# check small brute force
def brute(h,k):
    def saw(x):
        x=x-math.floor(x)
        return 0 if x==0 else x-0.5
    return sum(saw(i/k)*saw(h*i/k) for i in range(1,k))
for h,k in [(3,7),(5,13),(11,64+1)]:
    print(h,k,float(dedekind(h,k)),brute(h,k))
m=2**31
for name,a in [('RANDU',65539),('RANDM',6539)]:
    out=[]
    for k in range(1,301):
        b=pow(a,k,m)
        s=dedekind(b,m)
        r=12*m*s/((m-1)*(m-2))
        out.append((k,b,float(r)))
    out_sorted=sorted(out,key=lambda t:-abs(t[2]))
    print(name,'lag1..5', [ '%.4e'%t[2] for t in out[:5]], '1/a=%.4e'%(1/a))
    print(name,'largest |r_k|:', [(t[0],'%.3e'%t[2]) for t in out_sorted[:8]])
    print(name,'sqrt(2^29)*r:', [(t[0],'%.2f'%(t[2]*2**14.5)) for t in out_sorted[:5]])
