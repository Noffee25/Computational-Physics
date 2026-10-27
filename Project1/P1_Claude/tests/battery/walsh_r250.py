# Exact E[f(X) f(Y) f(X xor Y)] for independent uniform n-bit words via the
# Walsh-Hadamard transform: E = sum_S fhat(S)^3, fhat = WHT(f)/2^n.
import numpy as np
def fwht(a):
    a=a.copy(); h=1; n=len(a)
    while h<n:
        a=a.reshape(-1,2*h); x=a[:,:h].copy(); y=a[:,h:].copy()
        a[:,:h]=x+y; a[:,h:]=x-y; a=a.reshape(-1); h*=2
    return a
nb=22
w=np.arange(2**nb,dtype=np.float64)
x=(w+0.5)/2**nb      # cell midpoints -> converges to the continuum value as nb grows
for name,f in [('P1',np.sqrt(3)*(2*x-1)),('P2',np.sqrt(5)*(6*x*x-6*x+1))]:
    fh=fwht(f)/2**nb
    print(name,'E[f f f(xor)] =',(fh**3).sum(), ' E[f^2]=',(f*f).mean())
print('theory P1: -12^1.5/56 =', -12**1.5/56)
# variance of the P1 triple product under the XOR alternative (Monte Carlo)
rng=np.random.default_rng(1)
a=rng.integers(0,2**32,size=10**7,dtype=np.uint64); b=rng.integers(0,2**32,size=10**7,dtype=np.uint64); c=a^b
X,Y,Z=[(v.astype(np.float64)/2**32) for v in (a,b,c)]
t=np.sqrt(3)**3*(2*X-1)*(2*Y-1)*(2*Z-1)
print('MC mean %.5f sd %.4f'%(t.mean(),t.std()), ' SE at N=2^29: %.2e'%(t.std()/2**14.5))
t2=np.sqrt(5)**3*(6*X*X-6*X+1)*(6*Y*Y-6*Y+1)*(6*Z*Z-6*Z+1)
print('MC P2 mean %.5f sd %.4f'%(t2.mean(),t2.std()))
