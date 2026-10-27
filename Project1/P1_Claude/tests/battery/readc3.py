import numpy as np, sys
def readc3(fn):
    with open(fn,'rb') as f:
        K=np.fromfile(f,dtype=np.int32,count=1)[0]
        N=np.fromfile(f,dtype=np.int64,count=1)[0]
        C=np.fromfile(f,dtype=np.float64).reshape((K,2*K)).T   # Fortran order (2K,K)
    j=np.arange(1,K+1)
    cnt=N-np.maximum.outer(j,j)
    T1=C[:K,:]/cnt; T2=C[K:,:]/cnt
    Z1=T1*np.sqrt(cnt); Z2=T2*np.sqrt(cnt)
    return K,N,T1,T2,Z1,Z2
if __name__=='__main__':
    K,N,T1,T2,Z1,Z2=readc3(sys.argv[1])
    iu=np.triu_indices(K,1)
    for name,Z,T in [('T1',Z1,T1),('T2',Z2,T2)]:
        z=Z[iu]
        o=np.argsort(-np.abs(z))[:5]
        print(name,'N',N,'mean z %.3f sd z %.3f max|z| %.2f'%(z.mean(),z.std(),np.abs(z).max()))
        for t in o: print('   (j,k)=(%d,%d) z=%.1f T=%.5f'%(iu[0][t]+1,iu[1][t]+1,z[t],T[iu][t]))
