"""Independent direct P1 triangle assembly against the closed Fourier formula.

No original project implementation is imported. Periodic degrees of freedom
are identified modulo M and n, while triangles use unwrapped coordinates.
"""
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.polynomial.legendre import leggauss

def prediction(n, M, Delta=0.8):
    h = 1.0/n; eps=1.0/M
    theta=2*np.pi*h
    kappa=np.sin(theta/2)/h
    tau=8*np.sin(theta/4)**2/theta**2
    gam=kappa*kappa+(1-np.cos(theta/2))**2/(4*eps**2)
    return Delta**2/2 - 2*Delta**2*tau**2*kappa**2/gam

def direct(n, M, Delta=0.8, qorder=12):
    assert M%2==0 and n>=3
    m=M//2; N=2*m*n
    K=lil_matrix((N,N)); b=np.zeros(N)
    nodes,w=leggauss(qorder); nodes=(nodes+1)/2; w=w/2
    def idx(t,j,k):return (2*(j % m)+t)*n+(k % n)
    def P(j,k):return np.array([2*j/M,k/n],dtype=float)
    def Q(j,k):return np.array([(2*j+1)/M,(k+0.5)/n],dtype=float)
    for j in range(m):
        for k in range(n):
            pieces=[
                ([(0,j,k),(1,j,k),(0,j,k+1)],[P(j,k),Q(j,k),P(j,k+1)]),
                ([(0,j,k+1),(1,j,k),(1,j,k+1)],[P(j,k+1),Q(j,k),Q(j,k+1)]),
                ([(1,j,k),(0,j+1,k),(0,j+1,k+1)],[Q(j,k),P(j+1,k),P(j+1,k+1)]),
                ([(1,j,k),(0,j+1,k+1),(1,j,k+1)],[Q(j,k),P(j+1,k+1),Q(j,k+1)]),
            ]
            for ids,pts in pieces:
                loc=[idx(*i) for i in ids]; x0,x1,x2=pts
                B=np.stack([x1-x0,x2-x0],axis=1)
                area2=abs(np.linalg.det(B))
                grads=np.stack([-np.linalg.inv(B).T.sum(axis=1),*np.linalg.inv(B).T.T],axis=0)
                # safer: gradients = B^{-T} [-1,-1], B^{-T}[1,0], B^{-T}[0,1]
                bt=np.linalg.inv(B).T
                grads=np.stack([bt@np.array([-1,-1]),bt@np.array([1,0]),bt@np.array([0,1])])
                for a in range(3):
                    for c in range(3): K[loc[a],loc[c]]+=area2/2*(grads[a]@grads[c])
                source=0.
                for s,ws in zip(nodes,w):
                    for t,wt in zip(nodes,w):
                        u=s;v=(1-s)*t
                        yy=x0[1]+u*(x1[1]-x0[1])+v*(x2[1]-x0[1])
                        source+=ws*wt*(1-s)*np.cos(2*np.pi*yy)
                source*=area2*Delta
                for a in range(3): b[loc[a]]+=source*grads[a,1]
    # on a torus, K has the constant null direction: fix one degree of freedom
    KK=K.tocsr()[1:,1:]; z=spsolve(KK,b[1:]); return Delta**2/2-np.dot(b[1:],z)

if __name__=='__main__':
    for n,M in [(4,8),(6,18),(8,32),(10,50),(12,72)]:
        closed=prediction(n,M); measured=direct(n,M)
        err=abs(closed-measured)
        print(f'n={n:2d}, M={M:3d}, closed={closed:.12f}, direct={measured:.12f}, abs_err={err:.2e}')
        assert err<1e-9, f'Direct assembly does not agree with closed formula: {err}'
    print('PASS: independent direct P1 assembly matches the exact Fourier formula.')
