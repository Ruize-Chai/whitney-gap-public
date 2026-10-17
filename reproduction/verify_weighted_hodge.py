"""Independent weighted-metric P1/Hodge energy assembly on zig-zag torus.

Checks the central scalar-to-Hodge lift numerically without using Fourier
closed form in the stiffness or load assembly. SciPy sparse direct solve.
"""
import argparse
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def scalar_formula(n, M, Delta):
    h = 1/n; eps = 1/M; theta=2*np.pi*h
    kappa=np.sin(theta/2)/h
    tau=8*np.sin(theta/4)**2/theta**2
    gamma=kappa*kappa+(1-np.cos(theta/2))**2/(4*eps*eps)
    return Delta**2/2 - 2*Delta**2*tau*tau*kappa*kappa/gamma


def integrate_system(n, M, A=2., Delta=.8, nq=12, weighted=True):
    assert n >= 3 and M >= 4 and M%2 == 0 and A>abs(Delta)>0
    h=1/n; eps=1/M; m=M//2; N=2*m*n
    nodes, weights=leggauss(nq); nodes=.5*(nodes+1); weights=.5*weights
    rows=[];cols=[]; vals=[]; load=np.zeros(N)
    target_square=0.; weighted_period_x=0.; norm_harm_y=0.
    # For all triangles directly assemble arbitrary metric quadratic energy.
    def idx(t,j,k): return (2*(j % m)+t)*n+(k % n)
    def P(j,k): return np.array([2*j/M,k/n])
    def Q(j,k): return np.array([(2*j+1)/M,(k+.5)/n])
    for j in range(m):
        for k in range(n):
            cells = [
                ([(0,j,k),(1,j,k),(0,j,k+1)],[P(j,k),Q(j,k),P(j,k+1)]),
                ([(0,j,k+1),(1,j,k),(1,j,k+1)],[P(j,k+1),Q(j,k),Q(j,k+1)]),
                ([(1,j,k),(0,j+1,k),(0,j+1,k+1)],[Q(j,k),P(j+1,k),P(j+1,k+1)]),
                ([(1,j,k),(0,j+1,k+1),(1,j,k+1)],[Q(j,k),P(j+1,k+1),Q(j,k+1)]),
            ]
            for vi,pts in cells:
                ids = [idx(*v) for v in vi]
                X=np.asarray(pts,float)
                B=np.column_stack((X[1]-X[0],X[2]-X[0]))
                bt = np.linalg.inv(B).T
                grad=np.stack((bt@[-1.,-1.],bt@[1.,0.],bt@[0.,1.]))
                det=abs(np.linalg.det(B))
                wx=wy=ws=0.
                for a,wa in zip(nodes,weights):
                    for b,wb in zip(nodes,weights):
                        y=X[0,1]+a*(X[1,1]-X[0,1])+(1-a)*b*(X[2,1]-X[0,1])
                        d=Delta*np.cos(2*np.pi*y)
                        a_metric=1./(A+d)
                        w=wa*wb*(1-a)*det
                        wx += w*(1/a_metric if weighted else 1)
                        wy += w*(a_metric if weighted else 1)
                        ws += w*(a_metric*d*d if weighted else d*d)
                        if weighted:
                            weighted_period_x += w*(1/a_metric)
                            norm_harm_y += w*a_metric*(A+d)**2
                target_square+=ws
                kel=wx*np.outer(grad[:,0],grad[:,0])+wy*np.outer(grad[:,1],grad[:,1])
                # Source term f_a=∫a ψ' ∂_yφ_a; weighted integration separately.
                fr=0.
                for a,wa in zip(nodes,weights):
                    for b,wb in zip(nodes,weights):
                        y=X[0,1]+a*(X[1,1]-X[0,1])+(1-a)*b*(X[2,1]-X[0,1])
                        d=Delta*np.cos(2*np.pi*y)
                        am=1./(A+d)
                        w=wa*wb*(1-a)*det
                        fr+=w*(am if weighted else 1)*d
                for i,id_i in enumerate(ids):
                    load[id_i] += fr*grad[i,1]
                    for k2,id_k in enumerate(ids):
                        rows.append(id_i); cols.append(id_k); vals.append(kel[i,k2])
    K=coo_matrix((vals,(rows,cols)),shape=(N,N)).tocsr()
    z=np.zeros(N)
    z[1:]=spsolve(K[1:,1:],load[1:])
    e2=target_square-load@z
    kappa_g = min(A-abs(Delta),1/(A+abs(Delta)))
    nu_g = max(A+abs(Delta),1/(A-abs(Delta)))
    assert e2 >= -1e-9
    if not weighted: return e2
    # The full vertical harmonic representative has coefficients A dy + grad chi.
    # Integrate norm directly through independent weighted quadrature. Quadratic
    # identities yield this expression using an exact integration by parts of
    # y derivative and the mesh matrices.
    # ∫ a A^2 + 2a A z_y + ∫(a^-1 z_x² + a z_y²).
    # Directly calculate over physical cells (fresh integration, not Ritz identity).
    direct_norm=0.; direct_cross=0.; direct_xcross=0.
    for j in range(m):
        for k in range(n):
            cells = [
                ([(0,j,k),(1,j,k),(0,j,k+1)],[P(j,k),Q(j,k),P(j,k+1)]),
                ([(0,j,k+1),(1,j,k),(1,j,k+1)],[P(j,k+1),Q(j,k),Q(j,k+1)]),
                ([(1,j,k),(0,j+1,k),(0,j+1,k+1)],[Q(j,k),P(j+1,k),P(j+1,k+1)]),
                ([(1,j,k),(0,j+1,k+1),(1,j,k+1)],[Q(j,k),P(j+1,k+1),Q(j,k+1)]),
            ]
            for vi,pts in cells:
                X=np.asarray(pts,float)
                B=np.column_stack((X[1]-X[0],X[2]-X[0])); bt=np.linalg.inv(B).T
                grad=np.stack((bt@[-1.,-1.],bt@[1.,0.],bt@[0.,1.]))
                coeff=z[[idx(*v) for v in vi]]@grad
                det=abs(np.linalg.det(B))
                for a,wa in zip(nodes,weights):
                    for b,wb in zip(nodes,weights):
                        y=X[0,1]+a*(X[1,1]-X[0,1])+(1-a)*b*(X[2,1]-X[0,1])
                        d=Delta*np.cos(2*np.pi*y);am=1./(A+d)
                        w=wa*wb*(1-a)*det
                        direct_norm+=w*((coeff[0]**2)/am+am*(A+coeff[1])**2)
                        direct_cross+=w*am*(A+d)*(A+coeff[1])
                        direct_xcross+=w*(coeff[0]/am)
    direct_gap2=1-direct_cross**2/(norm_harm_y*direct_norm)
    identity_gap2=e2/(A+e2)
    return dict(n=n,M=M,e2=e2,scalarE=scalar_formula(n,M,Delta),gap2=identity_gap2,direct_gap2=direct_gap2,
                cross=direct_cross,vertical_norm=direct_norm,norm_y=norm_harm_y,dxnorm=weighted_period_x,
                xcross=direct_xcross, lower=kappa_g*scalar_formula(n,M,Delta),upper=nu_g*scalar_formula(n,M,Delta))

if __name__ == '__main__':
    for n,M in [(4,8),(6,18),(8,32),(10,50)]:
        r=integrate_system(n,M)
        eu=integrate_system(n,M,weighted=False)
        print(f"(n,M)=({n},{M}) E={r['scalarE']:.12f} E_direct={eu:.12f} "
              f"e_g^2={r['e2']:.12f} gap2={r['gap2']:.12f} direct_angle2={r['direct_gap2']:.12f} "
              f"muE={r['lower']:.12f} nuE={r['upper']:.12f} ")
        assert abs(eu-r['scalarE'])<1e-9
        assert r['lower']-1e-8 <= r['e2'] <= r['upper']+1e-8
        assert abs(r['direct_gap2']-r['gap2'])<1e-9
        assert abs(r['dxnorm']-2)<1e-9 and abs(r['norm_y']-2)<1e-9
        assert abs(r['cross']-2)<1e-9 and abs(r['xcross'])<1e-9
    print('PASS: weighted Ritz error, two-sided bound, and independently integrated principal-angle formula.')
