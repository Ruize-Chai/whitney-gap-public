from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(parents=True, exist_ok=True)

def E(n,M,delta=0.8):
    h=1/n; eps=1/M; theta=2*np.pi*h
    kap=np.sin(theta/2)/h
    tau=8*np.sin(theta/4)**2/theta**2
    gam=kap**2+(1-np.cos(theta/2))**2/(4*eps**2)
    return delta**2/2-2*delta**2*tau**2*kap**2/gam

ns=np.array([8,12,16,24,32,48,64,96,128,192,256])
fig,ax=plt.subplots(figsize=(6.2,3.5))
curves=[(1.5,'o','-',r'convergent: $\varepsilon\asymp h^{3/2}$'),
        (2.,'s','--',r'critical: $\varepsilon\asymp h^2$'),
        (2.5,'^','-.',r'locking: $\varepsilon\asymp h^{5/2}$')]
for p,mark,linestyle,label in curves:
    vals=[E(n,max(2,2*int(np.rint(n**p/2)))) for n in ns]
    ax.plot(ns,vals,color='0.10',linestyle=linestyle,marker=mark,markersize=3.2,
            markerfacecolor='white',markeredgewidth=0.75,linewidth=1.15,label=label)
ax.axhline(0.8**2/2,linestyle=':',linewidth=.95,color='0.50',label=r'locking limit $\Delta^2/2$')
ax.axhline(0.8**2*np.pi**2/(2*(16+np.pi**2)),linestyle=(0,(4,3)),linewidth=.90,color='0.53',label=r'critical limit $\Lambda_1$')
ax.set_xscale('log',base=2)
ax.set_xlabel(r'$n=1/h$');ax.set_ylabel(r'scalar minimum $\mathsf{E}_{h,\varepsilon}$')
ax.set_ylim(-0.012,.36)
ax.grid(alpha=.15,linewidth=.55)
ax.legend(loc='upper right',frameon=False,fontsize=7.7)
fig.tight_layout(pad=.6)
fig.savefig(OUT/'fig_scalar_regimes.pdf')
fig.savefig(OUT/'fig_scalar_regimes.png',dpi=240)
plt.close(fig)

M=8;n=5
fig,ax=plt.subplots(figsize=(6.3,2.5))
for j in range(M//2):
 for k in range(n):
  P=lambda jj,kk:(2*jj/M,kk/n)
  Q=lambda jj,kk:((2*jj+1)/M,(kk+.5)/n)
  tris=[(P(j,k),Q(j,k),P(j,k+1)),(P(j,k+1),Q(j,k),Q(j,k+1)),
        (Q(j,k),P(j+1,k),P(j+1,k+1)),(Q(j,k),P(j+1,k+1),Q(j,k+1))]
  for tri in tris:
   a=np.array([*tri,tri[0]])
   ax.plot(a[:,0],a[:,1],color='0.23',linewidth=.58)
ax.set_xlim(0,1);ax.set_ylim(0,1);ax.set_aspect('equal',adjustable='box')
ax.set_xlabel('$x$');ax.set_ylabel('$y$')
ax.set_xticks([0,.25,.5,.75,1]);ax.set_yticks([0,.5,1]);ax.tick_params(labelsize=8)
fig.tight_layout(pad=.7)
fig.savefig(OUT/'fig_zigzag_mesh.pdf')
plt.close(fig)
