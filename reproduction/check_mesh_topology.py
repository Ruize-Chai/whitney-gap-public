"""Combinatorial check of the periodic zig-zag triangulation.

Verifies unique simplices for M>=4, pairwise edge incidence, and
Euler characteristic of the torus on representative meshes.
"""
from collections import Counter

def mesh_counts(n,M):
    assert M%2==0 and n>=3
    m=M//2
    def P(j,k):return (0,j%m,k%n)
    def Q(j,k):return (1,j%m,k%n)
    faces=[]
    for j in range(m):
        for k in range(n):
            faces.extend([
                (P(j,k),Q(j,k),P(j,k+1)),
                (P(j,k+1),Q(j,k),Q(j,k+1)),
                (Q(j,k),P(j+1,k),P(j+1,k+1)),
                (Q(j,k),P(j+1,k+1),Q(j,k+1)),
            ])
    vertices={v for f in faces for v in f}
    edges=Counter(tuple(sorted((f[i],f[(i+1)%3]))) for f in faces for i in range(3))
    unique_faces=len({frozenset(f) for f in faces})
    no_degenerate=all(len(set(f))==3 for f in faces)
    return dict(n=n,M=M,V=len(vertices),E=len(edges),F=len(faces),Euler=len(vertices)-len(edges)+len(faces),
                simple_faces=unique_faces==len(faces),edge_incidence_two=all(c==2 for c in edges.values()),nondegenerate=no_degenerate)

if __name__=='__main__':
    for n,M in [(3,2),(3,4),(4,4),(5,6),(8,32)]:
        x=mesh_counts(n,M);print(x)
        if M>=4:assert x['Euler']==0 and x['simple_faces'] and x['edge_incidence_two'] and x['nondegenerate']
