import numpy as np
import cilmesh
import cuadratura
import fs_gf

X1, Y1 = 3., 3.
X2, Y2 = 6., 7. 

R = ((X2-X1)**2 + (Y2-Y1)**2)**0.5
tx = (X1-X2)/R
ty = (Y1-Y2)/R

nx = -ty
ny =  tx

X0, Y0 = 0., 0.

n_gsch = 5

t = nodos_gauss_chebyshev(n_gsch)
Xg, Yg = nodos_gauss_panel(t, X1, Y1, X2, Y2)

print t
print Xg
print Yg


Txxx, Txxy, Tyxx, Tyxy, Txyx, Txyy, Tyyx, Tyyy = gf_stress_fs(Xg,Yg,X0,Y0)

bxx, bxy, byx, byy = beta_comp_gs(t,Xg,Yg, X0, Y0, nx, ny, R)

print
print bxx
print bxy
print byx
print byy

