import D2_standard as ops
import numpy as np
import scipy.sparse.linalg as spsplg
import scipy.sparse as spsp
from typing import NamedTuple

# Class storing 1D SBP operators
class SBP_1D(NamedTuple):
    H: spsp._arrays.csc_array
    HI: spsp._arrays.csc_array
    D1: spsp._arrays.csc_array
    xvec: np.ndarray
    e_l: spsp._arrays.csc_array
    e_r: spsp._arrays.csc_array
    h: float
    m: int
    
# Class storing 2D SBP operators
class SBP_2D(NamedTuple):
    H: spsp._arrays.csc_array
    HI: spsp._arrays.csc_array
    Dx: spsp._arrays.csc_array
    Dy: spsp._arrays.csc_array
    H_w: spsp._arrays.csc_array
    H_e: spsp._arrays.csc_array
    H_s: spsp._arrays.csc_array
    H_n: spsp._arrays.csc_array
    e_w: spsp._arrays.csc_array
    e_e: spsp._arrays.csc_array
    e_s: spsp._arrays.csc_array
    e_n: spsp._arrays.csc_array
    n_w: list
    n_e: list
    n_s: list
    n_n: list
    t_w: list
    t_e: list
    t_s: list
    t_n: list
    N: int
    m: list

# Construct 1D central SBP operators.
# Input:
#   m - number of grid points (integer)
#   lims - domain limits (tuple)
#   order - order of accuracy of operators
# 
# Output:
# SBP_1D - object containing the operators
def get_1D_sbp(m,lims,order):
    xl,xr = lims
    xvec,h = np.linspace(xl,xr,m,retstep=True)
    if order == 2:
        H,HI,D1,D2,e_l,e_r,d1_l,d1_r = ops.sbp_cent_2nd(m,h)
    elif order == 4:
        H,HI,D1,D2,e_l,e_r,d1_l,d1_r = ops.sbp_cent_4th(m,h)
    elif order == 6:
        H,HI,D1,D2,e_l,e_r,d1_l,d1_r = ops.sbp_cent_6th(m,h)
            
    H = spsp.csc_array(H)
    HI = spsp.csc_array(HI)
    D1 = spsp.csc_array(D1)
    e_l = spsp.csc_array(e_l)
    e_r = spsp.csc_array(e_r)
    return SBP_1D(H,HI,D1,xvec,e_l,e_r,h,m)

# Merge two SBP_1D objects into one, with interface conditions 
# imposed using the embedding method.
# Input:
#   SBP_1 - object containing the operators of the left block
#   SBP_2 - object containing the operators of the right block
# 
# Output:
# SBP_1D - object containing the merged operators
def merge_blocks(SBP_1,SBP_2):
    m1 = SBP_1.m
    m2 = SBP_2.m
    
    N = m1 + m2 - 1
    
    E = spsp.eye(N, format='csc')
    E = spsp.vstack((E[:m1,:],E[m1-1,:],E[m1:,:]))
    
    H_plus = spsp.bmat([[SBP_1.H,None],[None,SBP_2.H]],format='csc')
    D1_plus = spsp.bmat([[SBP_1.D1,None],[None,SBP_2.D1]],format='csc')
    
    H = E.T@H_plus@E
    H = H.todia()
    HI = spsp.diags(1/H.data[0])
    D1 = HI@E.T@H_plus@D1_plus@E
        
    e_l = np.zeros(N)
    e_l[0] = 1
    e_r = np.zeros(N)
    e_r[-1] = 1
    
    e_l = spsp.csc_array(e_l)
    e_r = spsp.csc_array(e_r)
    
    xvec = np.concatenate((SBP_1.xvec[:-1],SBP_2.xvec),0)
    h = min(SBP_1.h,SBP_2.h)
    
    return SBP_1D(H,HI,D1,xvec,e_l,e_r,h,N)

# Build 2D SBP operators from two 1D SBP operators
# Input:
#   SBPx - object containing the operators in the x-direction
#   SBPy - object containing the operators in the y-direction
# 
# Output:
# SBP_2D - object containing the 2D operators
def build_2D_SBP(SBPx,SBPy):
    Imx = spsp.eye(SBPx.m)
    Imy = spsp.eye(SBPy.m)
    
    e_w = spsp.kron(SBPx.e_l,Imy,format='csc')
    e_e = spsp.kron(SBPx.e_r,Imy,format='csc')
    e_s = spsp.kron(Imx,SBPy.e_l,format='csc')
    e_n = spsp.kron(Imx,SBPy.e_r,format='csc')
    
    Dx = spsp.kron(SBPx.D1,Imy)
    Dy = spsp.kron(Imx,SBPy.D1)

    H = spsp.kron(SBPx.H,SBPy.H)
    HI = spsp.kron(SBPx.HI,SBPy.HI)
    
    ones_xvec = np.ones(SBPx.m)
    ones_yvec = np.ones(SBPy.m)
    n_w = [-1*ones_yvec,0*ones_yvec]
    n_e = [1*ones_yvec,0*ones_yvec]
    n_s = [0*ones_xvec,-1*ones_xvec]
    n_n = [0*ones_xvec,1*ones_xvec]
    t_w = [-n_w[1],n_w[0]]
    t_e = [-n_e[1],n_e[0]]
    t_s = [-n_s[1],n_s[0]]
    t_n = [-n_n[1],n_n[0]]
    
    H_w = SBPy.H
    H_e = SBPy.H
    H_s = SBPx.H
    H_n = SBPx.H
    
    # Check SBP properties
    assert (np.max(np.max(np.abs(H@Dx + Dx.T@H - (e_w.T@H_w@(spsp.diags(n_w[0])@e_w) + e_e.T@H_e@(spsp.diags(n_e[0])@e_e) + e_s.T@H_s@(spsp.diags(n_s[0])@e_s) + e_n.T@H_n@(spsp.diags(n_n[0])@e_n))))) < 1e-12), "The SBP-property of Dx does not hold"
    assert (np.max(np.max(np.abs(H@Dy + Dy.T@H - (e_w.T@H_w@(spsp.diags(n_w[1])@e_w) + e_e.T@H_e@(spsp.diags(n_e[1])@e_e) + e_s.T@H_s@(spsp.diags(n_s[1])@e_s) + e_n.T@H_n@(spsp.diags(n_n[1])@e_n))))) < 1e-12), "The SBP-property of Dy does not hold"
    
    return SBP_2D(H,HI,Dx,Dy,H_w,H_e,H_s,H_n,e_w,e_e,e_s,e_n,n_w,n_e,n_s,n_n,t_w,t_e,t_s,t_n,SBPx.m*SBPy.m,[SBPx.m,SBPy.m])

# Build 2D curvilinear SBP operators from 2D SBP operators on cartesian grid
# Input:
#   SBP_ref - object containing the 2D operators of cartesian grid
#   coords - list [x,y] where x and y contains the grid points of curvilinear grid
# 
# Output:
# SBP_2D - object containing the 2D curvilinear SBP operators
def build_SBP_curv(SBP_ref,coords):
    x,y = coords
    x_u = SBP_ref.Dx*x
    x_v = SBP_ref.Dy*x
    y_u = SBP_ref.Dx*y
    y_v = SBP_ref.Dy*y

    J = x_u*y_v - x_v*y_u
    Ji = 1/J
    
    assert (all(J > 0)),"Non-positive Jacobian, probably bad mapping."

    Xu = spsp.diags(x_u)
    Xv = spsp.diags(x_v)
    Yu = spsp.diags(y_u)
    Yv = spsp.diags(y_v)
    
    e_w = SBP_ref.e_w
    e_e = SBP_ref.e_e
    e_s = SBP_ref.e_s
    e_n = SBP_ref.e_n
    
    # Scaling factors
    s_w = np.sqrt((e_w@x_v)**2 + (e_w@y_v)**2)
    s_e = np.sqrt((e_e@x_v)**2 + (e_e@y_v)**2)
    s_s = np.sqrt((e_s@x_u)**2 + (e_s@y_u)**2)
    s_n = np.sqrt((e_n@x_u)**2 + (e_n@y_u)**2)

    # x- and y-derivatives
    Dx = 0.5*spsp.diags(Ji)*(Yv@SBP_ref.Dx + SBP_ref.Dx@Yv - Yu@SBP_ref.Dy - SBP_ref.Dy@Yu)
    Dy = 0.5*spsp.diags(Ji)*(Xu@SBP_ref.Dy + SBP_ref.Dy@Xu - Xv@SBP_ref.Dx - SBP_ref.Dx@Xv)
    
    # Inner products
    H = spsp.diags(J)@SBP_ref.H
    H = H.todia()
    Hi = spsp.diags(1/H.data[0])

    H_w = SBP_ref.H_w@spsp.diags(s_w)
    H_e = SBP_ref.H_e@spsp.diags(s_e)
    H_s = SBP_ref.H_s@spsp.diags(s_s)
    H_n = SBP_ref.H_n@spsp.diags(s_n)

    # Normal and tangent vectors
    nu_w = [-1,0];
    nu_e = [1,0];
    nu_s = [0,-1];
    nu_n = [0,1];

    K11 = y_v/J
    K12 = -y_u/J
    K21 = -x_v/J
    K22 = y_v/J

    n_w_1 = spsp.diags(1/s_w)@e_w@(spsp.diags(J)@(K11*nu_w[0] + K12*nu_w[1]))
    n_w_2 = spsp.diags(1/s_w)@e_w@(spsp.diags(J)@(K21*nu_w[0] + K22*nu_w[1]))
    n_e_1 = spsp.diags(1/s_e)@e_e@(spsp.diags(J)@(K11*nu_e[0] + K12*nu_e[1]))
    n_e_2 = spsp.diags(1/s_e)@e_e@(spsp.diags(J)@(K21*nu_e[0] + K22*nu_e[1]))
    n_s_1 = spsp.diags(1/s_s)@e_s@(spsp.diags(J)@(K11*nu_s[0] + K12*nu_s[1]))
    n_s_2 = spsp.diags(1/s_s)@e_s@(spsp.diags(J)@(K21*nu_s[0] + K22*nu_s[1]))
    n_n_1 = spsp.diags(1/s_n)@e_n@(spsp.diags(J)@(K11*nu_n[0] + K12*nu_n[1]))
    n_n_2 = spsp.diags(1/s_n)@e_n@(spsp.diags(J)@(K21*nu_n[0] + K22*nu_n[1]))

    n_w = [n_w_1,n_w_2]
    n_e = [n_e_1,n_e_2]
    n_s = [n_s_1,n_s_2]
    n_n = [n_n_1,n_n_2]
    t_w = [-n_w_2,n_w_1]
    t_e = [-n_e_2,n_e_1]
    t_s = [-n_s_2,n_s_1]
    t_n = [-n_n_2,n_n_1]

    # Check SBP properties
    assert (np.max(np.max(np.abs(H@Dx + Dx.T@H - (e_w.T@H_w@(spsp.diags(n_w[0])@e_w) + e_e.T@H_e@(spsp.diags(n_e[0])@e_e) + e_s.T@H_s@(spsp.diags(n_s[0])@e_s) + e_n.T@H_n@(spsp.diags(n_n[0])@e_n))))) < 1e-12), "The SBP-property of Dx does not hold"
    assert (np.max(np.max(np.abs(H@Dy + Dy.T@H - (e_w.T@H_w@(spsp.diags(n_w[1])@e_w) + e_e.T@H_e@(spsp.diags(n_e[1])@e_e) + e_s.T@H_s@(spsp.diags(n_s[1])@e_s) + e_n.T@H_n@(spsp.diags(n_n[1])@e_n))))) < 1e-12), "The SBP-property of Dy does not hold"
    
    return SBP_2D(H,Hi,Dx,Dy,H_w,H_e,H_s,H_n,e_w,e_e,e_s,e_n,n_w,n_e,n_s,n_n,t_w,t_e,t_s,t_n,SBP_ref.N,SBP_ref.m)
    
# Construct coordinates of curvilinear grid
# Input:
#   m - list [mx,my] of grid points in each dimension of reference grid
#   xlims - list [xl,xr], boundaries in x-direction of reference domain
#   ylims - list [yl,yr], boundaries in y-direction of reference domain
# 
# Output:
# [x,y] - x and y are arrays containing the grid points of the curvilinear grid
def get_physical_coords(m,xlims,ylims):
    N = m[0]*m[1]
    xvec = np.linspace(xlims[0],xlims[1],m[0])
    yvec = np.linspace(ylims[0],ylims[1],m[1])
    [XI,ETA] = np.meshgrid(xvec,yvec)
    xi = XI.T.reshape(N)
    eta = ETA.T.reshape(N)
    
    # Offset reference domain boundaries by sin-function
    amp = 0.1
    q_x = amp*np.sin(2*np.pi*xvec)
    q_y = amp*np.sin(4*np.pi*yvec)
    
    X_left = xlims[0] + np.kron(np.ones(m[0]),q_y)
    X_right = xlims[1] + np.kron(np.ones(m[0]),q_y)
    Y_left = ylims[0] + np.kron(q_x,np.ones(m[1]))
    Y_right = ylims[1] + np.kron(q_x,np.ones(m[1]))
    
    # Physical grid coordinates
    x = X_left + (X_right - X_left)*(xi - xlims[0])/(xlims[1] - xlims[0])
    y = Y_left + (Y_right - Y_left)*(eta - ylims[0])/(ylims[1] - ylims[0])
    
    return [x,y]
