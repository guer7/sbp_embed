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
    e_w: spsp._arrays.csc_array
    e_e: spsp._arrays.csc_array
    e_s: spsp._arrays.csc_array
    e_n: spsp._arrays.csc_array
    xvec: np.ndarray
    yvec: np.ndarray
    SBPx: SBP_1D
    SBPy: SBP_1D
    N: int

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
    HI = spsplg.inv(H.tocsc())
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
    
    e_w = spsp.kron(SBPx.e_l,Imy)
    e_e = spsp.kron(SBPx.e_r,Imy)
    e_s = spsp.kron(Imx,SBPy.e_l)
    e_n = spsp.kron(Imx,SBPy.e_r)
    
    Dx = spsp.kron(SBPx.D1,Imy)
    Dy = spsp.kron(Imx,SBPy.D1)

    H = spsp.kron(SBPx.H,SBPy.H)
    HI = spsp.kron(SBPx.HI,SBPy.HI)
    return SBP_2D(H,HI,Dx,Dy,e_w,e_e,e_s,e_n,SBPx.xvec,SBPy.xvec,SBPx,SBPy,SBPx.m*SBPy.m)