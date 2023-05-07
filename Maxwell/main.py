#!/usr/bin/env python3
import numpy as np
import scipy.sparse.linalg as spsplg
import scipy.sparse as spsp
import matplotlib.pyplot as plt
from matplotlib import cm
from utils import get_1D_sbp, merge_blocks, build_2D_SBP

# Simulates Maxwell's equations in 2D with 4 blocks.
# Input:
#   m - number of grid points per dimension and block (integer)
#   order - order of accuracy of SBP operators (integer)
#   animate - if plotting simulation and divergence over time (boolean)
# 
# Output:
# div_end - norm of divergence at final time
def run_maxwell(m,order,animate=False):
    
    # Final time
    T = 1 
        
    # Block boundaries
    xl = -1
    xi = 0
    xr = 1
    yl = -1
    yi = 0
    yr = 1

    print("Running Maxwell simulation with m = %d and order = %d." % (m,order))
    
    # PDE coefficients
    A = np.array([[0,0,0],[0,0,-1],[0,-1,0]])
    B = np.array([[0,1,0],[1,0,0],[0,0,0]])
    
    eps1 = 1.0
    eps2 = 1.0
    def eps(x,y):
        return eps1 + (eps2 - eps1)*np.logical_and(y >= yi,x >= xi)
    
    mu1 = 1.0
    mu2 = 1.0
    def mu(x,y):
        return mu1 + (mu2 - mu1)*np.logical_and(y >= yi,x >= xi)
    
    # 1D SBP operators in each block
    SBP_left = get_1D_sbp(m,[xl,xi],order)
    SBP_right = get_1D_sbp(m,[xi,xr],order)
    SBP_bot = get_1D_sbp(m,[yl,yi],order)
    SBP_top = get_1D_sbp(m,[yi,yr],order)
    
    # 1D merged operators using the embedding method
    SBPx = merge_blocks(SBP_left,SBP_right);
    SBPy = merge_blocks(SBP_bot,SBP_top);
    
    # 2D SBP operators
    SBP = build_2D_SBP(SBPx,SBPy)
    mx = SBP.SBPx.m
    my = SBP.SBPy.m
    
    [X,Y] = np.meshgrid(SBP.xvec,SBP.yvec)
    x = X.T.reshape(SBP.N)
    y = Y.T.reshape(SBP.N)

    I3 = np.eye(3)
    e0 = I3[0,:]
    e2 = I3[2,:]
    
    # Projection method for the boundary conditions
    Lw = spsp.kron(e2,SBP.e_w)
    Le = spsp.kron(e2,SBP.e_e)
    Ls = spsp.kron(e0,SBP.e_s)
    Ln = spsp.kron(e0,SBP.e_n)
    L = spsp.vstack((Lw,Le,Ls,Ln),format='csc')
    CI = spsp.diags(np.concatenate((1./eps(x,y),1./mu(x,y),1./eps(x,y)),axis=0),format='csc')
    HI = spsp.kron(I3,SBP.HI)@CI
    P = spsp.eye(3*SBP.N) - HI@L.T@spsplg.inv(L@HI@L.T)@L
    
    # RHS operator
    D = P@CI@(spsp.kron(A,SBP.Dx) + spsp.kron(B,SBP.Dy))@P
    
    # Time stepping
    dt_try = 0.1*min(SBP.SBPx.h,SBP.SBPy.h)
    mt = int(np.ceil(T/dt_try) + 1)
    tvec,dt = np.linspace(0,T,mt,retstep=True)
    
    # Initial data
    rstar = 0.1
    x0 = 0
    y0 = 0
    H0 = np.exp(-(x - x0)**2/rstar**2 - (y - y0)**2/rstar**2)
    Ex0 = 0*H0
    Ey0 = 0*H0
    
    v = np.concatenate((Ex0,H0,Ey0),axis=0)
    D_DIV = spsp.hstack((SBP.Dx,SBP.Dy),format='csc')
    div = np.zeros(mt)
    if animate:
        fig, ax = plt.subplots(nrows=1,ncols=3, subplot_kw={"projection": "3d"}, figsize=(20,6))
        ax[0].set_zlim(-1,1)
        plt.xlabel("x")
        plt.ylabel("y")
        ax[1].set_zlim(-1,1)
        plt.xlabel("x")
        plt.ylabel("y")
        ax[2].set_zlim(-1,1)
        plt.xlabel("x")
        plt.ylabel("y")
        
        title = plt.title("t = " + str(0))
        srf0 = ax[0].plot_surface(X,Y,Ex0.reshape((my,mx)),cmap=cm.coolwarm,vmin=-0.4,vmax=0.4)
        srf1 = ax[1].plot_surface(X,Y,H0.reshape((my,mx)),cmap=cm.coolwarm,vmin=-0.4,vmax=0.4)
        srf2 = ax[2].plot_surface(X,Y,Ey0.reshape((my,mx)),cmap=cm.coolwarm,vmin=-0.4,vmax=0.4)
        plt.draw()
    
    t = 0
    for tidx in range(mt):
        # RK4
        k1 = dt*D@(v)
        k2 = dt*D@(v + 0.5*k1)
        k3 = dt*D@(v + 0.5*k2)
        k4 = dt*D@(v + k3)
        v = v + 1/6*(k1 + 2*k2 + 2*k3 + k4)
        t = t + dt
        
        Ex = v[0:SBP.N]
        Hz = v[SBP.N:2*SBP.N]
        Ey = v[2*SBP.N:3*SBP.N]
        
        E = np.concatenate((Ex,Ey),axis=0)
        div_curr = D_DIV@E
        div[tidx] = np.sqrt(div_curr.T@SBP.H@div_curr)
        
        if animate and (tidx % 5 == 0 or tidx == mt-1):
            srf0.remove()
            srf1.remove()
            srf2.remove()
            srf0 = ax[0].plot_surface(X,Y,Ex.reshape((my,mx)),cmap=cm.coolwarm,vmin=-0.4,vmax=0.4)
            srf1 = ax[1].plot_surface(X,Y,Hz.reshape((my,mx)),cmap=cm.coolwarm,vmin=-0.4,vmax=0.4)
            srf2 = ax[2].plot_surface(X,Y,Ey.reshape((my,mx)),cmap=cm.coolwarm,vmin=-0.4,vmax=0.4)
            title.set_text("t = {:.2f}".format(t))
            plt.draw()
            plt.pause(1e-2)
            
        if (tidx % round(mt/10) == 0):
            print("%d/%d time steps done." % (tidx+1,mt))
    
    if animate:        
        plt.figure()
        plt.semilogy(tvec,div)
        plt.show()
        
    return div[-1]


if __name__ == "__main__":
    mvec = np.array([21,51,101,201])
    order_vec = np.array([2,4,6])
    errvec = np.zeros((mvec.size,order_vec.size))
    for order_idx,order in enumerate(order_vec):
        for m_idx,m in enumerate(mvec):
            errvec[m_idx,order_idx] = run_maxwell(m,order)
        
    q = np.zeros((mvec.size,order_vec.size))
    for order_idx,order in enumerate(order_vec):
        for m_idx,m in enumerate(mvec[1:]):
            q[m_idx+1,order_idx] = -np.log(errvec[m_idx,order_idx]/errvec[m_idx+1,order_idx])/np.log(mvec[m_idx]/mvec[m_idx+1])
    
    
    for order_idx,order in enumerate(order_vec):
        print("--- Order: %d ---" % order)
        print("m\terr\t\tq")
        for idx in range(mvec.size):
            print("%d\t%f\t%f" % (mvec[idx],np.log10(errvec[idx,order_idx]),q[idx,order_idx]))
    