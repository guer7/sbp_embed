#!/usr/bin/env python3
import numpy as np
import scipy.sparse.linalg as spsplg
import scipy.sparse as spsp
import matplotlib.pyplot as plt
from matplotlib import cm
from utils import get_1D_sbp, merge_blocks, build_2D_SBP, build_SBP_curv, get_physical_coords

# Simulates Maxwell's equations in 2D with 4 blocks.
# Input:
#   m - number of grid points per dimension and block (integer)
#   order - order of accuracy of SBP operators (integer)
#   animate - if plotting simulation and divergence over time (boolean)
# 
# Output:
# div_end - norm of error at final time
def run_maxwell(m,order,animate=False): 
    # m = 41
    # order = 4
    # animate = False
    # Final time
    T = 1
        
    # Block boundaries
    xlims = [-1,0,1]
    ylims = [-1,0,1]
    
    print("Running Maxwell simulation with m = %d and order = %d." % (m,order))
    
    # PDE coefficients
    A = np.array([[0,0,0],[0,0,-1],[0,-1,0]])
    B = np.array([[0,1,0],[1,0,0],[0,0,0]])
    
    # Analytical solutions
    b = 3
    c = 4
    def H_exact(x,y,t):
        return np.cos(b*x - t*np.sqrt(b**2 + c**2) + c*y)
    
    def Ex_exact(x,y,t):
        return -1/np.sqrt(b**2 + c**2)*c*np.cos(b*x - t*np.sqrt(b**2 + c**2) + c*y)
    
    def Ey_exact(x,y,t):
        return 1/np.sqrt(b**2 + c**2)*b*np.cos(b*x - t*np.sqrt(b**2 + c**2) + c*y)
    
    # 1D SBP operators in each block
    SBP_left = get_1D_sbp(m,xlims[:2],order)
    SBP_right = get_1D_sbp(m,xlims[1:],order)
    SBP_bot = get_1D_sbp(m,ylims[:2],order)
    SBP_top = get_1D_sbp(m,xlims[1:],order)
    
    # 1D merged operators using the embedding method
    SBPx = merge_blocks(SBP_left,SBP_right);
    SBPy = merge_blocks(SBP_bot,SBP_top);
    
    # 2D SBP operators on reference grid
    SBP_ref = build_2D_SBP(SBPx,SBPy)
    
    # 2D SBP operators on physical grid
    x,y = get_physical_coords(SBP_ref.m,xlims[0:3:2],ylims[0:3:2])
    SBP = build_SBP_curv(SBP_ref,[x,y])
    
    I3 = np.eye(3)
    e1 = np.array([0,1,0])
    
    # Projection method for the boundary conditions
    # Remove corner points of e_w and e_e to avoid linear dependence
    e_bound = spsp.vstack((SBP.e_w[1:-1,:],SBP.e_e[1:-1,:],SBP.e_s,SBP.e_n),format='csc')
    L = spsp.kron(e1,e_bound)
    
    x_bound = e_bound@x
    y_bound = e_bound@y
    def g(t):
        return H_exact(x_bound,y_bound,t)
    
    H = spsp.kron(I3,SBP.H)
    HI = spsp.kron(I3,SBP.HI)
    Lplus = HI@L.T@spsplg.inv(spsp.csc_array((L@HI@L.T)))
    P = spsp.eye(3*SBP.N) - Lplus@L
    
    # RHS function
    Dtilde = spsp.kron(A,SBP.Dx) + spsp.kron(B,SBP.Dy)
    D = P@Dtilde@P
    S = P*Dtilde*Lplus
    def rhs(v,t):
        return D@v + S@g(t)
    
    # Time stepping
    dt_try = 0.1*SBP_left.h
    mt = int(np.ceil(T/dt_try) + 1)
    tvec,dt = np.linspace(0,T,mt,retstep=True)
    
    # Initial data
    H0 = H_exact(x,y,0)
    Ex0 = Ex_exact(x,y,0)
    Ey0 = Ey_exact(x,y,0)
    
    w = P@np.concatenate((Ex0,H0,Ey0),axis=0)
    D_DIV = spsp.hstack((SBP.Dx,SBP.Dy),format='csc')
    div = np.zeros(mt)
    E0 = np.concatenate((Ex0,Ey0),axis=0)
    div_curr = D_DIV@E0
    div[0] = np.sqrt(div_curr.T@SBP.H@div_curr)
    if animate:
        fig, ax = plt.subplots(nrows=1,ncols=3, subplot_kw={"projection": "3d"}, figsize=(20,6))
        ax[0].set_zlim(-1,1)
        ax[0].set_xlabel("x")
        ax[0].set_ylabel("y")
        ax[0].title.set_text("Ex, t = " + str(0))
        ax[1].set_zlim(-1,1)
        ax[1].set_xlabel("x")
        ax[1].set_ylabel("y")
        ax[1].title.set_text("H")
        ax[2].set_zlim(-1,1)
        ax[2].set_xlabel("x")
        ax[2].set_ylabel("y")
        ax[2].title.set_text("Ey")
        
        X = np.reshape(x,SBP.m).T
        Y = np.reshape(y,SBP.m).T
        
        srf0 = ax[0].plot_surface(X,Y,Ex0.reshape(SBP.m).T,cmap=cm.coolwarm)
        srf1 = ax[1].plot_surface(X,Y,H0.reshape(SBP.m).T,cmap=cm.coolwarm,vmin=-0.4,vmax=0.4)
        srf2 = ax[2].plot_surface(X,Y,Ey0.reshape(SBP.m).T,cmap=cm.coolwarm,vmin=-0.4,vmax=0.4)
        
        plt.draw()
    
    t = 0
    for tidx in range(mt-1):
        # RK4
        k1 = dt*rhs(w,t)
        k2 = dt*rhs(w + 0.5*k1,t + 0.5*dt)
        k3 = dt*rhs(w + 0.5*k2,t + 0.5*dt)
        k4 = dt*rhs(w + k3,t + dt)
        w = w + 1/6*(k1 + 2*k2 + 2*k3 + k4)
        t = t + dt
        
        v = w + Lplus@g(t)
        
        Ex = v[0:SBP.N]
        Hz = v[SBP.N:2*SBP.N]
        Ey = v[2*SBP.N:3*SBP.N]
        
        E = np.concatenate((Ex,Ey),axis=0)
        div_curr = D_DIV@E
        div[tidx+1] = np.sqrt(div_curr.T@SBP.H@div_curr)
        
        if animate and (tidx % 5 == 0 or tidx == mt-2):
            srf0.remove()
            srf1.remove()
            srf2.remove()
            srf0 = ax[0].plot_surface(X,Y,Ex.reshape(SBP.m),cmap=cm.coolwarm)
            srf1 = ax[1].plot_surface(X,Y,Hz.reshape(SBP.m),cmap=cm.coolwarm,vmin=-0.4,vmax=0.4)
            srf2 = ax[2].plot_surface(X,Y,Ey.reshape(SBP.m),cmap=cm.coolwarm,vmin=-0.4,vmax=0.4)
            ax[0].title.set_text("Ex, t = {:.2f}".format(t))
            plt.draw()
            plt.pause(1e-2)
            
        if (tidx % round(mt/10) == 0):
            print("%d/%d time steps done." % (tidx+1,mt))
    
    if animate:        
        plt.figure()
        plt.xlabel("t")
        plt.ylabel("Divergence")
        plt.semilogy(tvec,div)
        plt.show()
    
    
    v_exact = np.concatenate((Ex_exact(x,y,T),H_exact(x,y,T),Ey_exact(x,y,T)),axis=0)
    
    err_diff = v - v_exact
    err = np.sqrt(err_diff.T@H@err_diff)
    return err


if __name__ == "__main__":
    mvec = np.array([21,61,101])
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
    