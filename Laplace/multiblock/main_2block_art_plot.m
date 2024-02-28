% This script reproduces numerical results from FIXME.
% 
% Solves the acoustic wave equation in 2D on a complex domain. Lower block
% is discretized using SBP GL operators on quadrilaterals. The upper block
% is discretized using traditional SBP operators on one equdistiant block.
% The two blocks are coupled together using the SBP-P-SAT method with
% interpolation operators from
% 
% J.E. Kozdon, L.C. Wilcox, Stable coupling of nonconforming, high-order 
% finite difference methods, SIAM J. Sci. Comput. (2016) A923–A952, 
% https://doi.org/10.1137/15M1022823.
% 
clear
close all

tic

animate = 1;
ts_order = 4; % time stepping order

% Accuracy of spatial operators
order_bot = 5;
order_top = 4;

% Wave speeds
c_top = 1500;
c_bot = 1800;

% Grid points per block and dimension
m_top = 151;
m_bot = order_bot + 1;

addpath('../sbplib')
addpath('../utils/')
addpath('../interpolation/')

fprintf("Building rhs matrices...")
% Create grid and operators in top region
dom_top = multiblock.domain.Rectangle([0,100],[100,0]);
G_top_ = dom_top.getGrid([m_top,m_top]);
G_top = G_top_.grids{1};
lapl_top = scheme.LaplaceCurvilinearNoPars(G_top,order_top);
N_top = lapl_top.size;

% Create grid and operators in bot region
mesh_filename = "../mesh_files/seabottom_h4/seabottom.mesh";
dom_bot = multiblock.HOHQMesh(mesh_filename);
G_bot = dom_bot.getGrid(m_bot,'GL',order_bot);
lapl_bot_plus = multiblock.DiffOp(@scheme.LaplaceCurvilinearNoPars,G_bot,order_bot,'hybrid',{@sbp.D2GaussLob,1});
lapl_bot = build_embedded_ops(lapl_bot_plus,G_bot); % Build embedding operator
N_bot = lapl_bot.Ntot;

Ntot = N_top + N_bot;
points = [G_top.points;lapl_bot.points];

% Build interpolation operators
[Pb2w,Pw2b] = get_gg_interpolation(G_bot,G_top,lapl_bot.bops.interface.H,lapl_top.H_s,order_bot,order_top);

% Couple the two regions together using SBP-P-SAT and impose boundary conditons
SAT_IC = -[sparse(N_top,N_top),sparse(N_top,N_bot);
    c_top^2*lapl_bot.HI*lapl_bot.bops.interface.e*lapl_bot.bops.interface.H*Pw2b*lapl_top.d_s',c_bot^2*lapl_bot.HI*lapl_bot.bops.interface.e*lapl_bot.bops.interface.H*lapl_bot.bops.interface.d'];
DL = [c_top^2*lapl_top.D,sparse(N_top,N_bot);
    sparse(N_bot,N_top),c_bot^2*lapl_bot.DL] + SAT_IC;

H = [lapl_top.H,sparse(N_top,N_bot);
    sparse(N_bot,N_top),lapl_bot.H];
HI = inv(H);

SAT_BC_D_top = -lapl_top.Hi*lapl_top.e_w*lapl_top.H_w*lapl_top.d_w' - lapl_top.Hi*lapl_top.e_e*lapl_top.H_e*lapl_top.d_e';
SAT_BC_D_bot = -inv(lapl_bot.H)*lapl_bot.bops.outer.e*lapl_bot.bops.outer.H*lapl_bot.bops.outer.d' - inv(lapl_bot.H)*lapl_bot.bops.inner.e*lapl_bot.bops.inner.H*lapl_bot.bops.inner.d';
SAT_BC_D = [c_top^2*SAT_BC_D_top,sparse(N_top,N_bot);
    sparse(N_bot,N_top),c_bot^2*SAT_BC_D_bot];

SAT_BC_E_top = -lapl_top.Hi*lapl_top.e_w*lapl_top.H_w*lapl_top.e_w' - lapl_top.Hi*lapl_top.e_e*lapl_top.H_e*lapl_top.e_e';
SAT_BC_E_bot = -inv(lapl_bot.H)*lapl_bot.bops.outer.e*lapl_bot.bops.outer.H*lapl_bot.bops.outer.e';
SAT_BC_E = [c_top*SAT_BC_E_top,sparse(N_top,N_bot);
    sparse(N_bot,N_top),c_bot*SAT_BC_E_bot];

Lbc = [lapl_top.e_n',sparse(size(lapl_top.e_n',1),N_bot)];
Lic = [lapl_top.e_s',-Pb2w*lapl_bot.bops.interface.e'];

L = [Lbc;Lic];
L = remove_deps(L,H);

Lplus = HI*L'*inv(L*HI*L');
P = speye(Ntot) - Lplus*L;
P = set2zero(P,1e-14,0);

% RHS matrices
A = P*(DL + SAT_BC_D)*P;
A = set2zero(A,1e-14,0);
B = P*SAT_BC_E*P;
B = set2zero(B,1e-14,0);

fprintf(" done!\n")

% Build discrete forcing function
pstar_try = [50,50];
sigma = 1/200;
t_s = 0.03;
f = @(t) 1/(sigma*sqrt(2*pi))*exp(-(t - t_s).^2/(2*sigma^2));
[~,pstar_idx] = min(sum(abs(points - pstar_try),2)); % find point closes to pstar_try
d = P*(get_unit_vec(pstar_idx,size(A,1))/H(pstar_idx,pstar_idx));
S = @(t) d*f(t);

% initial data
v = zeros(Ntot,1);
vt = zeros(Ntot,1);
v = P*v;
vt = P*vt;

specrad = abs(eigs(A,1));
dt_try = 0.2*2.8/sqrt(specrad);
tvec_save = [0.05,0.075,0.1,0.2];
Tend = tvec_save(end);
[dt,mt] = alignedTimestep(dt_try,Tend);
ts = time.RungekuttaSecondOrder(A,B,S,dt,0,v,vt,ts_order);

% Integrate in time
if animate
    % find indices to save in subplots
    save_ids = round(tvec_save/dt)+1;

    figure('units','pixels','pos',[313          21         932        1298])
    
    sp = 1;
    axs{1} = subplot(2,2,1);

    cminmax = [-6e-6,6e-6;-4e-6,4e-6;-4e-6,4e-6;-1e-6,1e-6];
    
    hold on
    box on
    fs = 18;
    xl = 0;
    xr = 100;
    yl = -100;
    yr = 100;
    zl = -3;
    zr = 3;
    v_top = v(1:N_top);
    v_bot = v(N_top+1:end);
    ph1 = multiblock.Surface(G_bot,lapl_bot.E*v_bot);
    ph1.set('Edgecolor','none')

    ph2 = multiblock.Surface(G_top_,v_top);
    ph2.set('Edgecolor','none')

    caxis(cminmax(sp,:))
    axis([xl,xr,yl,yr,zl,zr])
    xlabel('$x$','interpreter','latex')
    ylabel('$y$','interpreter','latex')
    set(gca,'Fontsize',fs)
    colormap(turbo)
    colorbar
    tstr = sprintf("%d",round(1000*ts.t));
    title("$t = " + tstr + "$ ms",'interpreter','latex')
    drawnow

    % time stepping loop
    for tidx = 1:mt
        ts.step();
        v = ts.getV();
        if mod(tidx,round(0.001/dt)) == 0 || any(tidx+1 == save_ids) % update plot
            v_top = v(1:N_top);
            v_bot = v(N_top+1:end);
            ph1.ZData = lapl_bot.E*v_bot;
            ph1.CData = lapl_bot.E*v_bot;
            ph2.ZData = v_top;
            ph2.CData = v_top;
            drawnow
            if any(tidx+1 == save_ids) && tidx ~= mt
                sp = sp + 1;
                axs{sp} = subplot(2,2,sp);
                hold on
                box on
                caxis(cminmax(sp,:))
                axis([xl,xr,yl,yr,zl,zr])
                xlabel('$x$','interpreter','latex')
                ylabel('$y$','interpreter','latex')
                set(gca,'Fontsize',fs)
                colormap(turbo)
                colorbar

                ph1 = multiblock.Surface(G_bot,lapl_bot.E*v_bot);
                ph1.set('Edgecolor','none')

                ph2 = multiblock.Surface(G_top_,v_top);
                ph2.set('Edgecolor','none')
                drawnow
            end
        end
        tstr = sprintf("%d",round(1000*ts.t));
        title("$t = " + tstr + "$ ms",'interpreter','latex')
    end
else
    tic
    ts.stepN(mt,true);
    time_elapsed = toc;
    v = ts.getV();
end

axs{3}.Position = axs{3}.Position + [0,0.04,0,0];
axs{4}.Position = axs{4}.Position + [0,0.04,0,0];
% exportgraphics(gcf,'t08.pdf','ContentType','vector')
exportgraphics(gcf,"multiblock.png",'Resolution',800)
