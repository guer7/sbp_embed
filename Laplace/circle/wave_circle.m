% Solves the acoustic wave equation in 2D on a circle domain using SBP GL
% operators with the embedding method, boundary optimized operators with SBP-P-SAT, or
% traditional operators with SBP-P-SAT.
%
% Inputs:
% m - Number of grid points per block and dimension
% order - Order of spatial SBP operators. m = order+1 for the GL operators
% method - Discretization method, "embed", "boundOpt", or "trad"
% animate - If plotting movie or not
% varargin, h - Element size for SBP GL method
%
% Output:
%
% err - Relative L2 error
% dofs - Degrees of freedom
% time_elapsed - Runtime of time stepping loop
% nonzeros - Number of non-zeros in the RHS matrix
% Nblocks - Number of blocks
%
function [err,dofs,specrad,time_elapsed,nonzeros,Nblocks] = wave_circle(m,order,method,animate,varargin)

tic
addpath('../sbplib')
addpath('../utils/')

ts_order = 4; % time stepping order
c = 1; % wave speed
Tend = 0.8; % final time

% Point source and exact solution
sigma = 1/25;
t_s = 0.3;
f = @(t) 1/(sigma*sqrt(2*pi))*exp(-(t - t_s).^2/(2*sigma^2));
uexact = @(r,t) uexact_func(r,t,t_s,sigma);

% Create grid
switch method
    case "embed"
        assert(nargin == 5)
        h = varargin{1};
        assert(m == order+1)
        opSet = @sbp.D2GaussLob;
        mesh_filename = "../mesh_files/circ_h" + num2str(h) + "/circle.mesh";
        dom = multiblock.HOHQMesh(mesh_filename);
        G = dom.getGrid(m,'GL',order);
    case "boundOpt"
        opSet = @sbp.D2Nonequidistant;
        dom = multiblock.domain.Circle();
        G = dom.getGrid(m,'boundaryopt',order,'accurate');
    case "trad"
        gridType = 'equidist';
        opSet = @sbp.D2Variable;
        dom = multiblock.domain.Circle();
        G = dom.getGrid(m,'equidist');
end
Nblocks = G.nBlocks;
points = G.points;

% Create global Laplace operator, including SAT for continuity of fluxes.
lapl_plus = multiblock.DiffOp(@scheme.LaplaceCurvilinearNoPars,G,order,'hybrid',{opSet,1});

% Impose continuity of the solution across interfaces and Dirichlet
% boundary conditions using the projection method.
e_bound = lapl_plus.getBoundaryOperator('e',G.boundaryGroups.all);
switch method
    case "embed"
        lapl = build_embedded_ops(lapl_plus,G);
        L = lapl.bops.outer.e';
        H = lapl.H;
        DL = lapl.DL;
        E = lapl.E;
        points = lapl.points;
    otherwise
        L = [e_bound';lapl_plus.L];
        H = lapl_plus.H;
        DL = lapl_plus.D;
        E = speye(G.N);
end

HI = inv(H);
L = remove_deps(L,H);
Lplus = HI*L'*inv(L*HI*L');
P = speye(size(H,1)) - Lplus*L;
A = c^2*P*DL*P;
A = set2zero(A,1e-14,0);

nonzeros = nnz(A);
dofs = size(A,1);

% Build discrete forcing function
pstar_try = [0,0];
[~,pstar_idx] = min(sum(abs(points - pstar_try),2)); % find point closes to [0,0]
pstar = full(points(pstar_idx,:));
d = P*(get_unit_vec(pstar_idx,size(A,1))/H(pstar_idx,pstar_idx));
S = @(t) d*f(t);

fprintf("Preprocessing time: %f seconds\n",toc)

% Define initial data
v = zeros(dofs,1);
vt = zeros(dofs,1);
v = P*v;
vt = P*vt;

% Define timestepper
specrad = abs(eigs(A,1));
dt_try = 0.2*2.8/sqrt(specrad);
[dt,mt] = alignedTimestep(dt_try,Tend);
ts = time.RungekuttaSecondOrder(A,[],S,dt,0,v,vt,ts_order);

% Integrate in time
if animate
    figure
    xl = -1.2;
    xr = 1.2;
    yl = -1.2;
    yr = 1.2;
    zmin = -0.05;
    zmax = 1.01;
    vplot = E*v;
    ph = multiblock.Surface(G,vplot);
    ph.set('Edgecolor','none')
    caxis([zmin,zmax])
    axis([xl,xr,yl,yr,zmin,zmax])
    colorbar

    pause
    for tidx = 1:mt
        ts.step();
        v = ts.getV();
        if mod(tidx,20) == 0
            vplot = E*v;
            ph.ZData = vplot;
            ph.CData = vplot;
            drawnow
        end
        title(ts.t)
    end
else
    tic
    ts.stepN(mt,true);
    time_elapsed = toc;
    v = ts.getV();
end

% Compute error
vexact = uexact(sqrt((pstar(1) - points(:,1)).^2 + (pstar(2) - points(:,2)).^2),Tend);
err_diff = vexact - v;
err_diff(pstar_idx) = 0; % remove point source
err = norm(err_diff)/norm(vexact);