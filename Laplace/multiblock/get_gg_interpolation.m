% Create interpolation operators between GL grid and equidistant grid.
function [Pv2u,Pu2v] = get_gg_interpolation(G_v,G_u,H_v,H_u,order_v,order_u)

% Vandermonde matrix
xtmp = 2*G_v.grids{1}.logic.x{1}-1;
V = legendrePolynomials(2/(xtmp(end) - xtmp(1))*(xtmp - xtmp(1))-1,order_v,0);

Nb_interf = length(G_v.boundaryGroups.interface);
Pf2b = inv(kron(eye(Nb_interf),V));

% Grid point distribution of GL grid 
bot_points = full(G_v.getBoundary(G_v.boundaryGroups.interface));
[~,I] = sort(bot_points(:,1)); % sort to make sure grid points are sequential
bot_points = bot_points(I,:);
xb_edges = full(bot_points([1:order_v+1:end,(order_v+1)*Nb_interf],1));
water_points = G_u.getBoundary('s');
water_points = water_points(:,1);
[~, ~, Pg2a, Pb2g, ~,~,~,~,~] = make_projection_g2g_hr(order_v, water_points, xb_edges);

G_water_m = G_u.size();
[~, Pag2f,~,~] = make_projection(G_water_m(1)-1,order_u);
Pag2f = Pag2f*kron(speye(G_water_m(1)-1),speye(order_u,order_v+1));

E = eye((order_v+1)*Nb_interf);
E = E(I,:);

% Build SBP preserving interpolation operators
Pv2u = Pag2f*Pg2a*Pb2g*Pf2b*E;
Pu2v = inv(E)*inv(H_v)*E*Pv2u'*H_u;