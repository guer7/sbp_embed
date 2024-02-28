% Build embedded Laplace operators
function lapl = build_embedded_ops(lapl_plus,G)

[E,idx_remove] = build_E(G);
Ntot = size(E,2);
points = G.points;
points(idx_remove,:) = [];
points = full(points);

% bops.e = struct();

for bidx = 1:numel(G.boundaryGroups.domain_bound_names)
    bound_name = G.boundaryGroups.domain_bound_names{bidx};
    bops.(bound_name).e = E'*lapl_plus.getBoundaryOperator('e',G.boundaryGroups.(bound_name));
    bops.(bound_name).d = E'*lapl_plus.getBoundaryOperator('d',G.boundaryGroups.(bound_name));
    bops.(bound_name).H = lapl_plus.getBoundaryQuadrature(G.boundaryGroups.(bound_name));
end

H = E'*lapl_plus.H*E;
HI = inv(H);
DL = inv(H)*E'*lapl_plus.H*lapl_plus.D*E;

Dx = inv(H)*E'*lapl_plus.H*lapl_plus.Dx*E;
Dy = inv(H)*E'*lapl_plus.H*lapl_plus.Dy*E;

lapl.H = H;
lapl.HI = HI;
lapl.DL = DL;
lapl.Dx = Dx;
lapl.Dy = Dy;
lapl.bops = bops;
lapl.points = points;
lapl.E = E;
lapl.idx_remove = idx_remove;
lapl.Ntot = Ntot;