% Get indices of specified boundary in a 2D block.
function ids_bound = get_bound_idx(bound,mx,my)
switch bound
    case 'w'
        e_l = sparse(mx,1); e_l(1) = 1;
        ew = kr(e_l,speye(my));
        [ids_bound,~] = find(ew);
    case 'e'
        e_r = sparse(mx,1); e_r(end) = 1;
        ee = kr(e_r,speye(my));
        [ids_bound,~] = find(ee);
    case 's'
        e_l = sparse(my,1); e_l(1) = 1;
        es = kr(speye(mx),e_l);
        [ids_bound,~] = find(es);
    case 'n'
        e_r = sparse(my,1); e_r(end) = 1;
        en = kr(speye(mx),e_r);
        [ids_bound,~] = find(en);
end