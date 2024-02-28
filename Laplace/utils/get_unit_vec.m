function ei = get_unit_vec(idx,m)
ei = sparse(m,1);
ei(idx) = 1;
end