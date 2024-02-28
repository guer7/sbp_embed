% Constructs embedding operator E and array of removed indices from "plus"
% vector from a grid object.
function [E,idx_remove] =  build_E(G)

N = G.grids{1}.size;
N = N(1);

Ntot = G.N;
Nblocks = G.nBlocks;
Npb = N*N;

% construct cell array of all points on the same coordinate
coupled_points = cell(Ntot,1);
for bidx = 1:Nblocks
    for nbidx = 1:Nblocks
        intf = G.connections{bidx,nbidx};
        if isempty(intf)
            continue
        end
        b_ids = Npb*(bidx-1) + get_bound_idx(intf{1},N,N);
        nb_ids = Npb*(nbidx-1) + get_bound_idx(intf{2},N,N);
        if intf{3} ~= intf{4}
            nb_ids = flip(nb_ids);
        end
        for idx = 1:numel(b_ids)
            coupled_points{b_ids(idx)}(end+1) = nb_ids(idx);
            coupled_points{nb_ids(idx)}(end+1) = b_ids(idx);
        end
    end
end

% gather all multiconnected points (repeat to get all of neighbours neighbours)
for idx = 1:Ntot
    nb_points = coupled_points{idx};
    for nb_idx = 1:numel(nb_points)
        nbnb_points = coupled_points{nb_points(nb_idx)};
        for nbnb_idx = 1:numel(nbnb_points)
            coupled_points{idx}(end+1) = nbnb_points(nbnb_idx);
        end
    end
    coupled_points{idx} = unique(coupled_points{idx});
end
for idx = 1:Ntot
    nb_points = coupled_points{idx};
    for nb_idx = 1:numel(nb_points)
        nbnb_points = coupled_points{nb_points(nb_idx)};
        for nbnb_idx = 1:numel(nbnb_points)
            coupled_points{idx}(end+1) = nbnb_points(nbnb_idx);
        end
    end
    coupled_points{idx} = unique(coupled_points{idx});
end

% remove duplicates
for idx = 1:Ntot
    nb_points = coupled_points{idx};
    if isempty(nb_points)
        continue
    end
    for nb_idx = 2:numel(nb_points)
        coupled_points{nb_points(nb_idx)} = [];
    end

end

% remove empties
coupled_points_tmp = {};
count = 1;
for idx = 1:Ntot
    nb_points = coupled_points{idx};
    if isempty(nb_points)
        continue
    else
        coupled_points_tmp{count} = nb_points;
        count = count + 1;
    end
end
coupled_points = coupled_points_tmp';

% construct list of indices that are removed
idx_mapping = (1:Ntot)';
points = G.points;

N_remove = 0;
idx_remove = [];
for cidx = 1:numel(coupled_points)
    for idx = coupled_points{cidx}(2:end)
        idx_mapping(idx) = coupled_points{cidx}(1);
        idx_remove = [idx_remove;idx];
        N_remove = N_remove + 1;
    end
    % make sure the indices corresponds to the same grid point
    if abs(sum(sum(diff(full(points(coupled_points{cidx},:)),1))))/abs(sum(sum(abs(full(points(coupled_points{cidx},:)))))) > 1e-14
        error("not good")
    end
end
N_small = Ntot - N_remove;

points_small = points;
points_small(idx_remove,:) = [];

idx_remove = flip(sort(idx_remove));
for idx = 1:numel(idx_remove)
    idx_mapping(idx_mapping > idx_remove(idx)) = idx_mapping(idx_mapping > idx_remove(idx)) - 1;
end

% Construct E. With v on the small grid, v_plus = v(idx_mapping) = E*v.
E = sparse(Ntot,N_small);
for idx = 1:Ntot
    E(idx,idx_mapping(idx)) = 1;
end