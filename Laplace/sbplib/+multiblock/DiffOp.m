classdef DiffOp < scheme.Scheme
    properties
        grid
        order
        diffOps
        D
        Dx
        Dy
        H
        Hi
        A
        L
        Lcell
        Pic

        blockmatrixDiv
    end

    methods
        function obj = DiffOp(doHand, g, order, interface_method, doParam, intfTypes)
            %  doHand -- may either be a function handle or a cell array of
            %            function handles for each grid. The function handle(s)
            %            should be on the form do = doHand(grid, order, ...)
            %            Additional parameters for each doHand may be provided in
            %            the doParam input.
            %       g -- a multiblock grid
            %   order -- integer specifying the order of accuracy
            % doParam -- may either be a cell array or a cell array of cell arrays
            %            for each block. If it is a cell array with length equal
            %            to the number of blocks then each element is sent to the
            %            corresponding function handle as extra parameters:
            %            doHand(..., doParam{i}{:}) Otherwise doParam is sent as
            %            extra parameters to all doHand: doHand(..., doParam{:})
            %
            % intfTypes (optional) -- nBlocks x nBlocks cell array of types for
            %                                 every interface.
            default_arg('interface_method', 'sat')
            default_arg('doParam', [])
            default_arg('intfTypes', cell(g.nBlocks(), g.nBlocks()) );

            [getHand, getParam] = parseInput(doHand, g, doParam);

            obj.order = order;
            nBlocks = g.nBlocks();

            % Create the diffOps for each block
            obj.diffOps = cell(1, nBlocks);
            for i = 1:nBlocks
                h = getHand(i);
                p = getParam(i);
                if ~iscell(p)
                    p = {p};
                end
                obj.diffOps{i} = h(g.grids{i}, order, p{:});
            end


            % Build the norm matrix
            H = cell(nBlocks, nBlocks);
            for i = 1:nBlocks
                H{i,i} = obj.diffOps{i}.H;
            end

            Dx = cell(nBlocks, nBlocks);
            for i = 1:nBlocks
                Dx{i,i} = obj.diffOps{i}.Dx;
            end

            Dy = cell(nBlocks, nBlocks);
            for i = 1:nBlocks
                Dy{i,i} = obj.diffOps{i}.Dy;
            end
            
%             obj.H = blockmatrix.toMatrix(H);
%             obj.Hi = inv(obj.H);

            % Build the coefficient matrix
%             A = cell(nBlocks, nBlocks);
%             for i = 1:nBlocks
%                 A{i,i} = spdiag(obj.diffOps{i}.a);
%             end
%             obj.A = blockmatrix.toMatrix(A);

            % Build the differentiation matrix
            Ns = zeros(nBlocks,1);
            for i = 1:nBlocks
                Ns(i) = length(obj.diffOps{i}.D);
            end
            obj.blockmatrixDiv = {Ns, Ns};
            D = blockmatrix.zero(obj.blockmatrixDiv);
            for i = 1:nBlocks
                D{i,i} = obj.diffOps{i}.D;
            end

            L_row_idx = 1;
            for i = 1:nBlocks
                for j = 1:nBlocks
                    intf = g.connections{i,j};
                    if isempty(intf)
                        continue
                    end
                    
                    [ii, ij] = obj.diffOps{i}.interface(intf{1}, obj.diffOps{j}, intf{2}, intfTypes{i,j},interface_method);
                    D{i,i} = D{i,i} + ii;
                    D{i,j} = D{i,j} + ij;

                    [jj, ji] = obj.diffOps{j}.interface(intf{2}, obj.diffOps{i}, intf{1}, intfTypes{i,j},interface_method);
                    D{j,j} = D{j,j} + jj;
                    D{j,i} = D{j,i} + ji;

                    switch interface_method
                        case 'hybrid'
                            L{L_row_idx,i} = obj.diffOps{i}.getBoundaryOperator('e',intf{1})';

                            % check orientation of grid points
                            my_bound_points = obj.diffOps{i}.getBoundaryOperator('e',intf{1})'*obj.diffOps{i}.grid.points;
                            neigh_bound_points = obj.diffOps{j}.getBoundaryOperator('e',intf{2})'*obj.diffOps{j}.grid.points;
                            if norm(my_bound_points - neigh_bound_points) > 1e-10
                                L{L_row_idx,j} = flip(-obj.diffOps{j}.getBoundaryOperator('e',intf{2})',1);
                            else
                                L{L_row_idx,j} = -obj.diffOps{j}.getBoundaryOperator('e',intf{2})';
                            end
%                             fprintf("i: %d, j: %d\n",i,j);
                            L_row_idx = L_row_idx + 1;
                        case 'projection'
                            L{L_row_idx,i} = obj.diffOps{i}.getBoundaryOperator('e',intf{1})';
                            L{L_row_idx,j} = -obj.diffOps{j}.getBoundaryOperator('e',intf{2})';
                            L_row_idx = L_row_idx + 1;
                            L{L_row_idx,i} = obj.diffOps{i}.a*obj.diffOps{i}.getBoundaryOperator('d',intf{1})';
                            L{L_row_idx,j} = obj.diffOps{j}.a*obj.diffOps{j}.getBoundaryOperator('d',intf{2})';
                            L_row_idx = L_row_idx + 1;
                        case {'sat','satStandard','satMartins'}
                            L = {};
                    end
                end
            end
%             obj.D = blockmatrix.toMatrix(D);
            obj.grid = g;

            cm_div = cumsum(obj.blockmatrixDiv{1});
            start_idx = [0;cm_div(1:end-1)] + 1;
            stop_idx = cm_div;
% 2
            obj.H = sparse(obj.size);
            obj.D = sparse(obj.size);
            obj.Dx = sparse(obj.size);
            obj.Dy = sparse(obj.size);
            for i = 1:g.nBlocks
                for j = 1:g.nBlocks
                    % H
                    if ~isempty(H{i,j})
                        obj.H(start_idx(i):stop_idx(i),start_idx(j):stop_idx(j)) = H{i,j};
                    end
                    if ~isempty(D{i,j})
                        Dtmp = D{i,j};
                        obj.D(start_idx(i):stop_idx(i),start_idx(j):stop_idx(j)) = Dtmp;
                    end
                    if ~isempty(Dx{i,j})
                        Dxtmp = Dx{i,j};
                        obj.Dx(start_idx(i):stop_idx(i),start_idx(j):stop_idx(j)) = Dxtmp;
                    end
                    if ~isempty(Dy{i,j})
                        Dytmp = Dy{i,j};
                        obj.Dy(start_idx(i):stop_idx(i),start_idx(j):stop_idx(j)) = Dytmp;
                    end
                end
            end
            obj.Hi = inv(obj.H);

            obj.Lcell = L;
            obj.L = blockmatrix.toMatrix(L);
            function [getHand, getParam] = parseInput(doHand, g, doParam)
                if ~isa(g, 'multiblock.Grid')
                    error('multiblock:DiffOp:DiffOp:InvalidGrid', 'Requires a multiblock grid.');
                end

                if iscell(doHand) && length(doHand) == g.nBlocks()
                    getHand = @(i)doHand{i};
                elseif isa(doHand, 'function_handle')
                    getHand = @(i)doHand;
                else
                    error('multiblock:DiffOp:DiffOp:InvalidGridDoHand', 'doHand must be a function handle or a cell array of length grid.nBlocks');
                end

                if isempty(doParam)
                    getParam = @(i){};
                    return
                end

                if ~iscell(doParam)
                    getParam = @(i)doParam;
                    return
                end

                % doParam is a non-empty cell-array

                if length(doParam) == g.nBlocks() && all(cellfun(@iscell, doParam))
                    % doParam is a cell-array of cell-arrays
                    getParam = @(i)doParam{i};
                    return
                end

                getParam = @(i)doParam;
            end
        end

        function ops = splitOp(obj, op)
            % Splits a matrix operator into a cell-matrix of matrix operators for
            % each grid.
            ops = sparse2cell(op, obj.NNN);
        end

        % Get a boundary operator specified by opName for the given boundary/BoundaryGroup
        function op = getBoundaryOperator(obj, opName, boundary)
            switch class(boundary)
                case 'cell'
                    blockId = boundary{1};
                    localOp = obj.diffOps{blockId}.getBoundaryOperator(opName, boundary{2});

                    div = {obj.blockmatrixDiv{1}, size(localOp,2)};
                    blockOp = blockmatrix.zero(div);
                    blockOp{blockId,1} = localOp;
                    op = blockmatrix.toMatrix(blockOp);
                    return
                case 'multiblock.BoundaryGroup'
                    op = sparse(size(obj.D,1),0);
                    for i = 1:length(boundary)
                        op = [op, obj.getBoundaryOperator(opName, boundary{i})];
                    end
                otherwise
                    error('Unknown boundary indentifier')
            end
        end

        function op = getBoundaryQuadrature(obj, boundary)
            switch class(boundary)
                case 'cell'
                    blockId = boundary{1};
                    op = obj.diffOps{blockId}.getBoundaryQuadrature(boundary{2});
                    return
                case 'multiblock.BoundaryGroup'
                    N = length(boundary);
                    H_bm = cell(N,N);
                    for i = 1:N
                        H_bm{i,i} = obj.getBoundaryQuadrature(boundary{i});
                    end
                    op = blockmatrix.toMatrix(H_bm);
                otherwise
                    error('Unknown boundary indentifier')
            end
        end

        function op = getBoundaryNormal(obj, opName, boundary)
            switch class(boundary)
                case 'cell'
                    blockId = boundary{1};
                    op = obj.diffOps{blockId}.getBoundaryNormal(opName, boundary{2});
                    return
                case 'multiblock.BoundaryGroup'
                    N = length(boundary);
                    N_bm = cell(N, N);
                    for i = 1:N
                        H_bm{i,i} = obj.getBoundaryNormal(opName, boundary{i});
                    end
                    op = blockmatrix.toMatrix(H_bm);
                otherwise
                    error('Unknown boundary indentifier')
            end
        end

        % Creates the closure and penalty matrix for a given boundary condition,
        %    boundary -- the name of the boundary on the form {id,name} where
        %                id is the number of a block and name is the name of a
        %                boundary of that block example: {1,'s'} or {3,'w'}. It
        %                can also be a boundary group
        function [closure, penalty] = boundary_condition(obj, boundary, type)
            switch class(boundary)
                case 'cell'
                    [closure, penalty] = obj.singleBoundaryCondition(boundary, type);
                case 'multiblock.BoundaryGroup'
                    [n,m] = size(obj.D);
                    closure = {sparse(n,m),sparse(n,m)};
                    penalty = sparse(n,0);
                    for i = 1:length(boundary)
                        [closurePart, penaltyPart] = obj.boundary_condition(boundary{i}, type);
                        closure{1} = closure{1} + closurePart{1};
                        closure{2} = closure{2} + closurePart{2};
                        penalty = [penalty, penaltyPart];
                    end
                otherwise
                    error('Unknown boundary indentifier')
            end

        end

        function [closure, penalty] = singleBoundaryCondition(obj, boundary, type)
            I = boundary{1};
            name = boundary{2};
            % Get the closure and penaly matrices
            [blockClosure, blockPenalty] = obj.diffOps{I}.boundary_condition(name, type);

            % Expand to matrix for full domain.
%             closure{1} = multiblock.local2globalClosure(blockClosure{1}, obj.blockmatrixDiv, I);
%             closure{2} = multiblock.local2globalClosure(blockClosure{2}, obj.blockmatrixDiv, I);

            
            cm_div = cumsum(obj.blockmatrixDiv{1});
            start_idx = [0;cm_div(1:end-1)] + 1;
            stop_idx = cm_div;
            closure{1} = sparse(obj.size,obj.size);
            closure{2} = sparse(obj.size,obj.size);
            closure{1}(start_idx(I):stop_idx(I),start_idx(I):stop_idx(I)) = blockClosure{1};
            closure{2}(start_idx(I):stop_idx(I),start_idx(I):stop_idx(I)) = blockClosure{2};


            penalty = multiblock.local2globalPenalty(blockPenalty, obj.blockmatrixDiv, I);
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            error('not implemented')
        end

        % Size returns the number of degrees of freedom
        function N = size(obj)
            N = 0;
            for i = 1:length(obj.diffOps)
                N = N + obj.diffOps{i}.size();
            end
        end
    end
end