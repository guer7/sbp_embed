% Creates a multiblock DefCurvilinear object from HOHQMesh output file.
classdef HOHQMesh < multiblock.DefCurvilinear

    methods
        function obj = HOHQMesh(mesh_filename)
            [nNodes,nEdges,nBlocks,bound_order,coords,edges_data,element_data,domain_bound_names] = parse_hohqmesh(mesh_filename);

            obs = cell(nBlocks,4);
            conn = cell(nBlocks, nBlocks);
            bg_all = {};
            bg_sepbound = cell(4,nBlocks);

            bound_names = {'s','e','n','w'};

            for bidx = 1:nBlocks
                n1_idx = element_data{bidx}{1}(1);
                n2_idx = element_data{bidx}{1}(2);
                n3_idx = element_data{bidx}{1}(3);
                n4_idx = element_data{bidx}{1}(4);

                % n1 - n2 edge
                edge_data = get_edge_data(edges_data,n1_idx,n2_idx);
                if edge_data(4) == 0 % is physical boundary
                        assert(edge_data(3) == bidx)
                        curr_element = element_data{bidx};
                        assert(any(1 == find(curr_element{2})))
                        start_idx = sum(curr_element{2}(1:1-1))*(bound_order+1) + 1;
                        stop_idx = start_idx + bound_order;
                        [g_fun,g_fun_deriv] = get_lgl_interp(curr_element{3}(start_idx:stop_idx,:),bound_order);
                        obs{bidx,1} = parametrization.Curve(g_fun,g_fun_deriv);
                        bg_all{end+1} = {bidx,'s'};
                        bg_sepbound{1,bidx} = {bidx,'s',curr_element{4}{1}};
                else
                    obs{bidx,1} =  parametrization.Curve.line(coords(n1_idx,:), coords(n2_idx,:));
                    sign_me = sign(edge_data(5));
                    sign_nb = sign(edge_data(6));
                    conn{edge_data(3),edge_data(4)} = {bound_names{abs(edge_data(5))},bound_names{abs(edge_data(6))},sign_me,sign_nb};
                end
    
                % n2 - n3 edge
                edge_data = get_edge_data(edges_data,n2_idx,n3_idx);
                if edge_data(4) == 0 % is physical boundary
                        assert(edge_data(3) == bidx)
                        curr_element = element_data{bidx};
                        assert(any(2 == find(curr_element{2})))
                        start_idx = sum(curr_element{2}(1:2-1))*(bound_order+1) + 1;
                        stop_idx = start_idx + bound_order;
                        points = curr_element{3}(start_idx:stop_idx,:);
                        if norm(points(1,:)' - obs{bidx,1}.g(1)) > 1e-12
                            points = flip(points,1);
                        end
                        [g_fun,g_fun_deriv] = get_lgl_interp(points,bound_order);
                        obs{bidx,2} = parametrization.Curve(g_fun,g_fun_deriv);
                        bg_all{end+1} = {bidx,'e'};
                        bg_sepbound{2,bidx} = {bidx,'e',curr_element{4}{2}};
                else
                    obs{bidx,2} =  parametrization.Curve.line(coords(n2_idx,:), coords(n3_idx,:));
                    sign_me = sign(edge_data(5));
                    sign_nb = sign(edge_data(6));
                    conn{edge_data(3),edge_data(4)} = {bound_names{abs(edge_data(5))},bound_names{abs(edge_data(6))},sign_me,sign_nb};
                end

                % n3 - n4 edge
                edge_data = get_edge_data(edges_data,n3_idx,n4_idx);
                if edge_data(4) == 0 % is physical boundary
                        assert(edge_data(3) == bidx)
                        curr_element = element_data{bidx};
                        assert(any(3 == find(curr_element{2})))
                        start_idx = sum(curr_element{2}(1:3-1))*(bound_order+1) + 1;
                        stop_idx = start_idx + bound_order;
                        points = curr_element{3}(start_idx:stop_idx,:);
                        if norm(points(1,:)' - obs{bidx,2}.g(1)) > 1e-12
                            points = flip(points,1);
                        end
                        [g_fun,g_fun_deriv] = get_lgl_interp(points,bound_order);
                        obs{bidx,3} = parametrization.Curve(g_fun,g_fun_deriv);
                        bg_all{end+1} = {bidx,'n'};
                        bg_sepbound{3,bidx} = {bidx,'n',curr_element{4}{3}};
                else
                    obs{bidx,3} =  parametrization.Curve.line(coords(n3_idx,:), coords(n4_idx,:));
                    sign_me = sign(edge_data(5));
                    sign_nb = sign(edge_data(6));
                    conn{edge_data(3),edge_data(4)} = {bound_names{abs(edge_data(5))},bound_names{abs(edge_data(6))},sign_me,sign_nb};
                end

                % n4 - n1 edge
                edge_data = get_edge_data(edges_data,n4_idx,n1_idx);
                if edge_data(4) == 0 % is physical boundary
                        assert(edge_data(3) == bidx)
                        curr_element = element_data{bidx};
                        assert(any(4 == find(curr_element{2})))
                        start_idx = sum(curr_element{2}(1:4-1))*(bound_order+1) + 1;
                        stop_idx = start_idx + bound_order;
                        points = curr_element{3}(start_idx:stop_idx,:);
                        if norm(points(1,:)' - obs{bidx,3}.g(1)) > 1e-12
                            points = flip(points,1);
                        end
                        [g_fun,g_fun_deriv] = get_lgl_interp(points,bound_order);
                        obs{bidx,4} = parametrization.Curve(g_fun,g_fun_deriv);
                        bg_all{end+1} = {bidx,'w'};
                        bg_sepbound{4,bidx} = {bidx,'w',curr_element{4}{4}};
                else
                    obs{bidx,4} =  parametrization.Curve.line(coords(n4_idx,:), coords(n1_idx,:));
                    sign_me = sign(edge_data(5));
                    sign_nb = sign(edge_data(6));
                    conn{edge_data(3),edge_data(4)} = {bound_names{abs(edge_data(5))},bound_names{abs(edge_data(6))},sign_me,sign_nb};
                end

                blockNames{1, bidx} = num2str(bidx);
            end
            boundaryGroups = struct();
            

            bgs = struct();
            bgs.all = {};

            for idx = 1:numel(domain_bound_names)
                bgs.(domain_bound_names{idx}) = {};
            end

            for bidx = 1:nBlocks
                for boundidx = 1:4
                    bound = bg_sepbound{boundidx,bidx};
                    if isempty(bound)
                        continue
                    end
                    bgs.all{end+1} = bound(1:2);
                    bgs.(bound{3}){end+1} = bound(1:2);
                end
            end

            boundaryGroups.all = multiblock.BoundaryGroup(bgs.all);
            for idx = 1:numel(domain_bound_names)
                boundaryGroups.(domain_bound_names{idx}) = multiblock.BoundaryGroup(bgs.(domain_bound_names{idx}));
            end
            boundaryGroups.domain_bound_names = domain_bound_names;

            % create each block using transfinite interpolation
            for bidx = 1:nBlocks
                blocks{bidx} = parametrization.Ti(obs{bidx, 1}, obs{bidx, 2}, obs{bidx, 3}, obs{bidx, 4});
            end

            obj = obj@multiblock.DefCurvilinear(blocks, conn, boundaryGroups, blockNames);


        end

        function ms = getGridSizes(obj, m)
            for idx = 1:obj.nBlocks
                ms{idx} = [m,m];
            end
        end

    end
end
