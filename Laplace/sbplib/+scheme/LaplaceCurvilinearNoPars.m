classdef LaplaceCurvilinearNoPars < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        m_tot
        h % Grid spacing
        dim % Number of spatial dimensions

        grid

        order % Order of accuracy for the approximation

        a

        % Inner products and operators for physical coordinates
        D % Laplace operator
        H, Hi % Inner product
        e_w, e_e, e_s, e_n
        d_w, d_e, d_s, d_n % Normal derivatives at the boundary
        H_w, H_e, H_s, H_n % Boundary inner products
        Dx, Dy % Physical derivatives
        M % Gradient inner product
        Pic

        n_w, n_e, n_s, n_n % Physical normals
        n1_w, n1_e, n1_s, n1_n
        n2_w, n2_e, n2_s, n2_n
        t1_w, t1_e, t1_s, t1_n
        t2_w, t2_e, t2_s, t2_n
        tangent_w, tangent_e, tangent_s, tangent_n % Physical tangents

        debugging

        % Metric coefficients
        J, Ji
        a11, a12, a22
        K
        x_u
        x_v
        y_u
        y_v
        s_w, s_e, s_s, s_n % Boundary integral scale factors

        % Inner product and operators for logical coordinates
        H_u, H_v % Norms in the x and y directions
        Hi_u, Hi_v
        Hu,Hv % Kroneckerd norms. 1'*Hx*v corresponds to integration in the x dir.
        Hiu, Hiv
        Du, Dv
        D2_u, D2_v
        Duu,Duv,Dvu,Dvv
        du_w, dv_w
        du_e, dv_e
        du_s, dv_s
        du_n, dv_n

        % Borrowing constants
        theta_M_u, theta_M_v
        theta_R_u, theta_R_v
        theta_H_u, theta_H_v

        % Temporary, only used for nonconforming interfaces but should be removed.
        lambda

        % Number of boundary points in minumum check
        bp

        % DERIVATIVES
        dd_w_dp, dH_w_dp
        dd_e_dp, dH_e_dp
        dd_s_dp, dH_s_dp
        dd_n_dp, dH_n_dp

        dlambda_dp
        dH_dp
        dHi_dp
        dD_dp
        dJ_dp
    end

    methods
        % Implements  a*div(b*grad(u)) as a SBP scheme
        % TODO: Implement proper H, it should be the real physical quadrature, the logic quadrature may be but in a separate variable (H_logic?)

        function obj = LaplaceCurvilinearNoPars(g, order, opSet, a)
            default_arg('opSet',@sbp.D2Variable);
            default_arg('a', 1);

            obj.debugging = false;

            % assert(isa(g, 'grid.Curvilinear'))
            if isa(a, 'function_handle')
                a = grid.evalOn(g, a);
                a = spdiag(a);
            end

            % Number of boundary points in minimum check
            switch order
                case 2
                    obj.bp = 2;
                case 4
                    obj.bp = 4;
                case 6
                    obj.bp = 7;
            end

            dim = 2;
            m = g.size();
            m_u = m(1);
            m_v = m(2);
            m_tot = g.N();

            % 1D operators
            ops_u = opSet(m_u, {0, 1}, order);
            ops_v = opSet(m_v, {0, 1}, order);

            h_u = ops_u.h;
            h_v = ops_v.h;

            I_u = speye(m_u);
            I_v = speye(m_v);

            D1_u = ops_u.D1;
            obj.D2_u = ops_u.D2;
            H_u =  ops_u.H;
            Hi_u = ops_u.HI;
            e_l_u = ops_u.e_l;
            e_r_u = ops_u.e_r;
            d1_l_u = ops_u.d1_l;
            d1_r_u = ops_u.d1_r;

            D1_v = ops_v.D1;
            obj.D2_v = ops_v.D2;
            H_v =  ops_v.H;
            Hi_v = ops_v.HI;
            e_l_v = ops_v.e_l;
            e_r_v = ops_v.e_r;
            d1_l_v = ops_v.d1_l;
            d1_r_v = ops_v.d1_r;

            % Misc.
            obj.m = m;
            obj.h = [h_u h_v];
            obj.order = order;
            obj.grid = g;
            obj.dim = dim;

            % Logical operators
            Du = kr(D1_u,I_v);
            Dv = kr(I_u,D1_v);
            obj.Du = Du;
            obj.Dv = Dv;
            obj.Hu  = kr(H_u,I_v);
            obj.Hv  = kr(I_u,H_v);
            obj.Hiu = kr(Hi_u,I_v);
            obj.Hiv = kr(I_u,Hi_v);
            obj.H_u = H_u;
            obj.H_v = H_v;
            obj.Hi_u = Hi_u;
            obj.Hi_v = Hi_v;

            e_w  = kr(e_l_u,I_v);
            e_e  = kr(e_r_u,I_v);
            e_s  = kr(I_u,e_l_v);
            e_n  = kr(I_u,e_r_v);
            obj.e_w = e_w;
            obj.e_e = e_e;
            obj.e_s = e_s;
            obj.e_n = e_n;
            obj.du_w = kr(d1_l_u,I_v);
            obj.dv_w = (e_w'*Dv)';
            obj.du_e = kr(d1_r_u,I_v);
            obj.dv_e = (e_e'*Dv)';
            obj.du_s = (e_s'*Du)';
            obj.dv_s = kr(I_u,d1_l_v);
            obj.du_n = (e_n'*Du)';
            obj.dv_n = kr(I_u,d1_r_v);

            obj = obj.build_D(g.points());
        end

        function [Duu,Dvv] = build_D2_kron(obj,ind,a1,a2)
            Duu = sparse(obj.grid.N);
            Dvv = sparse(obj.grid.N);

            for i = 1:obj.m(2)
                D = obj.D2_u(a1(ind(:,i)));
                p = ind(:,i);
                Duu(p,p) = D;
            end

            for i = 1:obj.m(1)
                D = obj.D2_v(a2(ind(i,:)));
                p = ind(i,:);
                Dvv(p,p) = D;
            end
        end

        function obj = build_D(obj,coords)
            x = coords(:,1);
            y = coords(:,2);

            x_u = obj.Du*x;
            x_v = obj.Dv*x;
            y_u = obj.Du*y;
            y_v = obj.Dv*y;

            J = x_u.*y_v - x_v.*y_u;
            a11 =  1./J .* (x_v.^2  + y_v.^2);
            a12 = -1./J .* (x_u.*x_v + y_u.*y_v);
            a22 =  1./J .* (x_u.^2  + y_u.^2);
            lambda = 1/2 * (a11 + a22 - sqrt((a11-a22).^2 + 4*a12.^2));

            K = cell(obj.dim, obj.dim);
            K{1,1} = y_v./J;
            K{1,2} = (-y_u./J);
            K{2,1} = (-x_v./J);
            K{2,2} = (x_u./J);

            obj.x_u = x_u;
            obj.x_v = x_v;
            obj.y_u = y_u;
            obj.y_v = y_v;

            % Assemble full operators
            L_12 = spdiag(a12);
            Duv = obj.Du*L_12*obj.Dv;
            Dvu = obj.Dv*L_12*obj.Du;

            %             Duu = sparse(m_tot);
            %             Dvv = sparse(m_tot);
            ind = grid.funcToMatrix(obj.grid, 1:obj.grid.N);
            [Duu,Dvv] = obj.build_D2_kron(ind,a11,a22);

            %             for i = 1:m_v
            %                 b_a11 = b*a11;
            %                 D = D2_u(b_a11(ind(:,i)));
            %                 p = ind(:,i);
            %                 Duu(p,p) = D;
            %             end
            %
            %             for i = 1:m_u
            %                 b_a22 = b*a22;
            %                 D = D2_v(b_a22(ind(i,:)));
            %                 p = ind(i,:);
            %                 Dvv(p,p) = D;
            %             end


            % Physical operators
            obj.J = spdiag(J);
            obj.Ji = spdiag(1./J);

%             Dx = obj.Ji*(spdiag(y_v)*obj.Du - spdiag(y_u)*obj.Dv);
%             Dy = obj.Ji*(-spdiag(x_v)*obj.Du + spdiag(x_u)*obj.Dv);

            obj.D = obj.Ji*(Duu + Duv + Dvu + Dvv);
            obj.H = obj.J*kr(obj.H_u,obj.H_v);
            obj.Hi = obj.Ji*kr(obj.Hi_u,obj.Hi_v);

            obj.e_w = obj.e_w;
            obj.e_e = obj.e_e;
            obj.e_s = obj.e_s;
            obj.e_n = obj.e_n;

            %% normal derivatives
            I_w = ind(1,:);
            I_e = ind(end,:);
            I_s = ind(:,1);
            I_n = ind(:,end);

            a11_w = spdiag(a11(I_w));
            a12_w = spdiag(a12(I_w));
            a11_e = spdiag(a11(I_e));
            a12_e = spdiag(a12(I_e));
            a22_s = spdiag(a22(I_s));
            a12_s = spdiag(a12(I_s));
            a22_n = spdiag(a22(I_n));
            a12_n = spdiag(a12(I_n));

            s_w = sqrt((obj.e_w'*x_v).^2 + (obj.e_w'*y_v).^2);
            s_e = sqrt((obj.e_e'*x_v).^2 + (obj.e_e'*y_v).^2);
            s_s = sqrt((obj.e_s'*x_u).^2 + (obj.e_s'*y_u).^2);
            s_n = sqrt((obj.e_n'*x_u).^2 + (obj.e_n'*y_u).^2);

            obj.d_w = -1*(spdiag(1./s_w)*(a11_w*obj.du_w' + a12_w*obj.dv_w'))';
            obj.d_e =    (spdiag(1./s_e)*(a11_e*obj.du_e' + a12_e*obj.dv_e'))';
            obj.d_s = -1*(spdiag(1./s_s)*(a22_s*obj.dv_s' + a12_s*obj.du_s'))';
            obj.d_n =    (spdiag(1./s_n)*(a22_n*obj.dv_n' + a12_n*obj.du_n'))';

            obj.Dx = spdiag( y_v./J)*obj.Du + spdiag(-y_u./J)*obj.Dv;
            obj.Dy = spdiag(-x_v./J)*obj.Du + spdiag( x_u./J)*obj.Dv;

            %% Boundary inner products
            obj.H_w = obj.H_v*spdiag(s_w);
            obj.H_e = obj.H_v*spdiag(s_e);
            obj.H_s = obj.H_u*spdiag(s_s);
            obj.H_n = obj.H_u*spdiag(s_n);

            obj.a11 = a11;
            obj.a12 = a12;
            obj.a22 = a22;
            obj.s_w = spdiag(s_w);
            obj.s_e = spdiag(s_e);
            obj.s_s = spdiag(s_s);
            obj.s_n = spdiag(s_n);


            nu_w = [-1,0];
            nu_e = [1,0];
            nu_s = [0,-1];
            nu_n = [0,1];
            
            obj.n_w = cell(2,1);
            obj.n_e = cell(2,1);
            obj.n_s = cell(2,1);
            obj.n_n = cell(2,1);
            
            % Compute normal and rotate (exactly!) 90 degrees counter-clockwise to get tangent
            n_w_1 = (1./s_w).*obj.e_w'*(J.*(K{1,1}*nu_w(1) + K{1,2}*nu_w(2)));
            n_w_2 = (1./s_w).*obj.e_w'*(J.*(K{2,1}*nu_w(1) + K{2,2}*nu_w(2)));
            obj.n_w{1} = spdiag(n_w_1);
            obj.n1_w = spdiag(n_w_1);
            obj.n_w{2} = spdiag(n_w_2);
            obj.n2_w = spdiag(n_w_2);
            obj.tangent_w = {-obj.n_w{2}, obj.n_w{1}};
            obj.t1_w = -obj.n_w{2};
            obj.t2_w = obj.n_w{1};
            
            n_e_1 = (1./s_e).*obj.e_e'*(J.*(K{1,1}*nu_e(1) + K{1,2}*nu_e(2)));
            n_e_2 = (1./s_e).*obj.e_e'*(J.*(K{2,1}*nu_e(1) + K{2,2}*nu_e(2)));
            obj.n_e{1} = spdiag(n_e_1);
            obj.n1_e = spdiag(n_e_1);
            obj.n_e{2} = spdiag(n_e_2);
            obj.n2_e = spdiag(n_e_2);
            obj.tangent_e = {-obj.n_e{2}, obj.n_e{1}};
            obj.t1_e = -obj.n_e{2};
            obj.t2_e = obj.n_e{1};
            
            n_s_1 = (1./s_s).*obj.e_s'*(J.*(K{1,1}*nu_s(1) + K{1,2}*nu_s(2)));
            n_s_2 = (1./s_s).*obj.e_s'*(J.*(K{2,1}*nu_s(1) + K{2,2}*nu_s(2)));
            obj.n_s{1} = spdiag(n_s_1);
            obj.n1_s = spdiag(n_s_1);
            obj.n_s{2} = spdiag(n_s_2);
            obj.n2_s = spdiag(n_s_2);
            obj.tangent_s = {-obj.n_s{2}, obj.n_s{1}};
            obj.t1_s = -obj.n_s{2};
            obj.t2_s = obj.n_s{1};            
            
            n_n_1 = (1./s_n).*obj.e_n'*(J.*(K{1,1}*nu_n(1) + K{1,2}*nu_n(2)));
            n_n_2 = (1./s_n).*obj.e_n'*(J.*(K{2,1}*nu_n(1) + K{2,2}*nu_n(2)));
            obj.n_n{1} = spdiag(n_n_1);
            obj.n1_n = spdiag(n_n_1);
            obj.n_n{2} = spdiag(n_n_2);
            obj.n2_n = spdiag(n_n_2);
            obj.tangent_n = {-obj.n_n{2}, obj.n_n{1}};
            obj.t1_n = -obj.n_n{2};
            obj.t2_n = obj.n_n{1};

            % Temporary
            K{1,1} = spdiag(K{1,1});
            K{1,2} = spdiag(K{1,2});
            K{2,1} = spdiag(K{2,1});
            K{2,2} = spdiag(K{2,2});
            obj.K = K;
            obj.lambda = lambda;
        end


        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj, boundary, type, parameter)
            default_arg('type','neumann');
            default_arg('parameter', []);

            e               = obj.getBoundaryOperator('e', boundary);
            d               = obj.getBoundaryOperator('d', boundary);
            H_b             = obj.getBoundaryQuadrature(boundary);
            s_b             = obj.getBoundaryScaling(boundary);
            m               = obj.getBoundaryNumber(boundary);

            K = obj.K;
            J = obj.J;
            Hi = obj.Hi;

            switch type
                % Dirichlet boundary condition
                case {'D','d','dirichlet'}
                    tuning = 1.0;

                    sigma = 0*b;
                    for i = 1:obj.dim
                        sigma = sigma + b*J*K{i,m}*K{i,m};
                    end

                    % Minimum check on sigma
                    mx = obj.m(1);
                    my = obj.m(2);
                    bp = obj.bp;
                    sigma_mat = reshape(diag(sigma), my, mx);
                    switch boundary
                        case 'w'
                            sigma_min = min(sigma_mat(:,1:bp), [], 2);
                            sigma_mat = repmat(sigma_min, 1, mx);
                        case 'e'
                            sigma_min = min(sigma_mat(:,end-bp+1:end), [], 2);
                            sigma_mat = repmat(sigma_min, 1, mx);
                        case 's'
                            sigma_min = min(sigma_mat(1:bp,:), [], 1);
                            sigma_mat = repmat(sigma_min, my, 1);
                        case 'n'
                            sigma_min = min(sigma_mat(end-bp+1:end,:), [], 1);
                            sigma_mat = repmat(sigma_min, my, 1);
                    end
                    sigma_min = sigma_mat(:);
                    sigma_min = spdiag(sigma_min);

                    % Window
                    sigma = e'*sigma*e;
                    sigma_min = e'*sigma_min*e;

                    sigma = sigma/s_b;
                    sigma_min = sigma_min/s_b;

                    tau = tuning*(1/th_R*sigma/sigma_min*sigma + obj.dim/th_H*sigma);

                    closure = {Hi*d*H_b*e' ...
                        -Hi*e*tau*H_b*e',sparse(obj.grid.N(),obj.grid.N())};

                    penalty = -Hi*d*b_b*H_b ...
                        +Hi*e*tau*H_b;


                    % Neumann boundary condition. Note that the penalty is for du/dn and not b*du/dn.
                case {'N','n','neumann'}
                    tau1 = -1;
                    tau2 = 0;
                    tau = (tau1*e + tau2*d)*H_b;

                    closure =  {Hi*tau*d',sparse(obj.grid.N(),obj.grid.N())};
                    penalty = -Hi*tau;

                case {'O','o','outflow'}

                    tau = -e*H_b;

                    closure =  {obj.Hi*tau*d',obj.Hi*tau*e'};
                    penalty = -obj.Hi*tau;



                    %                     tau1 = -1;
                    %                     tau2 = 0;
                    %                     tau = (tau1*e + tau2*d)*H_b;

                    %                     closure =  {obj.a*obj.Hi*tau*d',sparse(obj.grid.N(),obj.grid.N())};
                    %                     penalty = -obj.a*obj.Hi*tau;
                    % Unknown, boundary condition
                otherwise
                    error('No such boundary condition: type = %s',type);
            end
        end

        % type     Struct that specifies the interface coupling.
        %          Fields:
        %          -- tuning:           penalty strength, defaults to 1.2
        %          -- interpolation:    type of interpolation, default 'none'
        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary,type,method)

            defaultType.tuning = 1.2;
            defaultType.interpolation = 'none';
            default_struct('type', defaultType);
            default_arg('method', 'sat');

            switch method
                case 'hybrid'
                    [closure, penalty] = interfaceHybrid(obj,boundary,neighbour_scheme,neighbour_boundary);
                case 'projection'
                    closure = 0;
                    penalty = 0;
                otherwise
                    error('Unknown method: %s ', method)
            end
        end

        function [closure, penalty] = interfaceHybrid(obj,boundary,neighbour_scheme,neighbour_boundary)
            e_u = obj.getBoundaryOperator('e', boundary);
            d_u = obj.getBoundaryOperator('d', boundary);
            H_b_u = obj.getBoundaryQuadrature(boundary);

            e_v = neighbour_scheme.getBoundaryOperator('e', neighbour_boundary);
            d_v = neighbour_scheme.getBoundaryOperator('d', neighbour_boundary);

            % check orientation of grid points
            my_bound_points = e_u'*obj.grid.points;
            neigh_bound_points = e_v'*neighbour_scheme.grid.points;
            if norm(my_bound_points - neigh_bound_points) > 1e-10
                d_v = flip(d_v,2);
            else
                d_v = d_v;
            end

            tau_u = -e_u*H_b_u;
            closure = 0.5*obj.Hi*tau_u*d_u';
            penalty = 0.5*obj.Hi*tau_u*d_v';

        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op        -- string
        % boundary  -- string
        function o = getBoundaryOperator(obj, op, boundary)
            if obj.debugging
                assertIsMember(op, {'e', 'd'})
                assertIsMember(boundary, {'w', 'e', 's', 'n'})
            end
            o = obj.([op, '_', boundary]);
        end

        function o = getBoundaryNormal(obj, op, boundary)
            assertIsMember(op, {'n1', 'n2', 't1', 't2'})
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            o = obj.([op, '_', boundary]);
        end

        % Returns square boundary quadrature matrix, of dimension
        % corresponding to the number of boundary points
        %
        % boundary -- string
        function H_b = getBoundaryQuadrature(obj, boundary)
            if obj.debugging
                assertIsMember(boundary, {'w', 'e', 's', 'n'})
            end

            H_b = obj.(['H_', boundary]);
        end

        % Returns square boundary quadrature scaling matrix, of dimension
        % corresponding to the number of boundary points
        %
        % boundary -- string
        function s_b = getBoundaryScaling(obj, boundary)
            if obj.debugging
                assertIsMember(boundary, {'w', 'e', 's', 'n'})
            end

            s_b = obj.(['s_', boundary]);
        end

        % Returns the coordinate number corresponding to the boundary
        %
        % boundary -- string
        function m = getBoundaryNumber(obj, boundary)
            if obj.debugging
                assertIsMember(boundary, {'w', 'e', 's', 'n'})
            end
            switch boundary
                case {'w', 'e'}
                    m = 1;
                case {'s', 'n'}
                    m = 2;
            end
        end

        % Returns the indices of the boundary points in the grid matrix
        % boundary -- string
        function I = getBoundaryIndices(obj, boundary)
            if obj.debugging
                assertIsMember(boundary, {'w', 'e', 's', 'n'})
            end

            ind = grid.funcToMatrix(obj.grid, 1:prod(obj.m));
            switch boundary
                case 'w'
                    I = ind(1,:);
                case 'e'
                    I = ind(end,:);
                case 's'
                    I = ind(:,1)';
                case 'n'
                    I = ind(:,end)';
            end
        end

        function N = size(obj)
            N = prod(obj.m);
        end
    end
end
