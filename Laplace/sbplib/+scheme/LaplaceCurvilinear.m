classdef LaplaceCurvilinear < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid

        order % Order accuracy for the approximation

        a,b % Parameters of the operator


        % Inner products and operators for physical coordinates
        D % Laplace operator
        H, Hi % Inner product
        e_w, e_e, e_s, e_n
        d_w, d_e, d_s, d_n % Normal derivatives at the boundary
        H_w, H_e, H_s, H_n % Boundary inner products
        Dx, Dy % Physical derivatives
        M % Gradient inner product

        % Metric coefficients
        J, Ji
        a11, a12, a22
        x_u
        x_v
        y_u
        y_v

        % Inner product and operators for logical coordinates
        H_u, H_v % Norms in the x and y directions
        Hi_u, Hi_v
        Hu,Hv % Kroneckerd norms. 1'*Hx*v corresponds to integration in the x dir.
        Hiu, Hiv
        du_w, dv_w
        du_e, dv_e
        du_s, dv_s
        du_n, dv_n
        gamm_u, gamm_v
        lambda

    end

    methods
        % Implements  a*div(b*grad(u)) as a SBP scheme
        % TODO: Implement proper H, it should be the real physical quadrature, the logic quadrature may be but in a separate variable (H_logic?)

        function obj = LaplaceCurvilinear(g ,order, a, b, opSet)
            default_arg('opSet',@sbp.D2Variable);
            default_arg('a', 1);
            default_arg('b', 1);

            if b ~=1
                error('Not implemented yet')
            end

            % assert(isa(g, 'grid.Curvilinear'))
            if isa(a, 'function_handle')
                a = grid.evalOn(g, a);
                a = spdiag(a);
            end

            m = g.size();
            m_u = m(1);
            m_v = m(2);
            m_tot = g.N();

            h = g.scaling();
            h_u = h(1);
            h_v = h(2);


            % 1D operators
            ops_u = opSet(m_u, {0, 1}, order);
            ops_v = opSet(m_v, {0, 1}, order);

            I_u = speye(m_u);
            I_v = speye(m_v);

            D1_u = ops_u.D1;
            D2_u = ops_u.D2;
            H_u =  ops_u.H;
            Hi_u = ops_u.HI;
            e_l_u = ops_u.e_l;
            e_r_u = ops_u.e_r;
            d1_l_u = ops_u.d1_l;
            d1_r_u = ops_u.d1_r;

            D1_v = ops_v.D1;
            D2_v = ops_v.D2;
            H_v =  ops_v.H;
            Hi_v = ops_v.HI;
            e_l_v = ops_v.e_l;
            e_r_v = ops_v.e_r;
            d1_l_v = ops_v.d1_l;
            d1_r_v = ops_v.d1_r;


            % Logical operators
            Du = kr(D1_u,I_v);
            Dv = kr(I_u,D1_v);
            obj.Hu  = kr(H_u,I_v);
            obj.Hv  = kr(I_u,H_v);
            obj.Hiu = kr(Hi_u,I_v);
            obj.Hiv = kr(I_u,Hi_v);
            obj.H_u = H_u;
            obj.H_v = H_v;

            e_w  = kr(e_l_u,I_v);
            e_e  = kr(e_r_u,I_v);
            e_s  = kr(I_u,e_l_v);
            e_n  = kr(I_u,e_r_v);
            obj.du_w = kr(d1_l_u,I_v);
            obj.dv_w = (e_w'*Dv)';
            obj.du_e = kr(d1_r_u,I_v);
            obj.dv_e = (e_e'*Dv)';
            obj.du_s = (e_s'*Du)';
            obj.dv_s = kr(I_u,d1_l_v);
            obj.du_n = (e_n'*Du)';
            obj.dv_n = kr(I_u,d1_r_v);


            % Metric coefficients
            coords = g.points();
            x = coords(:,1);
            y = coords(:,2);

            x_u = Du*x;
            x_v = Dv*x;
            y_u = Du*y;
            y_v = Dv*y;

            J = x_u.*y_v - x_v.*y_u;
            a11 =  1./J .* (x_v.^2  + y_v.^2);
            a12 = -1./J .* (x_u.*x_v + y_u.*y_v);
            a22 =  1./J .* (x_u.^2  + y_u.^2);
            lambda = 1/2 * (a11 + a22 - sqrt((a11-a22).^2 + 4*a12.^2));

            obj.x_u = x_u;
            obj.x_v = x_v;
            obj.y_u = y_u;
            obj.y_v = y_v;


            % Assemble full operators
            L_12 = spdiag(a12);
            Duv = Du*L_12*Dv;
            Dvu = Dv*L_12*Du;

            Duu = sparse(m_tot);
            Dvv = sparse(m_tot);
            ind = grid.funcToMatrix(g, 1:m_tot);

            for i = 1:m_v
                D = D2_u(a11(ind(:,i)));
                p = ind(:,i);
                Duu(p,p) = D;
            end

            for i = 1:m_u
                D = D2_v(a22(ind(i,:)));
                p = ind(i,:);
                Dvv(p,p) = D;
            end


            % Physical operators
            obj.J = spdiag(J);
            obj.Ji = spdiag(1./J);

            obj.D = obj.Ji*a*(Duu + Duv + Dvu + Dvv);
            obj.H = obj.J*kr(H_u,H_v);
            obj.Hi = obj.Ji*kr(Hi_u,Hi_v);

            obj.e_w = e_w;
            obj.e_e = e_e;
            obj.e_s = e_s;
            obj.e_n = e_n;

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

            s_w = sqrt((e_w'*x_v).^2 + (e_w'*y_v).^2);
            s_e = sqrt((e_e'*x_v).^2 + (e_e'*y_v).^2);
            s_s = sqrt((e_s'*x_u).^2 + (e_s'*y_u).^2);
            s_n = sqrt((e_n'*x_u).^2 + (e_n'*y_u).^2);

            obj.d_w = -1*(spdiag(1./s_w)*(a11_w*obj.du_w' + a12_w*obj.dv_w'))';
            obj.d_e =    (spdiag(1./s_e)*(a11_e*obj.du_e' + a12_e*obj.dv_e'))';
            obj.d_s = -1*(spdiag(1./s_s)*(a22_s*obj.dv_s' + a12_s*obj.du_s'))';
            obj.d_n =    (spdiag(1./s_n)*(a22_n*obj.dv_n' + a12_n*obj.du_n'))';

            obj.Dx = spdiag( y_v./J)*Du + spdiag(-y_u./J)*Dv;
            obj.Dy = spdiag(-x_v./J)*Du + spdiag( x_u./J)*Dv;

            %% Boundary inner products
            obj.H_w = H_v*spdiag(s_w);
            obj.H_e = H_v*spdiag(s_e);
            obj.H_s = H_u*spdiag(s_s);
            obj.H_n = H_u*spdiag(s_n);

            % Misc.
            obj.m = m;
            obj.h = [h_u h_v];
            obj.order = order;
            obj.grid = g;

            obj.a = a;
            obj.b = b;
            obj.a11 = a11;
            obj.a12 = a12;
            obj.a22 = a22;
            obj.lambda = lambda;

            obj.gamm_u = h_u*ops_u.borrowing.M.d1;
            obj.gamm_v = h_v*ops_v.borrowing.M.d1;
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

            e = obj.getBoundaryOperator('e', boundary);
            d = obj.getBoundaryOperator('d', boundary);
            H_b = obj.getBoundaryQuadrature(boundary);
            gamm = obj.getBoundaryBorrowing(boundary);
            switch type
                % Dirichlet boundary condition
                case {'D','d','dirichlet'}
                    tuning = 1.2;
                    % tuning = 20.2;

                    b1 = gamm*obj.lambda./obj.a11.^2;
                    b2 = gamm*obj.lambda./obj.a22.^2;

                    tau1 = tuning * spdiag(-1./b1 - 1./b2);
                    tau2 = 1;

                    tau = (tau1*e + tau2*d)*H_b;

                    closure =  obj.a*obj.Hi*tau*e';
                    penalty = -obj.a*obj.Hi*tau;


                % Neumann boundary condition
                case {'N','n','neumann'}
                    tau1 = -1;
                    tau2 = 0;
                    tau = (tau1*e + tau2*d)*H_b;

                    closure =  obj.a*obj.Hi*tau*d';
                    penalty = -obj.a*obj.Hi*tau;


                % Unknown, boundary condition
                otherwise
                    error('No such boundary condition: type = %s',type);
            end
        end

        % type     Struct that specifies the interface coupling.
        %          Fields:
        %          -- tuning:           penalty strength, defaults to 1.2
        %          -- interpolation:    type of interpolation, default 'none'
        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary,type)

            defaultType.tuning = 1.2;
            defaultType.interpolation = 'none';
            default_struct('type', defaultType);

            switch type.interpolation
            case {'none', ''}
                [closure, penalty] = interfaceStandard(obj,boundary,neighbour_scheme,neighbour_boundary,type);
            case {'op','OP','glue'}
                [closure, penalty] = interfaceNonConforming(obj,boundary,neighbour_scheme,neighbour_boundary,type);
            otherwise
                error('Unknown type of interpolation: %s ', type.interpolation);
            end
        end

        function [closure, penalty] = interfaceStandard(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            tuning = type.tuning;

            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            e_u    = obj.getBoundaryOperator('e', boundary);
            d_u    = obj.getBoundaryOperator('d', boundary);
            H_b_u = obj.getBoundaryQuadrature(boundary);
            I_u = obj.getBoundaryIndices(boundary);
            gamm_u = obj.getBoundaryBorrowing(boundary);

            e_v    = neighbour_scheme.getBoundaryOperator('e', neighbour_boundary);
            d_v    = neighbour_scheme.getBoundaryOperator('d', neighbour_boundary);
            H_b_v = neighbour_scheme.getBoundaryQuadrature(neighbour_boundary);
            I_v = neighbour_scheme.getBoundaryIndices(neighbour_boundary);
            gamm_v = neighbour_scheme.getBoundaryBorrowing(neighbour_boundary);

            u = obj;
            v = neighbour_scheme;

            b1_u = gamm_u*u.lambda(I_u)./u.a11(I_u).^2;
            b2_u = gamm_u*u.lambda(I_u)./u.a22(I_u).^2;
            b1_v = gamm_v*v.lambda(I_v)./v.a11(I_v).^2;
            b2_v = gamm_v*v.lambda(I_v)./v.a22(I_v).^2;

            tau1 = -1./(4*b1_u) -1./(4*b1_v) -1./(4*b2_u) -1./(4*b2_v);
            tau1 = tuning * spdiag(tau1);
            tau2 = 1/2;

            sig1 = -1/2;
            sig2 = 0;

            tau = (e_u*tau1 + tau2*d_u)*H_b_u;
            sig = (sig1*e_u + sig2*d_u)*H_b_u;

            closure = obj.a*obj.Hi*( tau*e_u' + sig*d_u');
            penalty = obj.a*obj.Hi*(-tau*e_v' + sig*d_v');
        end

        function [closure, penalty] = interfaceNonConforming(obj,boundary,neighbour_scheme,neighbour_boundary,type)

            % TODO: Make this work for curvilinear grids
            warning('LaplaceCurvilinear: Non-conforming grid interpolation only works for Cartesian grids.');

            % User can request special interpolation operators by specifying type.interpOpSet
            default_field(type, 'interpOpSet', @sbp.InterpOpsOP);
            interpOpSet = type.interpOpSet;
            tuning = type.tuning;


            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            e_u    = obj.getBoundaryOperator('e', boundary);
            d_u    = obj.getBoundaryOperator('d', boundary);
            H_b_u  = obj.getBoundaryQuadrature(boundary);
            I_u    = obj.getBoundaryIndices(boundary);
            gamm_u = obj.getBoundaryBorrowing(boundary);

            e_v    = neighbour_scheme.getBoundaryOperator('e', neighbour_boundary);
            d_v    = neighbour_scheme.getBoundaryOperator('d', neighbour_boundary);
            H_b_v  = neighbour_scheme.getBoundaryQuadrature(neighbour_boundary);
            I_v    = neighbour_scheme.getBoundaryIndices(neighbour_boundary);
            gamm_v = neighbour_scheme.getBoundaryBorrowing(neighbour_boundary);


            % Find the number of grid points along the interface
            m_u = size(e_u, 2);
            m_v = size(e_v, 2);

            Hi = obj.Hi;

            u = obj;
            v = neighbour_scheme;

            a_u = u.a;
            a_v = v.a;

            b1_u = gamm_u*u.lambda(I_u)./u.a11(I_u).^2;
            b2_u = gamm_u*u.lambda(I_u)./u.a22(I_u).^2;
            b1_v = gamm_v*v.lambda(I_v)./v.a11(I_v).^2;
            b2_v = gamm_v*v.lambda(I_v)./v.a22(I_v).^2;

            tau_u = -1./(4*b1_u) -1./(4*b2_u);
            tau_v = -1./(4*b1_v) -1./(4*b2_v);

            tau_u = tuning * spdiag(tau_u);
            tau_v = tuning * spdiag(tau_v);
            beta_u = tau_v;

            % Build interpolation operators
            switch type.interpolation
            case {'op','OP','MC','mc'}
                intOps = interpOpSet(m_u, m_v, obj.order, neighbour_scheme.order);
                Iu2v = intOps.Iu2v;
                Iv2u = intOps.Iv2u;
            case 'glue'
                optype = type.optype; % TODO: Should be stored by the scheme?
                xu = obj.getBoundaryGrid(boundary);
                xv = neighbour_scheme.getBoundaryGrid(neighbour_boundary);
                hu = obj.getBoundaryGridSpacing(boundary);
                hv = neighbour_scheme.getBoundaryGridSpacing(neighbour_boundary);
                glueOps = interpOpSet(xu, xv, hu, hv, obj.order, neighbour_scheme.order, optype, optype);
                Iu2v = glueOps.Iu2v;
                Iv2u = glueOps.Iv2u;
            otherwise
                error();
            end

            closure = a_u*Hi*e_u*tau_u*H_b_u*e_u' + ...
                      a_v*Hi*e_u*H_b_u*Iv2u.bad*beta_u*Iu2v.good*e_u' + ...
                      a_u*1/2*Hi*d_u*H_b_u*e_u' + ...
                      -a_u*1/2*Hi*e_u*H_b_u*d_u';

            penalty = -a_u*Hi*e_u*tau_u*H_b_u*Iv2u.good*e_v' + ...
                      -a_v*Hi*e_u*H_b_u*Iv2u.bad*beta_u*e_v' + ...
                      -a_u*1/2*Hi*d_u*H_b_u*Iv2u.good*e_v' + ...
                      -a_v*1/2*Hi*e_u*H_b_u*Iv2u.bad*d_v';

        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op        -- string
        % boundary  -- string
        function o = getBoundaryOperator(obj, op, boundary)
            assertIsMember(op, {'e', 'd'})
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            o = obj.([op, '_', boundary]);
        end

        % Returns square boundary quadrature matrix, of dimension
        % corresponding to the number of boundary points
        %
        % boundary -- string
        function H_b = getBoundaryQuadrature(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            H_b = obj.(['H_', boundary]);
        end

        % Returns the indices of the boundary points in the grid matrix
        % boundary -- string
        function I = getBoundaryIndices(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

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

        % Returns borrowing constant gamma
        % boundary -- string
        function gamm = getBoundaryBorrowing(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            switch boundary
                case {'w','e'}
                    gamm = obj.gamm_u;
                case {'s','n'}
                    gamm = obj.gamm_v;
            end
        end

        function xb = getBoundaryGrid(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})
            if isa(obj.grid,'grid.Cartesian')
                lg = obj.grid;
            elseif isa(obj.grid,'grid.Curvilinear')
                lg = obj.grid.logicalGrid();
            else
                error('Grid type not supported');
            end
            switch boundary
            case {'w','e'}
                xb = lg.x{2};
            case {'s','n'}
                xb = lg.x{1};
            end
        end

        function hb = getBoundaryGridSpacing(obj, boundary)
            h = obj.grid.scaling();
            switch boundary
            case {'w','e'}
                hb = h(2);
            case {'s','n'}
                hb = h(1);
            end
        end

        function N = size(obj)
            N = prod(obj.m);
        end
    end
end
