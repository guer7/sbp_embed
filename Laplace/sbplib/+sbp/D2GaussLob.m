classdef D2GaussLob < sbp.OpSet
    properties
        D1 % SBP operator approximating first derivative
        H % Norm matrix
        HI % H^-1
        Q % Skew-symmetric matrix
        e_l % Left boundary operator
        e_r % Right boundary operator
        D2 % SBP operator for second derivative
        M % Norm matrix, second derivative
        d1_l % Left boundary first derivative
        d1_r % Right boundary first derivative
        m % Number of grid points.
        h % Step size
        x % grid
        borrowing % Struct with borrowing limits for different norm matrices
    end

    methods
        function obj = D2GaussLob(m,lim,order)

            x_l = lim{1};
            x_r = lim{2};
            L = x_r-x_l;
            obj.h = L/(m-1);
%             obj.x = linspace(x_l,x_r,m)';

            switch order

                case 3
                    assert(m == 4)
                    [obj.D1, obj.H, obj.x, obj.e_l, obj.e_r] = ...
                        sbp.implementations.d2_gausslob_3(lim);
                case 4
                    assert(m == 5)
                    [obj.D1, obj.H, obj.x, obj.e_l, obj.e_r] = ...
                        sbp.implementations.d2_gausslob_4(lim);
%                     2
                case 5
                    assert(m == 6)
                    [obj.D1, obj.H, obj.x, obj.e_l, obj.e_r] = ...
                        sbp.implementations.d2_gausslob_5(lim);
%                     2
                case 6
                    assert(m == 7)
                    [obj.D1, obj.H, obj.x, obj.e_l, obj.e_r] = ...
                        sbp.implementations.d2_gausslob_6(lim);
%                     2
                case 7
                    assert(m == 8)
                    [obj.D1, obj.H, obj.x, obj.e_l, obj.e_r] = ...
                        sbp.implementations.d2_gausslob_7(lim);
%                     2
                case 9
                    assert(m == 10)
                    [obj.D1, obj.H, obj.x, obj.e_l, obj.e_r] = ...
                        sbp.implementations.d2_gausslob_9(lim);
%                     2
                otherwise
                    error('Invalid operator order %d.',order);
            end
            obj.borrowing.H11 = obj.H(1,1)/obj.h; % First element in H/h,
            obj.m = m;
            obj.D2 = @(a) obj.D1*spdiag(a)*obj.D1;
            obj.d1_l = obj.D1'*obj.e_l;
            obj.d1_r = obj.D1'*obj.e_r;
            obj.HI = inv(obj.H);
            obj.M = [];
        end
        function str = string(obj)
            str = [class(obj) '_' num2str(obj.order)];
        end
    end


end





