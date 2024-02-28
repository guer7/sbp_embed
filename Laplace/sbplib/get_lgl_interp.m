function [g_fun,g_fun_deriv] = get_lgl_interp(xydata,order)
xdata = xydata(:,1);
ydata = xydata(:,2);

j = (0:order)';
t = -cos(j*pi/order);

P = cgl_vand(order,t);
% dP=cgl_vand_der(order,x,P);

weights_x = P\xdata;
weights_y = P\ydata;

% teval = linspace(-1,1,101)';

% Peval = cgl_vand(order,teval);
% yeval =

    function x = xfun(t)
        Peval = cgl_vand(order,t);
        x = Peval*weights_x;
    end

    function dx = xprimefun(t)
        Peval = cgl_vand(order,t);
        dPeval = cgl_vand_der(order,t,Peval);
        dx = dPeval*weights_x;
    end

    function y = yfun(t)
        Peval = cgl_vand(order,t);
        y = Peval*weights_y;
    end

    function dy = yprimefun(t)
        Peval = cgl_vand(order,t);
        dPeval = cgl_vand_der(order,t,Peval);
        dy = dPeval*weights_y;
    end

    function v = g_fun_(t)
        v = [xfun(t)';yfun(t)'];
    end

    function v = g_fun_deriv_(t)
        v = [xprimefun(t)';yprimefun(t)'];
    end

    g_fun = @(t) g_fun_(2*t-1);
    g_fun_deriv = @(t) g_fun_deriv_(2*t-1);


end