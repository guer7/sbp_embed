function [a,b,c,s,varargout] = rk_parameters(order)

switch order
    case 3
        s = 3;
        a = sparse(s,s);
        a(2,1) = 1/2;
        a(3,1) = -1;
        a(3,2) = 2;
        b = 1/6*[1; 4; 1];
        c = [0; 1/2; 1];
    case 4
        s = 4;
        a = zeros(s,s);
        a(2,1) = 1/2;
        a(3,2) = 1/2;
        a(4,3) = 1;
        b = 1/6*[1; 2; 2; 1];
        c = [0; 1/2; 1/2; 1];
    case 5
        s = 6;
        a = zeros(s,s);
        a(2,1) = 1/4;
        a(3,1) = 3/32;
        a(3,2) = 9/32;
        a(4,1) = 1932/2197;
        a(4,2) = -7200/2197;
        a(4,3) = 7296/2197;
        a(5,1) = 439/216;
        a(5,2) = -8;
        a(5,3) = 3680/513;
        a(5,4) = -845/4104;
        a(6,1) = -8/27;
        a(6,2) = 2;
        a(6,3) = -3544/2565;
        a(6,4) = 1859/4104;
        a(6,5) = -11/40;
        b = [16/135; 0; 6656/12825; 28561/56430; -9/50; 2/55];
        c = [0; 1/4; 3/8; 12/13; 1; 1/2];
    case 6
        s = 7;
        a = zeros(s,s);
        a(2,1) = 4/7;
        a(3,1) = 115/112; a(3,2) = -5/16;
        a(4,1) = 589/630; a(4,2) = 5/18; a(4,3) = -16/45;
        a(5,1) = 229/1200 - 29/6000*sqrt(5); a(5,2) = 119/240 - 187/1200*sqrt(5); a(5,3) = -14/75 + 34/375*sqrt(5); a(5,4) = -3/100*sqrt(5);
        a(6,1) = 71/2400 - 587/12000*sqrt(5); a(6,2) = 187/480 - 391/2400*sqrt(5); a(6,3) = -38/75 + 26/375*sqrt(5); a(6,4) = 27/80 - 3/400*sqrt(5); a(6,5) = (1+sqrt(5))/4;
        a(7,1) = -49/480 + 43/160*sqrt(5); a(7,2) = -425/96 + 51/32*sqrt(5); a(7,3) = 52/15 - 4/5*sqrt(5); a(7,4) = -27/16 + 3/16*sqrt(5); a(7,5) = 5/4 - 3/4*sqrt(5); a(7,6) = 5/2 - 1/2*sqrt(5);
        b = [1/12; 0; 0; 0; 5/12; 5/12; 1/12];
        c = [0; 4/7; 5/7; 6/7; (5-sqrt(5))/10; (5+sqrt(5))/10; 1];
    case 32
        s = 3;
        a = sparse(s,s);
        a(2,1) = 1/2;
        a(3,1) = -1;
        a(3,2) = 2;
        b = 1/6*[1; 4; 1];
        bStar = 1/2*[1; 0; 1];
        c = [0; 1/2; 1];

        varargout{1} = bStar;
    case 54
        s = 6;
        a = zeros(s,s);
        a(2,1) = 1/4;
        a(3,1) = 3/32;
        a(3,2) = 9/32;
        a(4,1) = 1932/2197;
        a(4,2) = -7200/2197;
        a(4,3) = 7296/2197;
        a(5,1) = 439/216;
        a(5,2) = -8;
        a(5,3) = 3680/513;
        a(5,4) = -845/4104;
        a(6,1) = -8/27;
        a(6,2) = 2;
        a(6,3) = -3544/2565;
        a(6,4) = 1859/4104;
        a(6,5) = -11/40;
        b = [16/135; 0; 6656/12825; 28561/56430; -9/50; 2/55];
        bStar = [25/216; 0; 1408/2565; 2197/4104; -1/5; 0];
        c = [0; 1/4; 3/8; 12/13; 1; 1/2];

        varargout{1} = bStar;
    otherwise
        error("Order not implemented.")
end