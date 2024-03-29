function [H, HI, D1, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = d4_variable_hollow_6(m,h)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 6:te ordn. SBP Finita differens         %%%
    %%% operatorer med diagonal norm            %%%
    %%% Extension to variable koeff             %%%
    %%%                                         %%%
    %%% H           (Normen)                    %%%
    %%% D1=H^(-1)Q  (approx f?rsta derivatan)   %%%
    %%% D2          (approx andra derivatan)    %%%
    %%% D2=HI*(R+C*D*S                          %%%
    %%%                                         %%%
    %%% R=-D1'*H*C*D1-RR                        %%%
    %%%                                         %%%
    %%% RR ?r dissipation)                      %%%
    %%% Dissipationen uppbyggd av D4:           %%%
    %%% DI=D4*B*H*D4                            %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    BP = 9;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    % Norm
    Hv = ones(m,1);
    Hv(1:6) = [13649/43200,12013/8640,2711/4320,5359/4320,7877/8640, 43801/43200];
    Hv(m-5:m) = rot90(Hv(1:6),2);
    Hv = h*Hv;
    H = spdiag(Hv, 0);
    HI = spdiag(1./Hv, 0);


    % Boundary operators
    e_l = sparse(m,1);
    e_l(1) = 1;
    e_r = rot90(e_l, 2);

    d1_l = sparse(m,1);
    d1_l(1:5) = [-25/12, 4, -3, 4/3, -1/4]/h;
    d1_r = -rot90(d1_l, 2);

    d2_l = sparse(m,1);
    d2_l(1:5) = [0.35e2/0.12e2 -0.26e2/0.3e1 0.19e2/0.2e1 -0.14e2/0.3e1 0.11e2/0.12e2;]/h^2;
    d2_r = rot90(d2_l, 2);

    d3_l = sparse(m,1);
    d3_l(1:5) = [-5/2 9 -12 7 -3/2]/h^3;
    d3_r = -rot90(d3_l, 2);


    % First derivtive
    x1=0.70127127127127;


    D1=(1/60*diag(ones(m-3,1),3)-9/60*diag(ones(m-2,1),2)+45/60*diag(ones(m-1,1),1)-45/60*diag(ones(m-1,1),-1)+9/60*diag(ones(m-2,1),-2)-1/60*diag(ones(m-3,1),-3));



    D1(1:6,1:9)=[-21600/13649, 43200/13649*x1-7624/40947, -172800/13649*x1+ ...
    	     715489/81894, 259200/13649*x1-187917/13649, -172800/13649* ...
    	     x1+735635/81894, 43200/13649*x1-89387/40947, 0, 0, 0; ...
    	     -8640/12013*x1+7624/180195, 0, 86400/12013*x1-57139/12013, ...
    	     -172800/12013*x1+745733/72078, 129600/12013*x1-91715/12013, ...
    	     -34560/12013*x1+240569/120130, 0, 0, 0; ...
             17280/2711*x1-715489/162660, -43200/2711*x1+57139/5422, 0, ...
             86400/2711*x1-176839/8133, -86400/2711*x1+242111/10844, ...
             25920/2711*x1-182261/27110, 0, 0, 0; ...
             -25920/5359*x1+187917/53590, 86400/5359*x1-745733/64308, ...
             -86400/5359*x1+176839/16077, 0, 43200/5359*x1-165041/32154, ...
             -17280/5359*x1+710473/321540, 72/5359, 0, 0; ...
             34560/7877*x1-147127/47262, -129600/7877*x1+91715/7877, ...
             172800/7877*x1-242111/15754, -86400/7877*x1+165041/23631, ...
             0, 8640/7877*x1, -1296/7877, 144/7877, 0; ...
             -43200/43801*x1+89387/131403, 172800/43801*x1-240569/87602,...
             -259200/43801*x1+182261/43801, 172800/43801*x1-710473/262806, ...
             -43200/43801*x1, 0, 32400/43801, -6480/43801, 720/43801];
    D1(m-5:m,m-8:m)=rot90( -D1(1:6,1:9),2);
    D1=D1/h;


    % Second derivative
    nBP = 9;
    M = sparse(m,m);
    coeffs = load('sbplib/+sbp/+implementations/coeffs_d2_variable_6.mat');

    function D2 = D2_fun(c)
        M_l = zeros(nBP, coeffs.nBPC);
        M_r = zeros(nBP, coeffs.nBPC);

        for i=1:coeffs.nBPC
            M_l = M_l + coeffs.C_l{i}*c(i);
            M_r = M_r + coeffs.C_r{i}*c(m-coeffs.nBPC+i);
        end

        M(1:nBP, 1:coeffs.nBPC) = M_l;
        M(m-nBP+1:m, m-coeffs.nBPC+1:m) = M_r;

        D2 = M/h^2;
    end
    D2 = @D2_fun;

    % Fourth derivative, 1th order accurate at first 8 boundary points (still
    % yield 5th order convergence if stable: for example u_tt=-u_xxxx
    stencil = [7/240, -2/5, 169/60, -122/15, 91/8, -122/15, 169/60, -2/5, 7/240];
    diags = -4:4;
    M4 = stripeMatrix(stencil, diags, m);

    M4_U = [
        0.1394226315049e13/0.367201486080e12 -0.1137054563243e13/0.114750464400e12 0.16614189027367e14/0.1836007430400e13 -0.1104821700277e13/0.306001238400e12 0.1355771086763e13/0.1836007430400e13 -0.27818686453e11/0.459001857600e12 -0.40671054239e11/0.1836007430400e13 0.5442887371e10/0.306001238400e12;
        -0.1137054563243e13/0.114750464400e12 0.70616795535409e14/0.2570410402560e13 -0.173266854731041e15/0.6426026006400e13 0.28938615291031e14/0.2570410402560e13 -0.146167361863e12/0.71400288960e11 0.2793470836571e13/0.12852052012800e14 0.6219558097e10/0.428401733760e12 -0.7313844559e10/0.166909766400e12;
        0.16614189027367e14/0.1836007430400e13 -0.173266854731041e15/0.6426026006400e13 0.378613061504779e15/0.12852052012800e14 -0.9117069604217e13/0.642602600640e12 0.632177582849e12/0.233673672960e12 -0.1057776382577e13/0.6426026006400e13 0.443019868399e12/0.4284017337600e13 -0.3707981e7/0.2318191200e10;
        -0.1104821700277e13/0.306001238400e12 0.28938615291031e14/0.2570410402560e13 -0.9117069604217e13/0.642602600640e12 0.5029150721885e13/0.514082080512e12 -0.5209119714341e13/0.1285205201280e13 0.12235427457469e14/0.12852052012800e14 -0.13731270505e11/0.64260260064e11 0.2933596129e10/0.40800165120e11;
        0.1355771086763e13/0.1836007430400e13 -0.146167361863e12/0.71400288960e11 0.632177582849e12/0.233673672960e12 -0.5209119714341e13/0.1285205201280e13 0.14871726798559e14/0.2570410402560e13 -0.7504337615347e13/0.1606506501600e13 0.310830296467e12/0.171360693504e12 -0.55284274391e11/0.183600743040e12;
        -0.27818686453e11/0.459001857600e12 0.2793470836571e13/0.12852052012800e14 -0.1057776382577e13/0.6426026006400e13 0.12235427457469e14/0.12852052012800e14 -0.7504337615347e13/0.1606506501600e13 0.106318657014853e15/0.12852052012800e14 -0.14432772918527e14/0.2142008668800e13 0.58102695589e11/0.22666758400e11;
        -0.40671054239e11/0.1836007430400e13 0.6219558097e10/0.428401733760e12 0.443019868399e12/0.4284017337600e13 -0.13731270505e11/0.64260260064e11 0.310830296467e12/0.171360693504e12 -0.14432772918527e14/0.2142008668800e13 0.27102479467823e14/0.2570410402560e13 -0.1216032192203e13/0.153000619200e12;
        0.5442887371e10/0.306001238400e12 -0.7313844559e10/0.166909766400e12 -0.3707981e7/0.2318191200e10 0.2933596129e10/0.40800165120e11 -0.55284274391e11/0.183600743040e12 0.58102695589e11/0.22666758400e11 -0.1216032192203e13/0.153000619200e12 0.20799922829107e14/0.1836007430400e13;
    ];

    M4(1:8,1:8) = M4_U;
    M4(m-7:m,m-7:m) = rot90(  M4_U ,2 );
    M4 = M4/h^3;



    D4=HI*(M4 - e_l*d3_l'+e_r*d3_r' + d1_l*d2_l'-d1_r*d2_r');
end
