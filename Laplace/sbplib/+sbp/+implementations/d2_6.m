function [H, HI, D1, D2, e_1, e_m, M, Q, S_1, S_m] = d2_6(m,h)
    
    BP = 6;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

    H=speye(m,m);
    H(1:6,1:6)=diag([13649/43200,12013/8640,2711/4320,5359/4320,7877/8640, ...
             43801/43200]);
    H(m-5:m,m-5:m)=rot90(diag([13649/43200,12013/8640, ...
                2711/4320,5359/4320,7877/8640,43801/43200]),2);

    H=H*h;
    HI=inv(H);


    % D1 har en fri parameter x1.
    % Ett optimerat varde ger x1=0.70127127127127 = 331/472
    x1=0.70127127127127;

    diags   = -3:3;
    stencil = [-1/60,3/20,-3/4,0,3/4,-3/20,1/60];
    D1 = stripeMatrix(stencil, diags, m);

    D1(1:6,1:9)=[-21600/13649, 43200/13649*x1-7624/40947, -172800/13649*x1+ ...
             715489/81894, 259200/13649*x1-187917/13649, -172800/13649* ...
             x1+735635/81894, 43200/13649*x1-89387/40947, 0, 0, 0; ...
             -8640/12013*x1+7624/180195, 0, 86400/12013*x1-57139/12013, ...
             -172800/12013*x1+745733/72078, 129600/12013*x1-91715/12013, ...
             -34560/12013*x1+240569/120130, 0, 0, 0; 17280/2711*x1-715489/162660, -43200/2711*x1+57139/5422, 0, 86400/2711*x1-176839/8133, -86400/2711*x1+242111/10844, 25920/2711*x1-182261/27110, 0, 0, 0; -25920/5359*x1+187917/53590, 86400/5359*x1-745733/64308, -86400/5359*x1+176839/16077, 0, 43200/5359*x1-165041/32154, -17280/5359*x1+710473/321540, 72/5359, 0, 0; 34560/7877*x1-147127/47262, -129600/7877*x1+91715/7877, 172800/7877*x1-242111/15754, -86400/7877*x1+165041/23631, 0, 8640/7877*x1, -1296/7877, 144/7877, 0; -43200/43801*x1+89387/131403, 172800/43801*x1-240569/87602, -259200/43801*x1+182261/43801, 172800/43801*x1-710473/262806, -43200/43801*x1, 0, 32400/43801, -6480/43801, 720/43801];
    D1(m-5:m,m-8:m)=rot90( -D1(1:6,1:9),2);
    D1=D1/h;

    Q=H*D1 + 1/2*(e_1*e_1') - 1/2*(e_m*e_m');

    %D2=(2*diag(ones(m-3,1),3)-27*diag(ones(m-2,1),2)+270*diag(ones(m-1,1),1)+270*diag(ones(m-1,1),-1)-27*diag(ones(m-2,1),-2)+2*diag(ones(m-3,1),-3)-490*diag(ones(m,1),0))/180;
    diags   = -3:3;
    stencil = 1/180*[2,-27,270,-490,270,-27,2];
    D2 = stripeMatrix(stencil, diags, m);
    
    D2(1:6,1:9)=[114170/40947, -438107/54596, 336409/40947, -276997/81894, 3747/13649, 21035/163788, 0, 0, 0;6173/5860, -2066/879, 3283/1758, -303/293, 2111/3516, -601/4395, 0, 0, 0;-52391/81330, 134603/32532, -21982/2711, 112915/16266, -46969/16266, 30409/54220, 0, 0, 0;68603/321540, -12423/10718, 112915/32154, -75934/16077, 53369/21436, -54899/160770, 48/5359, 0, 0;-7053/39385, 86551/94524, -46969/23631, 53369/15754, -87904/23631, 820271/472620, -1296/7877, 96/7877, 0;21035/525612, -24641/131403, 30409/87602, -54899/131403, 820271/525612, -117600/43801, 64800/43801, -6480/43801, 480/43801];
    D2(m-5:m,m-8:m)=rot90( D2(1:6,1:9) ,2 );

    D2=D2/h^2;

    S_U=[-25/12, 4, -3, 4/3, -1/4]/h;
    S_1=sparse(1,m);
    S_1(1:5)=S_U;
    S_m=sparse(1,m);
    S_m(m-4:m)=fliplr(-S_U);


    M=-H*D2-e_1*S_1+e_m*S_m;
    S_1 = S_1';
    S_m = S_m';
end