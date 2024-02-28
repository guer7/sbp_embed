function [D1,H] = d1_noneq_8(N,h)

% N: Number of grid points
if(N<16)
    error('Operator requires at least 16 grid points');
end

% BP: Number of boundary points
BP = 8;

%%%% Norm matrix %%%%%%%%
P = sparse(BP,1);
%#ok<*NASGU>
P0 =  1.0758368078310e-01;
P1 =  6.1909685107891e-01;
P2 =  9.6971176519117e-01;
P3 =  1.1023441350947e+00;
P4 =  1.0244688965833e+00;
P5 =  9.9533550116831e-01;
P6 =  1.0008236941028e+00;
P7 =  9.9992060631812e-01;

for i = 0:BP-1
    P(i+1) = eval(['P' num2str(i)]);
end

H = ones(N,1);
H(1:BP) = P;
H(end-BP+1:end) = flip(P);
H = spdiags(h*H,0,N,N);
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Q matrix %%%%%%%%%%%
% interior stencil
order = 8;
d = [1/280,-4/105,1/5,-4/5,0,4/5,-1/5,4/105,-1/280];
d = repmat(d,N,1);
Q = spdiags(d,-order/2:order/2,N,N);

% Boundaries
Q0_0 = -5.0000000000000e-01;
Q0_1 =  6.7284756079369e-01;
Q0_2 = -2.5969732837062e-01;
Q0_3 =  1.3519390385721e-01;
Q0_4 = -6.9678474730984e-02;
Q0_5 =  2.6434024071371e-02;
Q0_6 = -5.5992311465618e-03;
Q0_7 =  4.9954552590464e-04;
Q0_8 =  0.0000000000000e+00;
Q0_9 =  0.0000000000000e+00;
Q0_10 =  0.0000000000000e+00;
Q0_11 =  0.0000000000000e+00;
Q1_0 = -6.7284756079369e-01;
Q1_1 =  0.0000000000000e+00;
Q1_2 =  9.4074021172233e-01;
Q1_3 = -4.0511642426516e-01;
Q1_4 =  1.9369192209331e-01;
Q1_5 = -6.8638079843479e-02;
Q1_6 =  1.3146457241484e-02;
Q1_7 = -9.7652615479254e-04;
Q1_8 =  0.0000000000000e+00;
Q1_9 =  0.0000000000000e+00;
Q1_10 =  0.0000000000000e+00;
Q1_11 =  0.0000000000000e+00;
Q2_0 =  2.5969732837062e-01;
Q2_1 = -9.4074021172233e-01;
Q2_2 =  0.0000000000000e+00;
Q2_3 =  9.4316393361096e-01;
Q2_4 = -3.5728039257451e-01;
Q2_5 =  1.1266686855013e-01;
Q2_6 = -1.8334941452280e-02;
Q2_7 =  8.2741521740941e-04;
Q2_8 =  0.0000000000000e+00;
Q2_9 =  0.0000000000000e+00;
Q2_10 =  0.0000000000000e+00;
Q2_11 =  0.0000000000000e+00;
Q3_0 = -1.3519390385721e-01;
Q3_1 =  4.0511642426516e-01;
Q3_2 = -9.4316393361096e-01;
Q3_3 =  0.0000000000000e+00;
Q3_4 =  8.7694387866575e-01;
Q3_5 = -2.4698058719506e-01;
Q3_6 =  4.7291642094198e-02;
Q3_7 = -4.0135203618880e-03;
Q3_8 =  0.0000000000000e+00;
Q3_9 =  0.0000000000000e+00;
Q3_10 =  0.0000000000000e+00;
Q3_11 =  0.0000000000000e+00;
Q4_0 =  6.9678474730984e-02;
Q4_1 = -1.9369192209331e-01;
Q4_2 =  3.5728039257451e-01;
Q4_3 = -8.7694387866575e-01;
Q4_4 =  0.0000000000000e+00;
Q4_5 =  8.1123946853807e-01;
Q4_6 = -2.0267150541446e-01;
Q4_7 =  3.8680398901392e-02;
Q4_8 = -3.5714285714286e-03;
Q4_9 =  0.0000000000000e+00;
Q4_10 =  0.0000000000000e+00;
Q4_11 =  0.0000000000000e+00;
Q5_0 = -2.6434024071371e-02;
Q5_1 =  6.8638079843479e-02;
Q5_2 = -1.1266686855013e-01;
Q5_3 =  2.4698058719506e-01;
Q5_4 = -8.1123946853807e-01;
Q5_5 =  0.0000000000000e+00;
Q5_6 =  8.0108544742793e-01;
Q5_7 = -2.0088756283071e-01;
Q5_8 =  3.8095238095238e-02;
Q5_9 = -3.5714285714286e-03;
Q5_10 =  0.0000000000000e+00;
Q5_11 =  0.0000000000000e+00;
Q6_0 =  5.5992311465618e-03;
Q6_1 = -1.3146457241484e-02;
Q6_2 =  1.8334941452280e-02;
Q6_3 = -4.7291642094198e-02;
Q6_4 =  2.0267150541446e-01;
Q6_5 = -8.0108544742793e-01;
Q6_6 =  0.0000000000000e+00;
Q6_7 =  8.0039405922650e-01;
Q6_8 = -2.0000000000000e-01;
Q6_9 =  3.8095238095238e-02;
Q6_10 = -3.5714285714286e-03;
Q6_11 =  0.0000000000000e+00;
Q7_0 = -4.9954552590464e-04;
Q7_1 =  9.7652615479254e-04;
Q7_2 = -8.2741521740941e-04;
Q7_3 =  4.0135203618880e-03;
Q7_4 = -3.8680398901392e-02;
Q7_5 =  2.0088756283071e-01;
Q7_6 = -8.0039405922650e-01;
Q7_7 =  0.0000000000000e+00;
Q7_8 =  8.0000000000000e-01;
Q7_9 = -2.0000000000000e-01;
Q7_10 =  3.8095238095238e-02;
Q7_11 = -3.5714285714286e-03;
for i = 1:BP
    for j = 1:BP
        Q(i,j) = eval(['Q' num2str(i-1) '_' num2str(j-1)]);
        Q(N+1-i,N+1-j) = -eval(['Q' num2str(i-1) '_' num2str(j-1)]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Difference operator %%
D1 = H\Q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%