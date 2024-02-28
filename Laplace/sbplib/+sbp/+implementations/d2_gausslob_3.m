% Central SBP on 4 gridpoints GL nodes (3rd in all points)
function [D1,H,x,e_l,e_r] = d2_gausslob_3(lim)
% Create operators for [-1,1] 
m=4;
d1 = 1 - sqrt(5)/5;
x=[-1 (d1-1) 1-d1 1]';

H=[0.1e1 / 0.6e1 0 0 0; 0 0.5e1 / 0.6e1 0 0; 0 0 0.5e1 / 0.6e1 0; 0 0 0 0.1e1 / 0.6e1;];
HI=inv(H);
Q=[0 0.5e1 / 0.24e2 * sqrt(0.5e1) + 0.5e1 / 0.24e2 -0.5e1 / 0.24e2 * sqrt(0.5e1) + 0.5e1 / 0.24e2 0.1e1 / 0.12e2; -0.5e1 / 0.24e2 * sqrt(0.5e1) - 0.5e1 / 0.24e2 0 0.5e1 / 0.12e2 * sqrt(0.5e1) -0.5e1 / 0.24e2 * sqrt(0.5e1) + 0.5e1 / 0.24e2; 0.5e1 / 0.24e2 * sqrt(0.5e1) - 0.5e1 / 0.24e2 -0.5e1 / 0.12e2 * sqrt(0.5e1) 0 0.5e1 / 0.24e2 * sqrt(0.5e1) + 0.5e1 / 0.24e2; -0.1e1 / 0.12e2 0.5e1 / 0.24e2 * sqrt(0.5e1) - 0.5e1 / 0.24e2 -0.5e1 / 0.24e2 * sqrt(0.5e1) - 0.5e1 / 0.24e2 0;];

e_l=zeros(m,1);e_l(1)=1;
e_r=zeros(m,1);e_r(m)=1;

D1=HI*(Q-1/2*e_l*e_l'+1/2*e_r*e_r') ;

% Rescale to lim 
L = lim{2} - lim{1};

x = 0.5*x + 0.5; % now between 0 and 1
x = lim{1} + L*x; % now between lim(1) and lim(2)

H = H*L/2;
D1 = D1*2/L;
