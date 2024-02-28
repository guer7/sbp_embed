% Upwind SBP on 5 gridpoints
function [D1,H,x,e_l,e_r] = d2_gausslob_4(lim)
% Create operators for [-1,1] 
m=5;
d1 = 1 - sqrt(21)/7;
x=[-1 (d1-1) 0 1-d1 1]';

H=[0.1e1 / 0.10e2 0 0 0 0; 0 0.49e2 / 0.90e2 0 0 0; 0 0 0.32e2 / 0.45e2 0 0; 0 0 0 0.49e2 / 0.90e2 0; 0 0 0 0 0.1e1 / 0.10e2;];
HI=inv(H);
Q=[0 0.49e2 / 0.120e3 + 0.7e1 / 0.120e3 * sqrt(0.21e2) -0.4e1 / 0.15e2 0.49e2 / 0.120e3 - 0.7e1 / 0.120e3 * sqrt(0.21e2) -0.1e1 / 0.20e2; -0.49e2 / 0.120e3 - 0.7e1 / 0.120e3 * sqrt(0.21e2) 0 0.28e2 / 0.135e3 * sqrt(0.21e2) -0.49e2 / 0.540e3 * sqrt(0.21e2) 0.49e2 / 0.120e3 - 0.7e1 / 0.120e3 * sqrt(0.21e2); 0.4e1 / 0.15e2 -0.28e2 / 0.135e3 * sqrt(0.21e2) 0 0.28e2 / 0.135e3 * sqrt(0.21e2) -0.4e1 / 0.15e2; -0.49e2 / 0.120e3 + 0.7e1 / 0.120e3 * sqrt(0.21e2) 0.49e2 / 0.540e3 * sqrt(0.21e2) -0.28e2 / 0.135e3 * sqrt(0.21e2) 0 0.49e2 / 0.120e3 + 0.7e1 / 0.120e3 * sqrt(0.21e2); 0.1e1 / 0.20e2 -0.49e2 / 0.120e3 + 0.7e1 / 0.120e3 * sqrt(0.21e2) 0.4e1 / 0.15e2 -0.49e2 / 0.120e3 - 0.7e1 / 0.120e3 * sqrt(0.21e2) 0;];

e_l=zeros(m,1);e_l(1)=1;
e_r=zeros(m,1);e_r(m)=1;

D1=HI*(Q-1/2*e_l*e_l'+1/2*e_r*e_r') ;

% Rescale to lim 
L = lim{2} - lim{1};

x = 0.5*x + 0.5; % now between 0 and 1
x = lim{1} + L*x; % now between lim(1) and lim(2)

H = H*L/2;
D1 = D1*2/L;
