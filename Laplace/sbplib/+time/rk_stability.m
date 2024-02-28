% Plot stability region of Runge-Kutta methods and print where the
% stability region intercepts the imagn
clear
close all

order = 4;

[a,b,~,s] = time.rk_parameters(order);

a = sym(a);
b = sym(b);

syms z

R = 1 + z*b'*inv(eye(s) - z*a)*ones(s,1);
R = simplify(R);

ruku = matlabFunction(abs(R));

figure
fcontour(ruku,[1,1],[-5, 2.6],[-5.2, 5.2],'k')
title(sprintf("Runge-Kutta stability region. Order: %d",order))
xlabel('Re')
ylabel('Im')
grid on
box on
axis([-5,3,-5,5])
hold off

rkroots(ruku)

function fcontour(f,levels,x_lim,y_lim,opt)
default_arg('opt','b')
x = linspace(x_lim(1),x_lim(2),10001);
y = linspace(y_lim(1),y_lim(2),10001);
[X,Y] = meshgrid(x,y);
mu = X+ 1i*Y;

z = f(mu);

contour(X,Y,z,levels,opt)

end


function rkroots(R)
% Roots for real evalues:
F = @(x)(real(R(x))-1);
real_x = fzero(F,-3);

% Roots for imaginary evalues:
F = @(x)(real(R(1i*x))-1);
imag_x1 = fzero(F,-3);
imag_x2 = fzero(F,3);

fprintf('Real x = %f\n',real_x)
fprintf('Imag x = %f\n',imag_x1)
fprintf('Imag x = %f\n',imag_x2)
end