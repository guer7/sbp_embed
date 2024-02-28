% Analytical solution to the 2D acoustic wave equation with a point source
% f(t) = 1/(sigma*sqrt(2*pi))*exp(-(t - t_s).^2/(2*sigma^2));
function ur = uexact_func(r,t,t_s,sigma)
fint = @(omega,r) exp(-(t - t_s - r*cosh(omega)).^2/(2*sigma^2));

ur = zeros(numel(r),1);
for idx = 1:numel(r)
    ur(idx) = integral(@(omega) fint(omega,r(idx)),0,100, 'AbsTol', 1e-14, 'RelTol', 1e-14);
end
ur = 1/(sigma*(2*pi)^(3/2))*ur;