clear
close all

m = 401;
lw = 4;
fs = 36;

x_l = -1;
x_r = 1;
x_i = 0;
y_l = 0;
y_r = 1;
amp = 0.05;
sine_curve = @(x) amp*sin(4*pi*x);

xi = linspace(0,1,m)';
eta = linspace(0,1,m)';

colors = [166,206,227;
 31,120,180;
 178,223,138;
 51,160,44]/255;

figure('pos',[595         401        1397         726])
hold on
box on

xfun = @(xi,eta) x_l + sine_curve(eta) + xi;
yfun = @(xi,eta) y_l + sine_curve(xi) + eta;
plot(xfun(0,eta),yfun(0,eta),'Color',colors(1,:),'Linewidth',lw) % west
plot(xfun(1,eta),yfun(1,eta),'Color',colors(3,:),'Linewidth',lw) % east
plot(xfun(xi,0),yfun(xi,0),'Color',colors(2,:),'Linewidth',lw) % south
plot(xfun(xi,1),yfun(xi,1),'Color',colors(4,:),'Linewidth',lw) % north

xfun = @(xi,eta) x_i + sine_curve(eta) + xi;
yfun = @(xi,eta) y_l + sine_curve(xi) + eta;
% plot(xfun(0,eta),yfun(0,eta),'Color',colors(3,:),'Linewidth',lw) % west
plot(xfun(1,eta),yfun(1,eta),'Color',colors(1,:),'Linewidth',lw) % east
plot(xfun(xi,0),yfun(xi,0),'Color',colors(2,:),'Linewidth',lw) % south
plot(xfun(xi,1),yfun(xi,1),'Color',colors(4,:),'Linewidth',lw) % north

xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
set(gca,'Fontsize',fs)
axis([x_l-amp-0.02,x_r+amp+0.02,y_l-amp-0.02,y_r+amp+0.02])

tb1 = annotation("textbox");
tb1.String = "$\partial \Omega_1$";
tb1.EdgeColor = 'none';
tb1.Interpreter = 'latex';
tb1.FontSize = fs;
tb1.Position = [0.3199    0.8080    0.1000    0.0400];

tb2 = annotation("textbox");
tb2.String = "$\partial \Omega_2$";
tb2.EdgeColor = 'none';
tb2.Interpreter = 'latex';
tb2.FontSize = fs;
tb2.Position = [0.6885    0.2480    0.1000    0.0400];
% 
tb3 = annotation("textbox");
tb3.String = "$\partial \Omega_3$";
tb3.EdgeColor = 'none';
tb3.Interpreter = 'latex';
tb3.FontSize = fs;
tb3.Position = [0.1643    0.5149    0.1000    0.0400];

tbI = annotation("textbox");
tbI.String = "$\partial \Omega_I$";
tbI.EdgeColor = 'none';
tbI.Interpreter = 'latex';
tbI.FontSize = fs;
tbI.Position = [0.53   0.5080    0.1000    0.0400];

tbl = annotation("textbox");
tbl.String = "$\Omega^{(l)}$";
tbl.EdgeColor = 'none';
tbl.Interpreter = 'latex';
tbl.FontSize = fs;
tbl.Position = [0.3199    0.5080    0.1000    0.0400];

tbr = annotation("textbox");
tbr.String = "$\Omega^{(r)}$";
tbr.EdgeColor = 'none';
tbr.Interpreter = 'latex';
tbr.FontSize = fs;
tbr.Position = [0.6885    0.5053    0.1000    0.0400];

exportgraphics(gcf,'multiblock_example_domain.pdf','ContentType','vector')


