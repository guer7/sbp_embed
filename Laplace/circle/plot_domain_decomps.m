% Plot block decomposition of unit circle.
clear
close all

addpath('../sbplib')
addpath('../utils/')

h = 0.1;
mesh_filename = "../mesh_files/circ_h" + num2str(h) + "/circle.mesh";

dom_circ = multiblock.domain.Circle();
dom_hohq_circ = multiblock.HOHQMesh(mesh_filename);

fs = 32;

figure('pos',[1233         438         943         899])
hold on
box on

for bidx = 1:dom_circ.nBlocks
    dom_circ.blockMaps{bidx}.plot(2,2);
end
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
set(gca,'Fontsize',fs)

exportgraphics(gcf,'circle_domain_decomp.pdf','ContentType','vector')

figure('pos',[1233         438         943         899])
hold on
box on

for bidx = 1:dom_hohq_circ.nBlocks
    dom_hohq_circ.blockMaps{bidx}.plot(2,2);
end
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
set(gca,'Fontsize',fs)

exportgraphics(gcf,'circle_hohq_domain_decomp.pdf','ContentType','vector')