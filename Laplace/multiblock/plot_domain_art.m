% Plot boundary conditions and block decomposition of domain
clear
close all

fs = 32;
order_bot = 5;
order_water = 4;
opSet = {@sbp.D2GaussLob,order_bot,order_bot+1,'GL'};
animate = 1;

c_water = 1;
c_bot = 1;

im_folder = "seabottom";
addpath('../sbplib')
addpath('../utils/')
order_bot = opSet{2};

mesh_filename = "../mesh_files/seabottom_h4/seabottom.mesh";
m_water = 151;

Tend = 0.6;
prob = "pointsource";
uexact = @(r,t) uexact_func(r,t,t_s,sigma);
utexact = @(x,y,t) 0*x;

dom_water = multiblock.domain.Rectangle([0,100],[100,0]);
G_water_ = dom_water.getGrid([m_water,m_water]);
G_water = G_water_.grids{1};

lapl_water = scheme.LaplaceCurvilinearNoPars(G_water,order_water);
N_water = lapl_water.size;
tic

% Read mesh and create "plus" operators
dom_bot = multiblock.HOHQMesh(mesh_filename);

figure('units','pixels','pos',[310          48         714        1231])
hold on
box on
dom_water.blockTi{1}.plot(2,2);
for bidx = 1:dom_bot.nBlocks
    dom_bot.blockMaps{bidx}.plot(2,2);
end
axis([0,100,-100,100])
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
set(gca,'Fontsize',fs)

exportgraphics(gcf,'multiblock_domain_decomp.pdf','ContentType','vector')

colors = [166,206,227;
 31,120,180;
 178,223,138;
 51,160,44]/255;

N = 100;
figure('units','pixels','pos',[310          48         714        1231])
hold on
box on
for bidx = 1:length(dom_bot.boundaryGroups.outer)
    bound = dom_bot.boundaryGroups.outer{bidx};
    block_obj = dom_bot.blockMaps{bound{1}};
    switch bound{2}
        case 's'
            block_obj.gs{1}.plot(N,colors(1,:));
        case 'e'
            block_obj.gs{2}.plot(N,colors(1,:));
        case 'n'
            block_obj.gs{3}.plot(N,colors(1,:));
        case 'w'
            block_obj.gs{4}.plot(N,colors(1,:));
    end
end
for bidx = 1:length(dom_bot.boundaryGroups.inner)
    bound = dom_bot.boundaryGroups.inner{bidx};
    block_obj = dom_bot.blockMaps{bound{1}};
    switch bound{2}
        case 's'
            block_obj.gs{1}.plot(N,colors(2,:));
        case 'e'
            block_obj.gs{2}.plot(N,colors(2,:));
        case 'n'
            block_obj.gs{3}.plot(N,colors(2,:));
        case 'w'
            block_obj.gs{4}.plot(N,colors(2,:));
    end
end
for bidx = 1:length(dom_bot.boundaryGroups.interface)
    bound = dom_bot.boundaryGroups.interface{bidx};
    block_obj = dom_bot.blockMaps{bound{1}};
    switch bound{2}
        case 's'
            block_obj.gs{1}.plot(N,colors(3,:));
        case 'e'
            block_obj.gs{2}.plot(N,colors(3,:));
        case 'n'
            block_obj.gs{3}.plot(N,colors(3,:));
        case 'w'
            block_obj.gs{4}.plot(N,colors(3,:));
    end
end

for bidx = 1:length(dom_water.boundaryGroups.all)
    bound = dom_water.boundaryGroups.all{bidx};
    block_obj = dom_water.blockTi{bound{1}};
    switch bound{2}
        case 's'
%             block_obj.gs{1}.plot(N,colors(1,:));
        case 'e'
            block_obj.gs{2}.plot(N,colors(1,:));
        case 'n'
            block_obj.gs{3}.plot(N,colors(4,:));
        case 'w'
            block_obj.gs{4}.plot(N,colors(1,:));
    end
end
axis([0,100,-100,100])
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
set(gca,'Fontsize',fs)

tb1 = annotation("textbox");
tb1.String = "$\partial \Omega_1$";
tb1.EdgeColor = 'none';
tb1.Interpreter = 'latex';
tb1.FontSize = fs;
tb1.Position = [0.4860    0.8662    0.1000    0.0400];

tb2 = annotation("textbox");
tb2.String = "$\partial \Omega_2$";
tb2.EdgeColor = 'none';
tb2.Interpreter = 'latex';
tb2.FontSize = fs;
tb2.Position = [0.5137    0.2874    0.1000    0.0400];

tb3 = annotation("textbox");
tb3.String = "$\partial \Omega_3$";
tb3.EdgeColor = 'none';
tb3.Interpreter = 'latex';
tb3.FontSize = fs;
tb3.Position = [0.1786    0.6968    0.1000    0.0400];

tbI = annotation("textbox");
tbI.String = "$\partial \Omega_I$";
tbI.EdgeColor = 'none';
tbI.Interpreter = 'latex';
tbI.FontSize = fs;
tbI.Position = [0.4902    0.5257    0.1000    0.0400];

exportgraphics(gcf,'multiblock_domain_bc.pdf','ContentType','vector')