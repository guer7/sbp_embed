% Plots and prints the numerical results from paper
% "Effcient multi-block discretization of the Laplacian on complex geometries"

clear
close all

do_comps = 1;
print_conv_table = 1;
plot_conv = 1;
plot_eff1 = 0;
plot_eff2 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Do simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_comps
    hvec_GL = [0.5,0.4,0.3,0.2,0.1,0.08,0.06];
    order_vec_GL = [5,7,9];
    op_GL = @sbp.D2GaussLob;

    for j = 1:numel(order_vec_GL)
        for i = 1:numel(hvec_GL)
            opSet = {op_GL,order_vec_GL(j),order_vec_GL(j)+1,'GL'};
            fprintf("--------- Embedding ----------\n")
            [errvec_GL(i,j),dofsvec_GL(i,j),specradvec_GL(i,j),time_elapsedvec_GL(i,j),nonzerosvec_GL(i,j),Nelementsvec_GL(i,j)] = wave_circle(order_vec_GL(j)+1,order_vec_GL(j),'embed',0,hvec_GL(i));
            fprintf("Order: %d, elements: %d, DOFS: %d, error: %e\n",order_vec_GL(j),Nelementsvec_GL(i,j),dofsvec_GL(i,j),errvec_GL(i,j))
        end
        [dofsvec_GL(:,j),I] = sort(dofsvec_GL(:,j));
        errvec_GL(:,j) = errvec_GL(I,j);
        specradvec_GL(:,j) = specradvec_GL(I,j);
        time_elapsedvec_GL(:,j) = time_elapsedvec_GL(I,j);
        nonzerosvec_GL(:,j) = nonzerosvec_GL(I,j);
    end

    mvec_trad = [41,61,81,101,121,161,201];
    order_vec_trad = [2,4,6];
    op_trad = @sbp.D2Variable;
    for j = 1:numel(order_vec_trad)
        for i = 1:numel(mvec_trad)
            opSet = {op_trad,order_vec_trad(j),mvec_trad(i),'equidist',order_vec_trad(j)};
            fprintf("--------- Traditional ----------\n")
            [errvec_trad(i,j),dofsvec_trad(i,j),specradvec_trad(i,j),time_elapsedvec_trad(i,j),nonzerosvec_trad(i,j)] = wave_circle(mvec_trad(i),order_vec_trad(j),'trad',0);
            fprintf("Order: %d, DOFS: %d, error: %e\n",order_vec_trad(j),dofsvec_trad(i,j),errvec_trad(i,j))
            if dofsvec_trad(i,j) > 1e5
                break
            end
        end
    end

    mvec_bopt = [61,81,101,121,161,201];
    order_vec_bopt = [8,10,12];
    op_bopt = @sbp.D2Nonequidistant;
    for j = 1:numel(order_vec_bopt)
        for i = 1:numel(mvec_bopt)
            opSet = {op_bopt,order_vec_bopt(j),mvec_bopt(i),'boundaryopt','accurate'};
            fprintf("--------- Boundary optimized ----------\n")
            [errvec_bopt(i,j),dofsvec_bopt(i,j),specradvec_bopt(i,j),time_elapsedvec_bopt(i,j),nonzerosvec_bopt(i,j)] = wave_circle(mvec_bopt(i),order_vec_bopt(j),'boundOpt',0);
            fprintf("Order: %d, DOFS: %d, error: %e\n",order_vec_bopt(j),dofsvec_bopt(i,j),errvec_bopt(i,j))
            if dofsvec_bopt(i,j) > 1e5
                break
            end
        end
    end

    save("meas_data.mat")
    %     return
else
    load("meas_data.mat")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = 5;
fs = 18;
ms = 20;

leg = {};
figure('pos',[1000         275        1176        1062])
hold on
box on
grid on
marker_types = {'*-','x--','d:','+-.'};
count = 1;
for j = 1:numel(order_vec_GL)
    plot(dofsvec_GL(:,j),errvec_GL(:,j),marker_types{count},'Color',[27,158,119]/255,'Linewidth',lw,'Markersize',ms)
    leg{end+1} = "Embedding, order: " + num2str(order_vec_GL(j));
    count = count + 1;
end

count = 1;
for j = 1:numel(order_vec_trad)

    plot(dofsvec_trad(:,j),errvec_trad(:,j),marker_types{count},'Color',[217,95,2]/255,'Linewidth',lw,'Markersize',ms)
    leg{end+1} = "Traditional, order: " + num2str(order_vec_trad(j));
    count = count + 1;
end

count = 1;
for j = 1:numel(order_vec_bopt)
    plot(dofsvec_bopt(:,j),errvec_bopt(:,j),marker_types{count},'Color',[117,112,179]/255,'Linewidth',lw,'Markersize',ms)
    leg{end+1} = "Boundary optimized, order: " + num2str(order_vec_bopt(j));
    count = count + 1;
end

xlabel('DOFs')
ylabel('$L_2$-error','interpreter','latex')
set(gca,'XScale','log','YScale','log')
set(gca,'Fontsize',fs)
legend(leg,'Location','southwest')
axis([2e2,2e5,1e-10,1e0])

exportgraphics(gca, "conv_plot.pdf", 'ContentType', 'vector');

leg = {};
figure('pos',[1000         275        1176        1062])
hold on
box on
grid on
marker_types = {'*-','x--','d:','+-.'};
count = 1;
for j = 1:numel(order_vec_GL)
    plot(time_elapsedvec_GL(:,j),errvec_GL(:,j),marker_types{count},'Color',[27,158,119]/255,'Linewidth',lw,'Markersize',ms)
    leg{end+1} = "Embedding, order: " + num2str(order_vec_GL(j));
    count = count + 1;
end

count = 1;
for j = 1:numel(order_vec_trad)
    plot(time_elapsedvec_trad(:,j),errvec_trad(:,j),marker_types{count},'Color',[217,95,2]/255,'Linewidth',lw,'Markersize',ms)
    leg{end+1} = "Traditional, order: " + num2str(order_vec_trad(j));
    count = count + 1;
end

count = 1;
for j = 1:numel(order_vec_bopt)
    plot(time_elapsedvec_bopt(:,j),errvec_bopt(:,j),marker_types{count},'Color',[117,112,179]/255,'Linewidth',lw,'Markersize',ms)
    leg{end+1} = "Boundary optimized, order: " + num2str(order_vec_bopt(j));
    count = count + 1;
end

xlabel('Runtime [s]')
ylabel('$L_2$-error','interpreter','latex')
set(gca,'XScale','log','YScale','log')
set(gca,'Fontsize',fs)
legend(leg,'Location','southwest')
% axis([1e2,2e5,1e-10,1e0])

exportgraphics(gca, "eff_plot.pdf", 'ContentType', 'vector');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Print tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
order_vec_GL = [5,7,9];
order_vec_trad = [2,4,6];
order_vec_bopt = [8,10,12];

for j = 1:numel(order_vec_GL)
    % for j = 1:numel(order_vec)
    qvec_GL = [];
    for i = 1:numel(hvec_GL)-1
        qvec_GL(i) = 2*log(errvec_GL(i,j)/errvec_GL(i+1,j))/log(dofsvec_GL(i,j)/dofsvec_GL(i+1,j));
    end
    % end
    qvec_GL = [0,qvec_GL];

    fprintf("\n--------- Order: %d -----------\n",order_vec_GL(j))
    fprintf("dofs\terr\tq\n")
    fprintf("%d\t%.2f\t-\n",dofsvec_GL(1,j),log10(errvec_GL(1,j)));
    for i = 2:numel(qvec_GL)
        fprintf("%d\t%.2f\t%.2f\n",dofsvec_GL(i,j),log10(errvec_GL(i,j)),qvec_GL(i));
    end
    fprintf("---------------------------------------------------\n")
end

% return
for j = 1:numel(order_vec_trad)
    qvec_trad = [];
    for i = 1:numel(mvec_trad)-1
        qvec_trad(i) = 2*log(errvec_trad(i,j)/errvec_trad(i+1,j))/log(dofsvec_trad(i,j)/dofsvec_trad(i+1,j));
    end
    % end
    qvec_trad = [0,qvec_trad];

    fprintf("\n--------- Traditional, order: %d -----------\n",order_vec_trad(j))
    fprintf("dofs\terr\tq\n")
    fprintf("%d\t%.2f\t-\n",dofsvec_trad(1,j),log10(errvec_trad(1,j)));
    for i = 2:numel(qvec_trad)
        fprintf("%d\t%.2f\t%.2f\n",dofsvec_trad(i,j),log10(errvec_trad(i,j)),qvec_trad(i));
    end
    fprintf("---------------------------------------------------\n")
end

for j = 1:numel(order_vec_bopt)
    qvec_bopt = [];
    for i = 1:numel(mvec_bopt)-1
        qvec_bopt(i) = 2*log(errvec_bopt(i,j)/errvec_bopt(i+1,j))/log(dofsvec_bopt(i,j)/dofsvec_bopt(i+1,j));
    end
    % end
    qvec_bopt = [0,qvec_bopt];

    fprintf("\n--------- Boundary optimized, order: %d -----------\n",order_vec_bopt(j))
    fprintf("dofs\terr\tq\n")
    fprintf("%d\t%.2f\t-\n",dofsvec_bopt(1,j),log10(errvec_bopt(1,j)));
    for i = 2:numel(qvec_bopt)
        fprintf("%d\t%.2f\t%.2f\n",dofsvec_bopt(i,j),log10(errvec_bopt(i,j)),qvec_bopt(i));
    end
    fprintf("---------------------------------------------------\n")
end