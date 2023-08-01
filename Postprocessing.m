precalculated = 0;
if precalculated
    Config;
    load('lambda_store_cost_function_store.mat','lambda_store','cost_function_store');
end

%%

lambda_init = lambda_store(:,:,1);
lambda_opt = lambda_store(:,:,end);

n_l = size(lambda_store,1);    

[Psi_store_init] = Psi_xt(u0,grid,par,V,lambda_init);
[Psi_store_init2] = Psi_xt(Psi_store_init(:,end),grid2,par,V,zeros(n_l,grid2.Nt));

[Psi_store_opt] = Psi_xt(u0,grid,par,V,lambda_opt);
[Psi_store_opt2] = Psi_xt(Psi_store_opt(:,end),grid2,par,V,zeros(n_l,grid2.Nt));



%%
pd = max(max(max(abs(Psi_store_opt).^2)),max(max(abs(Psi_store_opt2).^2)));

%%

figure(27)

subplot(2,2,1)
imagesc([grid2.t],grid.x,[abs(Psi_store_init2).^2])
caxis([0,pd])
ax = gca;
ax.FontSize=10;
ax.LabelFontSizeMultiplier = 1.5;
ax.TickLabelInterpreter='latex';
title('linear','FontSize',16,'Interpreter','latex')
ylabel('x, $\mu$ m','FontSize',16,'Interpreter','latex');



subplot(2,2,3)
imagesc([grid2.t],grid.x,[angle(Psi_store_init2)])
caxis([0,pd])
ax = gca;
ax.FontSize=10;
ax.LabelFontSizeMultiplier = 1.5;
ax.TickLabelInterpreter='latex';
xlabel('t, ms','FontSize',16,'Interpreter','latex');
ylabel('x, $\mu$ m','FontSize',16,'Interpreter','latex');


subplot(2,2,2)
imagesc([grid2.t],grid.x,[abs(Psi_store_opt2).^2])
caxis([0,pd])
ax = gca;
ax.FontSize=10;
ax.LabelFontSizeMultiplier = 1.5;
ax.TickLabelInterpreter='latex';
title('optimized','FontSize',16,'Interpreter','latex')

subplot(2,2,4)
imagesc([grid2.t],grid.x,[angle(Psi_store_opt2)])
caxis([0,pd])
ax = gca;
ax.FontSize=10;
ax.LabelFontSizeMultiplier = 1.5;
ax.TickLabelInterpreter='latex';
xlabel('t, ms','FontSize',16,'Interpreter','latex');


%%
if 1
%% plot density
    for i = 1:50:grid.Nt
        
        figure(26)
        subplot(2,2,1)
        plot(grid.x,abs(Psi_store_init(:,i)).^2,grid.x,den_T,'--',grid.x,den_0,'--')
        grid on
        ylim([0, pd])
        ax = gca;
        ax.FontSize=10;
        ax.LabelFontSizeMultiplier = 1.5;
        ax.TickLabelInterpreter='latex';
        title(sprintf('linear: t = %.3f ms',grid.t(i)),'FontSize',16,'Interpreter','latex')
        ylabel('$\rho$','FontSize',16,'Interpreter','latex');

        drawnow
        subplot(2,2,3)
        plot(grid.x,angle(Psi_store_init(:,i)))
        grid on
        ylim([-4, 4])
        ax = gca;
        ax.FontSize=10;
        ax.LabelFontSizeMultiplier = 1.5;
        ax.TickLabelInterpreter='latex';
        xlabel('x, $\mu$m','FontSize',16,'Interpreter','latex');
        ylabel('arg $\Psi$','FontSize',16,'Interpreter','latex');


        subplot(2,2,2)
        plot(grid.x,abs(Psi_store_opt(:,i)).^2,grid.x,den_T,'--',grid.x,den_0,'--')
        grid on
        ylim([0, pd])
        ax = gca;
        ax.FontSize=10;
        ax.LabelFontSizeMultiplier = 1.5;
        ax.TickLabelInterpreter='latex';
        title(sprintf('optimized: t = %.3f ms',grid.t(i)),'FontSize',16,'Interpreter','latex')
       
        drawnow
        subplot(2,2,4)
        plot(grid.x,angle(Psi_store_opt(:,i)) )
        grid on
        ylim([-4, 4])
        ax = gca;
        ax.FontSize=10;
        ax.LabelFontSizeMultiplier = 1.5;
        ax.TickLabelInterpreter='latex';
        xlabel('x, $\mu$m','FontSize',16,'Interpreter','latex');

    end
    
    
    
for i = 1:20:grid2.Nt
    
        figure(26)
        subplot(2,2,1)
        plot(grid2.x,abs(Psi_store_init2(:,i)).^2,grid2.x,den_T,'--',grid2.x,den_0,'--')
        grid on
        ylim([0, pd])
        ax = gca;
        ax.FontSize=10;
        ax.LabelFontSizeMultiplier = 1.5;
        ax.TickLabelInterpreter='latex';
        title(sprintf('linear: t = %.3f ms',grid2.t(i)),'FontSize',16,'Interpreter','latex')
        ylabel('$\rho$','FontSize',16,'Interpreter','latex');
        
        drawnow
        subplot(2,2,3)
        plot(grid2.x,angle(Psi_store_init2(:,i)))
        grid on
        ylim([-4, 4])
        ax = gca;
        ax.FontSize=10;
        ax.LabelFontSizeMultiplier = 1.5;
        ax.TickLabelInterpreter='latex';
        xlabel('x, $\mu$m','FontSize',16,'Interpreter','latex');
        ylabel('arg $\Psi$','FontSize',16,'Interpreter','latex');


        subplot(2,2,2)
        plot(grid2.x,abs(Psi_store_opt2(:,i)).^2,grid2.x,den_T,'--',grid2.x,den_0,'--')
        grid on
        ylim([0, pd])
        ax = gca;
        ax.FontSize=10;
        ax.LabelFontSizeMultiplier = 1.5;
        ax.TickLabelInterpreter='latex';
        title(sprintf('optimized: t = %.3f ms',grid2.t(i)),'FontSize',16,'Interpreter','latex')
        drawnow
        
        subplot(2,2,4)
        plot(grid2.x,angle(Psi_store_opt2(:,i)) )
        grid on
        ylim([-4, 4])
        ax = gca;
        ax.FontSize=10;
        ax.LabelFontSizeMultiplier = 1.5;
        ax.TickLabelInterpreter='latex';
        xlabel('x, $\mu$m','FontSize',16,'Interpreter','latex');


end

end

%%
figure(28);

minb = min(min(min(lambda_init)),min(min(lambda_opt)));
maxb = max(max(max(lambda_init)),max(max(lambda_opt)));
for i = 1:n_l
subplot(1,n_l,i)

plot(lambda_init(i,:),'--')
xlim([-inf inf])
ylim([minb, maxb])
hold on

plot(lambda_opt(i,:))
xlim([-inf inf])
ylim([minb, maxb])
ax = gca;
ax.FontSize=10;
ax.LabelFontSizeMultiplier = 1.5;
ax.TickLabelInterpreter='latex';
xlabel('t, ms','FontSize',16,'Interpreter','latex');
title(['$\lambda$', num2str(i)],'FontSize',16,'Interpreter','latex')


end
% drawnow



%%
fprintf('cost function reduction factor: %f\n',cost_function_store(1)/cost_function_store(end))



