function gradJ = gradJ(lambda,gradJ_par)

% lambda = 1 1 1 1 1    
%          2 2 2 2 2 
n_l = size(lambda,1);

grid = gradJ_par.grid;
pgrid = gradJ_par.pgrid;
grid2 = gradJ_par.grid2;
pgrid2 = gradJ_par.pgrid2;
par = gradJ_par.par;
V = gradJ_par.V;
dV_dlambda_fun = gradJ_par.dV_dlambda_fun;
pT_fun = gradJ_par.pT_fun;
int_vec_fun = gradJ_par.int_vec_fun;
u0 = gradJ_par.u0;


%%
[Psi_store] = Psi_xt(u0,grid,par,V,lambda);

%         figure(21)
%         plot(grid.x,abs(Psi_store(:,end)).^2)
%         set(gca,'XMinorGrid','on');
%         set(gca,'YMinorGrid','on');


% plambda = interp1(grid.t,lambda,pgrid.t,'linear','extrap');
plambda = zeros(n_l,pgrid.Nt);
for i = 1:n_l
   plambda(i,:) = interp1(grid.t,lambda(i,:),pgrid.t,'linear','extrap');
end

pPsi_store = interp2(grid.t_mesh,grid.x_mesh,Psi_store,pgrid.t_mesh,pgrid.x_mesh,"nearest",0);                       % psi_p = interp2(grid.t_mesh,grid.x_mesh,Psi_store(:,2:end),pgrid.t_mesh,pgrid.x_mesh,"nearest",0);

% pT = pT_fun(Psi_store(:,end));
pT = calc_pT(Psi_store(:,end),grid2,pgrid2,par,V,pT_fun,n_l);

[p_store] = p_xt(pT,pgrid,par,V,plambda,pPsi_store);

%%
int_vec = zeros(n_l,pgrid.Nt);
pdV_dlambda = zeros(pgrid.N,pgrid.Nt);


for i = 1:n_l
    for j = 1:pgrid.Nt
        pdV_dlambda(:,j) = dV_dlambda_fun(plambda(:,j),i);
    end
    int_vec(i,:) = int_vec_fun(pPsi_store,pdV_dlambda,p_store);
end

        figure(1)
        subplot(4,1,2)
        hold on
        plot(pgrid.t,trapz(grid.x,real(pdV_dlambda),1))
        subplot(4,1,4)
        hold on
        plot(pgrid.t,trapz(grid.x,real(conj(pPsi_store).*pdV_dlambda.*p_store),1))

        subplot(4,1,1)
        plot(pgrid.t,trapz(grid.x,real(conj(pPsi_store)),1))
        hold on
        subplot(4,1,3)
        plot(pgrid.t,trapz(pgrid.x,real(p_store),1))
        hold on

gradJ = par.gamma*plambda' + (pgrid.lap4_t\(int_vec'));
gradJ = gradJ./max(abs(gradJ),[],1);



end



function p0 = calc_pT(u0,grid2,pgrid2,par,V,pT_fun,n_l)

lambda = zeros(n_l,grid2.Nt);
[Psi_store] = Psi_xt(u0,grid2,par,V,lambda);
% plambda = interp1(grid2.t,lambda,pgrid2.t,'linear','extrap');
plambda = zeros(n_l,pgrid2.Nt);
pPsi_store = interp2(grid2.t_mesh,grid2.x_mesh,Psi_store,pgrid2.t_mesh,pgrid2.x_mesh,"nearest",0);                       % psi_p = interp2(grid.t_mesh,grid.x_mesh,Psi_store(:,2:end),pgrid.t_mesh,pgrid.x_mesh,"nearest",0);
pT = pT_fun(Psi_store(:,end));
[p_store] = p_xt(pT,pgrid2,par,V,plambda,pPsi_store);

p0 = p_store(:,1);


        figure(22)
        subplot(2,1,1)
        plot(grid2.x,abs(Psi_store(:,end)).^2)
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');
        subplot(2,1,2)
        plot(grid2.x,angle(Psi_store(:,end)))
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');




end
