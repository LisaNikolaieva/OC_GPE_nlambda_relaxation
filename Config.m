%%
SET_PARAMETERS;


%%
close all
Physical_parameters;


%% grid
grid.x = linspace(-grid.L/2,grid.L/2,grid.N)';
grid.dx = grid.x(2)-grid.x(1);
grid.lap4 = gallery( 'toeppen', grid.N, - 1,  16, - 30, 16, - 1 ) / ( 12 * grid.dx ^ 2 );
grid.wav = 2 * pi * ( 0 : grid.N - 1 )' / grid.N;
grid.ilap = - 2 * ( 1 - cos( grid.wav ) ) / grid.dx ^ 2;


%% define potential
SET_V;

%% calculate initial and final states

SET_u0_uT;

den_0 = abs(u0).^2;
den_T = abs(uT).^2;
figure
plot(grid.x,den_T)
hold on
plot(grid.x,0.1*max(den_T)*ones(size(den_T)),'--')
drawnow

c_s1 = sqrt(par.g*max(abs(uT).^2)/par.m); % estimated speed of sound for final state
c_s0 = sqrt(par.g*max(abs(u0).^2)/par.m); % estimated speed of sound for initial state




%% time grid
dt = min(0.1*grid.dx/max(c_s1,c_s0),grid.T/dim_lambda);
grid.Nt = ceil(grid.T/dt);
grid.t = linspace(0,grid.T,grid.Nt);
grid.dt = grid.t(2) - grid.t(1);
[grid.t_mesh,grid.x_mesh] = meshgrid(grid.t,grid.x);


%% pgrid
pgrid.L = grid.L;
pgrid.N = grid.N;
pgrid.x = grid.x;
pgrid.dx = grid.dx;
pgrid.lap4 = grid.lap4;
% pgrid.nc_factor = 5;
pgrid.T = grid.T;
pgrid.Nt = ceil(grid.Nt/pgrid.nc_factor);
pgrid.t = linspace(0,pgrid.T,pgrid.Nt);
pgrid.dt = pgrid.t(2) - pgrid.t(1);         
pgrid.lap4_t   =  gallery( 'toeppen', pgrid.Nt, - 1,  16, - 30, 16, - 1 ) / ( 12 * pgrid.dt ^ 2 );
[pgrid.t_mesh,pgrid.x_mesh] = meshgrid(pgrid.t,pgrid.x);


%% grid2 & pgrid2

grid2.L = grid.L;
grid2.N=grid.N;
grid2.x = grid.x;
grid2.dx = grid.dx;
grid2.lap4 = grid.lap4;
grid2.ilap = grid.ilap;

dt = 0.1*grid2.dx/max(c_s1,c_s0);
grid2.Nt = ceil((T_relaxation-grid.T)/dt);
grid2.t = linspace(grid.T,T_relaxation,grid2.Nt);
grid2.dt = grid2.t(2) - grid2.t(1);
[grid2.t_mesh,grid2.x_mesh] = meshgrid(grid2.t,grid2.x);

pgrid2.L = grid2.L;
pgrid2.N = grid2.N;
pgrid2.x = grid2.x;
pgrid2.dx = grid2.dx;
pgrid2.lap4 = grid2.lap4;
pgrid2.nc_factor = 5;
pgrid2.Nt = ceil(grid2.Nt/pgrid2.nc_factor);
pgrid2.t = linspace(grid.T,T_relaxation,pgrid2.Nt);
pgrid2.dt = pgrid2.t(2) - pgrid2.t(1);         
pgrid2.lap4_t   =  gallery( 'toeppen', pgrid2.Nt, - 1,  16, - 30, 16, - 1 ) / ( 12 * pgrid2.dt ^ 2 );
[pgrid2.t_mesh,pgrid2.x_mesh] = meshgrid(pgrid2.t,pgrid2.x);




%% cost_fun_fun 
ham       = @(lambda_t0)-0.5*grid.lap4/par.m+spdiag(V(lambda_t0));
mu_fun = @(lambda_t0,Psi_xt0) real(...
   sum(( -0.5*grid.lap4/par.m + spdiag(V(lambda_t0)+par.g*(abs(Psi_xt0).^2)) )*Psi_xt0.*conj(Psi_xt0) )...
   /sum(abs(Psi_xt0.^2)) );



energy_fun = @(Psi_xt0,lambda_t0) real(trapz(grid.x,conj(Psi_xt0).*(ham(lambda_t0)+par.g/2*spdiag(abs(Psi_xt0).^2))*Psi_xt0));
state_fun = @(Psi_xT) 1/2*(1-abs(trapz(grid.x,conj(uT).*Psi_xT)).^2);
energy_state_fun = @(Psi_xT) a*(energy_fun(Psi_xT,lambda_T)-energy_fun(uT,lambda_T)) + b*state_fun(Psi_xT) ;



switch cost
    case 'energy'
        phi = @(Psi_xT) energy_fun(Psi_xT,lambda_T)-energy_fun(uT,lambda_T);
        pT_fun = @(Psi_xT) -2*1i*(ham(lambda_T) - spdiag(mu_fun(lambda_T,Psi_xT)*ones(size(Psi_xT))) + par.g*spdiag(abs(Psi_xT).^2))*Psi_xT;
    case 'state'
        phi = @(Psi_xT) state_fun(Psi_xT);
        pT_fun = @(Psi_xT) 1i * trapz(grid.x,conj(uT).*Psi_xT)*uT;
    case 'energy_state'
        phi = @(Psi_xT) energy_state_fun(Psi_xT);
        pT_fun = @(Psi_xT) a*(-2*1i*(ham(lambda_T) - spdiag(mu_fun(lambda_T,Psi_xT)*ones(size(Psi_xT))) + par.g*spdiag(abs(Psi_xT).^2))*Psi_xT)...
            + b*(1i * trapz(grid.x,conj(uT).*Psi_xT)*uT);
end


stab = @(lambda_t) par.gamma/2*sum(abs(grid.dt)*trapz((diff(lambda_t,2)/abs(grid.dt)).^2,2));

cost_fun_fun = @(Psi_xT,lambda_t) [phi(Psi_xT) + stab(lambda_t); ...
                                                    phi(Psi_xT); ...
                                                stab(lambda_t)];


                                            
int_vec_fun = @(pPsi_xt,pdVdl_xt,p_xt) trapz(pgrid.x,real(conj(pPsi_xt).*pdVdl_xt.*p_xt),1);

%% initialize lambda
n_l = numel(lambda_0);
SET_lambda0;

%% pack parameters

workspace.grid = grid;
workspace.pgrid = pgrid;
workspace.grid2 = grid2;
workspace.pgrid2 = pgrid2;
workspace.par = par;
workspace.cost_fun_fun = cost_fun_fun;
workspace.pT_fun = pT_fun;
workspace.int_vec_fun = int_vec_fun;
workspace.V = V;
workspace.dV_dlambda_fun = dV_dlambda_fun;
workspace.u0 = u0;


cf_par = workspace;
gradJ_par = workspace;



%% Interpalation functions for Psi & p grids X <-> x 

X2x = @(X) X2x_matrix(X,grid,pgrid); % X2x:  X |-> x
x2X = @(x) x2X_matrix(x,grid,pgrid); % x2X:  x |-> X


function x = X2x_matrix(X,grid,pgrid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:            output:
% X = 1 1 1 1 1     x = 1 1 1 
%     2 2 2 2 2         2 2 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = pgrid.Nt;
n_l = size(X,1);

x = zeros(n_l,n);

for i = 1:n_l
   x(i,:) = interp1(grid.t,X(i,:),pgrid.t,'linear','extrap');
end


end


function X = x2X_matrix(x,grid,pgrid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:         output:
% x = 1 1 1      X = 1 1 1 1 1
%     2 2 2          2 2 2 2 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = grid.Nt;
n_l = size(x,1);

X = zeros(n_l,N);

for i = 1:n_l
   X(i,:) = interp1(pgrid.t,x(i,:),grid.t,'linear','extrap');
end


end
