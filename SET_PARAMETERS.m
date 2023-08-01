%% Frequently changed parameters

grid.N=4*512;
grid.L = 160/2;

grid.T = 0.01;          %0.01;%0.001; 0.01; 0.1
T_relaxation = 3.8;     %4*4.5;%4.4*1;

cost = 'energy_state';
a = 0.2;    % 0.18;     % energy
b = 0.8;    % 0.82;     % state


lambda_0 = [0;0];
lambda_T = [0;0];

dim_lambda = 1000;      % dimension of lambda space

pgrid.nc_factor = 5;    % smallifier for p grid
