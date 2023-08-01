close all;
Config;

%% Optimization

f_fun = @(X) cost_function(X,cf_par); % f:      X |-> R
gradf_fun = @(X) gradJ(X,gradJ_par);  % gradf:  X |-> x
X0 = lambda0;

n_steps = 3;

[lambda_store, cost_function_store] = find_min_BFGS(f_fun,gradf_fun,X0,X2x,x2X,n_steps);

%%
Postprocessing;