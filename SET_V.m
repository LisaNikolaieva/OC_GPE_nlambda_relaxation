Vmax        = 2e3;
sigma       = 2;
r0 = 40;
V0      =  (Vmax - Vmax * 1/2 * (1 + erf((grid.x + r0)/sigma)) + Vmax * 1/2 * ( 1+ erf((grid.x - r0)/sigma)));

Theta = (-1)*(grid.x<0)+(grid.x>0);
Theta2 = exp(-(grid.x).^2);
Theta2 = 1 - Theta2./max(abs(Theta2));

V      = @(lambda) Vmax*Theta*lambda(1) + Vmax*Theta2*lambda(2) + V0;

dV_dlambda_fun = @(lambda,lambda_indx) 1/(1e-3)*(V(lambda+[zeros(lambda_indx-1,1);5e-4;zeros(numel(lambda)-lambda_indx,1)])-V(lambda-[zeros(lambda_indx-1,1);5e-4;zeros(numel(lambda)-lambda_indx,1)]));
