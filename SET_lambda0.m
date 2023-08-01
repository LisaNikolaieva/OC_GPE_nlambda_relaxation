lambda0 = zeros(n_l,grid.Nt);
for i = 1:n_l
    lambda0(i,:) = -(grid.t-grid.T/2).^2+(grid.T/2)^2;
    lambda0(i,:) = lambda0(i,:)./max(abs(lambda0(i,:)));
end
