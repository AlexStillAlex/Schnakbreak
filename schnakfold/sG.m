function r=sG(p,u)
    u1=u(1:p.np); % solution component 1
    u2=u(p.np+1:2*p.np); % solution component 2
    par=u(p.nu+1:end); % parameters
    f1=par(1)-par(3)*u1+u1.^2.*u2; % F_1(u), see eqn (1) in tut
    f2=par(2)-u1.^2.*u2; % F_2(u)
    f=[f1;f2];
    K=kron([[1/(par(5))^2,0];[0,par(4)/(par(5)^2)]],p.mat.K); % assemble full FEM matrix
    r=K*[u1;u2]-p.mat.M*f; % the residual
end


