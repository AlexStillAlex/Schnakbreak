function Gu=sGjac(p,u)
    par=u(p.nu+1:end); % parameters
    n=p.np;
    [f1u,f1v,f2u,f2v]=njac(p,u,par); % the Jacobian, see below
    Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
       [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
   
   %par(5) = L, par(4) = D
    Gu=kron([[1/(par(5))^2,0];[0,par(4)/(par(5)^2)]],p.mat.K)-p.mat.M*Fu; % assemble the Jacobian. 
    %D = (1/L^2 0) (diffusion matrix)
        %(0 D/L^2)
end 
function [f1u,f1v,f2u,f2v]=njac(p,u,par) % Jacobian for Schnakenberg
    u1=u(1:p.np); % solution component 1
    u2=u(p.np+1:2*p.np); % solution component 2
    % entries of the jacobian
    f1u=-par(3)+2*u1.*u2; %schnak = a - cu +u^2v; b - u^2v
    f1v=u1.^2; 
    f2u=-2*u1.*u2; 
    f2v=-u1.^2-2; 
end