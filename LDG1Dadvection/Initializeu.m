function [U] = Initializeu(Env,M,N)
%function Initializeu output variable
%   U initial value
%   input variable
%   Env:element node values
%   M,N
%purpose:calculate initial value
% predistribution
U=zeros(3*N,1);
F=zeros(3*N,1);
% element loop
for i=1:N
    a=Env(i,1);
    b=Env(i,2);
    mid=(a+b)/2;
    
 % integrand   
f1=@(x)(abs(x-1));
f2=@(x)(abs(x-1).*(x-mid));
f3=@(x)(abs(x-1).*(x-mid).^2);

% integrate
quad1=quadrature(f1,a,b);
quad2=quadrature(f2,a,b);
quad3=quadrature(f3,a,b);
F(3*i-2:3*i)=[quad1,quad2,quad3]';
end
%[l,u]=lu(M);
U=M\F;
end


