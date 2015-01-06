function [x,step] = congradient(A,b,iv)
%function filename : congradient.m
%Purpose:
%To use the conjugate gradient method
%
%Record of revisions:
%Date       Programmer      Description of change
%=====    ==========     =================
%11/15     Artvigo               Original code 
%
%Define variables:
%   A:matrix of coefficient
%   b:constant matrix
%   iv:initialization vector
%   eps:precision
%   x:result
%   step:how many steps convergence

%exception  handling
if (nargin<3)
        error(message('MATLAB:congradient:NotEnoughInputs'));
elseif (nargin>3)
     error(message('MATLAB:congradient:TooMuchInputs'));
end

%initializaiton variable
r0 = b - A*iv;
p0 = r0;
alpha = dot(r0,r0)/dot(p0,A*p0);
x1 = iv + alpha*p0;
r1 = r0 - alpha*A*p0;
beta = dot(r1,r1)/dot(r0,r0);
p1 = r1 + beta*p0;

%loop start 
for n=(1:rank(A)-1)
    if (dot(p1,A*p1)<1.0e-20)
        break;
    end
    x0 = x1;
    p0 = p1;
    r0 = r1;
    
    alpha = dot(r0,r0)/dot(p0,A*p0);
    x1 = x0 + alpha*p0;
    r1 = r0 - alpha*A*p0;
    if(r1 == 1.0e-20)
        break;
    end
    beta = dot(r1,r1)/dot(r0,r0);
    p1 = r1 + beta*p0;
    n = n+1;
end
    alpha = dot(r1,r1)/dot(p1,A*p1);
    x = x1 + alpha*p1;
    step = n;
end



