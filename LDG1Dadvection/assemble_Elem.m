function [M,D,BR,BL] = assemble_Elem(Env,N,h)
%function assemble_Elem 
%   output variable
%   M:use basis function quadrature matrix
%   D:use basis function  derivative matrix
%   B:the value of basis function multiplication on boundary  
%   input variable
%   Env:element node values
%purpose:operator Assemble Element 

%predistribution
m = zeros(3,3);
d = zeros(3,3);
b = zeros(3,3);
bl = zeros(3,3);
br = zeros(3,3);

M=zeros(3*N,3*N);
D=zeros(3*N,3*N);
BR=zeros(3*N,3*N);
BL=zeros(3*N,3*N);

%N is element number，h is the value of element
%经过推到得到两个基函数相乘在单元内的右左端点的值为,只与h有关，无需遍历单元
br=[1,h/2,h^2/4;h/2,h^2/4,h^3/8;h^2/4,h^3/8,h^4/16];
bl=[1,h/2,h^2/4;-h/2,-h^2/4,-h^3/8;h^2/4,h^3/8,h^4/16];

%基函数的形式
% phi1=@(x)(1);
% phi2=@(x)(x-mid);
% phi3=@(x)((x-mid)^2);
% 求导后的基函数形式
% Dphi1=@(x)(0);
% Dphi2=@(x)(1);
% Dphi3=@(x)(2*(x-mid));

for i=1:N
    a=Env(i,1);
    b=Env(i,2);
    mid=(a+b)/2;
    % primary function multiplication
    f1=@(x)(1);
    f2=@(x)(x-mid);
    f3=@(x)((x-mid).^2);
    f4=@(x)((x-mid).^3);
    f5=@(x)((x-mid).^4);

    % primary function derivation multiplication
    y1=@(x)(1);
    y2=@(x)(x-mid);
    y3=@(x)((x-mid).^2);
    y4=@(x)(2*(x-mid));
    y5=@(x)(2*(x-mid).^2);
    y6=@(x)(2*(x-mid).^3);
    %calculate integ
    quad1 = quadrature(f1,a,b);
    quad2 = quadrature(f2,a,b);
    quad3 = quadrature(f3,a,b);
    quad4 = quadrature(f4,a,b);
    quad5 = quadrature(f5,a,b);
    %
    dquad1 = quadrature(y1,a,b);
    dquad2 = quadrature(y2,a,b);
    dquad3 = quadrature(y3,a,b);
    dquad4 = quadrature(y4,a,b);
    dquad5 = quadrature(y5,a,b);
    dquad6 = quadrature(y6,a,b);
    % operator matrix m
    m=[quad1,quad2,quad3;quad2,quad3,quad4;quad3,quad4,quad5];
    % operator matrix d
    d=[0,0,0;dquad1,dquad2,dquad3;dquad4,dquad5,dquad6];
    % assemble
    M(3*i-2:3*i,3*i-2:3*i)=m;
    D(3*i-2:3*i,3*i-2:3*i)=d;
    BR(3*i-2:3*i,3*i-2:3*i)=br;
    BL(3*i-2:3*i,3*i-2:3*i)=bl;
    end
end


