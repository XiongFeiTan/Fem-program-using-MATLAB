function [M,D,r_bl,r_br,q_bl,q_br] = assemble_Elem(Env,N,h)
%函数 Assemble_Elem 组装方程里的刚度矩阵
%   输入变量：Env单元节点值，N为单元数，h为步长
%   输出变量：M为基函数做内积的总刚矩阵，D为求导后的基函数做内积的总刚矩阵
%                    R_BR,R_BL关于R方程边界右端点的总刚矩阵,左端点的总刚矩阵

%预分配
m = zeros(3,3);
d = zeros(3,3);
bl = zeros(3,3);
br = zeros(3,3);
%总刚矩阵
M=zeros(3*N,3*N);
D=zeros(3*N,3*N);
% R_BR=zeros(3*N,3*N);
% R_BL=zeros(3*N,3*N);
% Q_BR=zeros(3*N,3*N);
% Q_BL=zeros(3*N,3*N);

%经过推到得到两个基函数相乘在单元内的右左端点的值，只与h有关，无需遍历单元
r_br=[1,h/2,h^2/4;h/2,h^2/4,h^3/8;h^2/4,h^3/8,h^4/16];
r_bl=[1,h/2,h^2/4;-h/2,-h^2/4,-h^3/8;h^2/4,h^3/8,h^4/16];

q_br=[1,-h/2,h^2/4;h/2,-h^2/4,h^3/8;h^2/4,-h^3/8,h^4/16];
q_bl=[1,-h/2,h^2/4;-h/2,h^2/4,-h^3/8;h^2/4,-h^3/8,h^4/16];

%基函数的形式
% phi1=@(x)(1);
% phi2=@(x)(x-mid);
% phi3=@(x)((x-mid)^2);
% 求导后的基函数形式
% Dphi1=@(x)(0);
% Dphi2=@(x)(1);
% Dphi3=@(x)(2*(x-mid));

%遍历单元，组装成总刚矩阵
for i=1:N
    %a为单元左端点
    a=Env(i,1);
    %b为单元右端点
    b=Env(i,2);
    mid=(a+b)/2;
    %两基函数的乘积有五种情况，分别如下
    f1=@(x)(1);
    f2=@(x)(x-mid);
    f3=@(x)((x-mid).^2);
    f4=@(x)((x-mid).^3);
    f5=@(x)((x-mid).^4);

    %基函数求导后的乘积有六种情况，分别如下
    y1=@(x)(1);
    y2=@(x)(x-mid);
    y3=@(x)((x-mid).^2);
    y4=@(x)(2*(x-mid));
    y5=@(x)(2*(x-mid).^2);
    y6=@(x)(2*(x-mid).^3);
    %利用求积函数进行求积运算
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
    % 单刚矩阵m
    m=[quad1,quad2,quad3;quad2,quad3,quad4;quad3,quad4,quad5];
    % 单刚矩阵d
    d=[0,0,0;dquad1,dquad2,dquad3;dquad4,dquad5,dquad6];
    % 组装
    M(3*i-2:3*i,3*i-2:3*i)=m;
    D(3*i-2:3*i,3*i-2:3*i)=d;
    end
end




