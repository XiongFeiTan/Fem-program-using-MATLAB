function [ql,Ak,xk]=guasslegendre(fun,a,b,n,tol)
% 高斯-勒让德数值积分
%
% 参数说明
% fun：积分表达式，可以是函数句柄、inline函数、匿名函数、字符串表达式，但是必须可以接受矢量输入
% a,b：积分上下限，注意积分区间太大会降低精度，此时建议使用复化求积公式，默认[-1 1]
% n：积分阶数，可以任意正整数，但是不建议设置过大，大不一定能得到更好的精度，默认7
% tol：积分精度，默认1e-6
% ql：积分结果
% Ak：求积系数
% xk：求积节点，满足ql=sum(Ak.*fun(xk))
%
% 举例说明
% fun=@(x)exp(x).*cos(x); % 必须可以接受矢量输入
% quadl(fun,0,pi) % 调用MATLAB内部积分函数检验
% [ql,Ak,xk]=guasslegendre(fun,0,pi)

if nargin==1
    a=-1;b=1;n=7;tol=1e-8;
elseif nargin==3
    n=7;tol=1e-8;
elseif nargin==4
    tol=1e-8;
elseif nargin==2|nargin>5
    error('The Number of Input Arguments Is Wrong!');
end
% 计算求积节点
syms x
p=sym2poly(diff((x^2-1)^(n+1),n+1))/(2^n*factorial(n));
tk=roots(p); % 求积节点
% 计算求积系数
Ak=zeros(n+1,1);
for i=1:n+1
    xkt=tk;
    xkt(i)=[];
    pn=poly(xkt);
    fp=@(x)polyval(pn,x)/polyval(pn,tk(i));
    Ak(i)=quadl(fp,-1,1,tol); % 求积系数
end
% 积分变量代换，将[a,b]变换到[-1,1]
xk=(b-a)/2*tk+(b+a)/2;
% 检验积分函数fun有效性
fun=fcnchk(fun,'vectorize');
% 计算变量代换之后积分函数的值
fx=fun(xk)*(b-a)/2;
% 计算积分值
ql=sum(Ak.*fx);