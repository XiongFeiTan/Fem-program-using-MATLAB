%%
a=0;
b=2;
N=24;
%网格
[h,Env] = MeshGen1D(a,b,N);
%装配局部矩阵
[M,D,R_BR,R_BL,Q_BR,Q_BL] = assemble_Elem(Env,N,h);
%初始化右端项
[Q] = Initializeu(Env,N);
%组装耦合矩阵
[A] = assemble(M,R_BR,R_BL,Q_BR,Q_BL,D,N);
%%
%求解线性方程组
X=zeros(6*N,1);
U=zeros(3*N,1);
R=zeros(3*N,1);

X=A\Q;
U(1:3*N,1)=X(1:3*N,1);
%%
x = linspace (a, b, N+1);

nu=zeros(N,1);
eu=zeros(N,1);

% 数值解 nu
for i=1:N
    nu(i,1)=U(3*i-2,1);
end

% 真解 eu
    for i= 1:N
        x(i)=(x(i)+x(i+1))/2;
        eu(i,1) = sin(2*pi*x(i));        
    end
% eu


% 误差估计
error=abs(nu-eu)


% L1_error=sum(abs(nu-eu))/N 
% L2_error=sqrt(sum((nu-eu).^2)/N)
L2_error=norm(eu-nu,2)
Lmax_error=max(abs(nu-eu))

% 画图
plot ( x(1:N),nu,'.b:', x(1:N),eu,'+r:');
% plot ( x(1:N),nu,'.b:'); 
% plot ( x(1:N),error,'+r:');
legend ( 'numrical solution ', 'exact    solution ') ;
grid on ;
u_min=min(min(nu));
u_max=max(max(nu));
axis([0 2 u_min-1,u_max+1]);
title ( sprintf ( ' NStep%d', N ) );
xlabel ( '------X------' ) ;
ylabel ( '------U------' );


