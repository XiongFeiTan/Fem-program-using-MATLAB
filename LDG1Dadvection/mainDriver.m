% step1:写出基函数
% step2:求出三个算子M,D,B
% step3:求出M的逆矩阵
% step4:形成常微分方程：具有初值的常微分方程，进行时间的离散。
%%
N=10;
a=0;
b=2;
%产生网格
[h,Env] = MeshGen1D(a,b,N);
%计算单元
[M,D,BR,BL] = assemble_Elem(Env,N,h);
M
%计算初值
[U] = Initializeu(Env,M,N);
%%
 %边界处理,跟循环的U无关的项
 UL=zeros(3*N,1);
 UL(1:3,1)=[1,1,1]';
 BL(1:3,1:3)=zeros(3,3);
 BL(1,1)=1;
 BL(2,2)=-h/2;
 BL(3,3)=h^2/4;
 %U
 %BL
%%
%时间步长
time=1;
dt=0.001;
%时间离散
for tstep=1:1000
    %边界处理，跟循环U有关的项
    for i=2:N
        UL(3*i-2:3*i,1)=U(3*i-5:3*i-3,1);
    end
   % UL
    %计算右端项
    %U
    [Ut] = assemble_RH(U,UL,M,BR,BL,D,N,h,dt);
    U=Ut;
end
%U;
%%
x = linspace (0, 2, N+1);
t = linspace (0.0,1.0,100) ;
nu=zeros(N,1);
eu=zeros(N,1);

% select numerical solution nu
for i=1:N
    nu(i,1)=U(3*i-2);
end
%nu

% eaxct solution 
    for i= 1:N
        x(i)=(x(i)+x(i+1))/2;
        if x(i)>time +1
            eu(i,1) = x(i)-time -1.0;
        elseif x(i)<=time+1&&x(i) >time
            eu(i,1) =1-x(i)+time;
        elseif x(i)<=time
            eu(i,1) =1;
        end
    end
%eu
%%
%   error estimate
error=abs(nu-eu)
L1_error=sum(abs(nu-eu))/N
L2_error=sqrt(sum((nu-eu).^2)/N)
L3_error=max(abs(nu-eu))
%  draw 
plot ( x(1:N),nu,'.b:', x(1:N),eu,'+r:');
legend ( 'numrical solution ', 'exact    solution ') ;
grid on ;
u_min=min(min(nu));
u_max=max(max(nu));
axis([0 2 u_min-1.0,u_max+1.0]);
title ( sprintf ( ' NStep%d,Time%f \n ', N,1) );
xlabel ( '------X------' ) ;
ylabel ( '------U------' );

