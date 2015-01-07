function [A] = assemble(M,r_bl,r_br,q_bl,q_br,D,N)
%函数 assemble_R_RH 组装R方程里的矩阵
%   输入变量：U初始值，N为单元数，h为步长，D为求导后的基函数做内积的总刚矩阵
%                   R_BR,R_BL关于R方程边界右端点的总刚矩阵,左端点的总刚矩阵
%   输出变量：A
%预分配
A=zeros(6*N,6*N);
c1=zeros(3*N,3*N);
c2=zeros(3*N,3*N);
c3=zeros(3*N,3*N);
%c2
for i=1:N-1
    c2(3*i-2:3*i,3*i-2:3*i)=D(3*i-2:3*i,3*i-2:3*i)+q_bl;
    c2(3*i-2:3*i,3*i+1:3*i+3)=-q_br;
end
c2(3*N-2:3*N,3*N-2:3*N)=D(3*N-2:3*N,3*N-2:3*N)+q_bl;
c2(3*N-2:3*N,1:3)=-q_br;

%c3
for i=2:N
    c3(3*i-2:3*i,3*i-2:3*i)=D(3*i-2:3*i,3*i-2:3*i)-r_br;
    c3(3*i-2:3*i,3*i-5:3*i-3)=r_bl;
end
c3(1:3,1:3)=D(1:3,1:3)-r_br;
c3(1:3,3*N-2:3*N)=r_bl;

%组装耦合矩阵
c1=M;
A(1:3*N,1:3*N)=c1;
A(3*N+1:6*N,3*N+1:6*N)=c1;
A(1:3*N,3*N+1:6*N)=c2;
A(3*N+1:6*N,1:3*N)=c3;

end
