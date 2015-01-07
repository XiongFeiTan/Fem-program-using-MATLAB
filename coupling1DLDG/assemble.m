function [A] = assemble(M,r_bl,r_br,q_bl,q_br,D,N)
%���� assemble_R_RH ��װR������ľ���
%   ���������U��ʼֵ��NΪ��Ԫ����hΪ������DΪ�󵼺�Ļ��������ڻ����ܸվ���
%                   R_BR,R_BL����R���̱߽��Ҷ˵���ܸվ���,��˵���ܸվ���
%   ���������A
%Ԥ����
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

%��װ��Ͼ���
c1=M;
A(1:3*N,1:3*N)=c1;
A(3*N+1:6*N,3*N+1:6*N)=c1;
A(1:3*N,3*N+1:6*N)=c2;
A(3*N+1:6*N,1:3*N)=c3;

end