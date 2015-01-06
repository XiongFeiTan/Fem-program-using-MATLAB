function x = SonePossion(a,b,N)
%% 
% 
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%Record of revisions:
%Date       Programmer      Description of change
%=====    ==========     =================
%11/18     Artvigo               Original code 
%
%Define variables:
%   A:matrix of coefficient
%   F:constant matrix
%   a:lower limit of integral
%   b:upper limit of integral
%   x:result

%exception  handling
if (nargin<3)
        error(message('NotEnoughInputs'));
elseif (nargin>3)
     error(message('TooMuchInputs'));
end

%% Assemble stiffness matrix
   A = local2Global(a,b,N);
   
%% Assemble the right hand side and Set up boundary conditions
%N=50;
%a=0;
%b=1;
h=(b-a)/N;
nodes=linspace(a,b,N+1);
F = zeros(N+1,1);
for i=1:N 
 quadvalue1 = integ(nodes(i),nodes(i+1),1);
 quadvalue2 = integ(nodes(i),nodes(i+1),2);
F(i,1) = F(i,1)+quadvalue1;
F(i+1,1) = F(i+1,1)+ quadvalue2;
end


%%
%Set up boundary conditions
F(1,1)=0;
F(N,1)=F(N,1)-A(N+1,N)*1;
F(N+1,1)=1;
%F
%set matrix conditions
A(1,1)=1;
A(1,2)=0;
A(2,1)=0;
A(N+1,N+1)=1;
A(N+1,N)=0;
A(N,N+1)=0;

%A

%% Solve the system of linear equations
x=FollowUp(A,F);
format ;