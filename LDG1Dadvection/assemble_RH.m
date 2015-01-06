function [Ut] = assemble_RH(U,UL,M,BR,BL,D,N,h,dt)
%function assemble_RH 
%   output variable
%   input variable
%   purpose:calculate the right hand sides
%predistribution
F=zeros(3*N,3*N);
%right hand side calculate
F=M*U+dt*((D-BR)*U+BL*UL);
%[l,u]=lu(M);
Ut=M\F;
end
