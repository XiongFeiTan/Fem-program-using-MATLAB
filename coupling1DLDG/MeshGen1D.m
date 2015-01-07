function [h,Env] = MeshGen1D(min,max,N)
%N £ºelement number
K = N+1; 
% element node coordinate
Env = zeros(N, 2);
NodeX = (1:K);
for i = 1:K
 NodeX(i) = (max-min)*(i-1)/N + min;
end
for k = 1:N
  Env(k,1) = NodeX(k); 
  Env(k,2) = NodeX(k+1);
end
h=(max-min)/N;
return
