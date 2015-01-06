N=50;
a=0;
b=1;
h=(b-a)/N;
nodes=linspace(a,b,N+1);
F = zeros(N+1,1);
for i=1:N 
 quadvalue1 = integ(nodes(i),nodes(i+1),1);
 quadvalue2 = integ(nodes(i),nodes(i+1),2);
F(i,1) = F(i,1)+quadvalue1;
F(i+1,1) = F(i+1,1)+ quadvalue2;
end
F;
