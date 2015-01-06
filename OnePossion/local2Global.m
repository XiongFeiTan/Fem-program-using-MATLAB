function M = local2Global(a,b,N)
%局部到全局，dprifun 这是对基函数求导求积后的值。
 h =((b-a)/N);
dprifun1= 1/h;
dprifun2= -1/h;
dprifun3=-1/h;	
dprifun4 =1/h;
M=sparse(N+1,N+1);
%循环的单元
		for i = 1:1: N
				M(i,i) =dprifun1+M(i,i);
				M(i,i+1) = dprifun2;
				M(i+1,i) = dprifun3;
				M(i+1,i+1) = dprifun4+M(i+1,i+1);
        end
% 处理边界      
M=full(M);
  



