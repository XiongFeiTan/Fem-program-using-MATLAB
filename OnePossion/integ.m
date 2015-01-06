function sum = integ(lim,sup,index)
%定义一个求积分函数

%有两个基函数：定义一个标记index 传递这个标记值进来
t = [-0.8611363,-0.3399810,0.3399810,0.8611363];
%这是右端项中间的，两次求导后的函数
g = -2;
w  = [0.3478548,0.6521452,0.6521452,0.3478548];
x = zeros(4,1);
f=zeros(4,1);
sum = 0;
if index == 1
    for i=1:4
        x(i) = ((sup-lim)/2)*t(i) + (sup + lim)/2;
        f(i) = g*((sup-x(i))/(sup-lim));
        sum = sum + w(i)*f(i);
    end    
    sum=(sup-lim)/2*sum;
else
     for i=1:4
        x(i) = ((sup-lim)/2)*t(i) + (sup + lim)/2;
        f(i) = g*((x(i)-lim)/(sup-lim));
        sum = sum + w(i)*f(i);  
     end
     sum=(sup-lim)/2*sum;
end  
format rat;
end
