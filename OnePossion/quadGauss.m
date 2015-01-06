function [quad1 quad2]= quadGauss(node1,node2)
%想办法把两个被积函数定义在这个函数里面它们分别接收传递过来的两个参数
%怎么定义这个被积函数呢
%g = (1 + (x-node1)/h)
%g = (1 + (x-node2)/h) 类似做法 基函数二也是这么定义
%可以在里面定义一个内部函数，通过传递进来的参数，分别求出这个函数表达式。

[f1,f2] = solvequad(node1,node2);
na = (node1-node2)/2;
nb = (node1+node2)/2;
quad1 = na*(0.2369269*subs(sym(f1),findsym(sym(f1)),na*0.9061793 + nb)+...
            0.2369269*subs(sym(f1),findsym(sym(f1)),-na*0.9061793 + nb)+...
            0.4786287*subs(sym(f1),findsym(sym(f1)),na*0.5384693 + nb)+...
             0.4786287*subs(sym(f1),findsym(sym(f1)),-na*0.5384693 + nb)+...
             0.5688889*subs(sym(f1),findsym(sym(f1)), nb));
quad2 = na*(0.2369269*subs(sym(f2),findsym(sym(f2)),na*0.9061793 + nb)+...
            0.2369269*subs(sym(f2),findsym(sym(f2)),-na*0.9061793 + nb)+...
            0.4786287*subs(sym(f2),findsym(sym(f2)),na*0.5384693 + nb)+...
             0.4786287*subs(sym(f2),findsym(sym(f2)),-na*0.5384693 + nb)+...
             0.5688889*subs(sym(f2),findsym(sym(f2)), nb));

function [f1,f2] = solvequad(node1,node2)
%这是利用上面的传进来的节点来计算函数的表达式
node=[node1,node2];
b=1;
a=0;
N=3;
h = (b-a)/N;
f1 = (1 - (x-node)/h);
f2 = (1 + (x-node)/h);
end

end