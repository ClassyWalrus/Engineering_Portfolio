nodes=20;
xstep=4/(nodes-1);
x=zeros(nodes,1);
for i=1:nodes
    if i==1
        x(1)=0;
    else
        x(i)=x(i-1)+xstep;
    end
end
x
y=sqrt(abs(4-(x-2).^2))
node=[x, y]