function P=cgl_vand(order,x)
x = reshape(x,[numel(x),1]);
P=zeros(numel(x),order+1);


P(:,1)=1;    P(:,2)=x;

for k=2:order
    P(:,k+1) = 2*x.*P(:,k)-P(:,k-1);
end
