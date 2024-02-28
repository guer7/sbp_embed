function dP=cgl_vand_der(order,x,P)

dP=zeros(numel(x),order+1);


dP(:,1)=0;    dP(:,2)=1;

for k=2:order
    dP(:,k+1) = 2*x.*dP(:,k) + 2*P(:,k) - dP(:,k-1);
end
