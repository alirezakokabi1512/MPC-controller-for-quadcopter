function [CC,dd,dupast] = constraints_mimo(Dumax,u_max,u_min,n_inputs,nu)
CC=[];
dim = nu*n_inputs; %matrices dimensions
CC(1:dim,1:dim)=eye(dim);
CC(dim+1:2*dim,1:dim) = -eye(dim);
s=0;
E=zeros(dim,dim);
for j=1:nu
    E(1+s:n_inputs*j,1:n_inputs) = eye(n_inputs);
    s=s+n_inputs;
end
p=n_inputs+1;
for counter = 1:nu-1
    i=counter;
    for k = (n_inputs*counter)+1:n_inputs:dim-n_inputs+1
        E(k:n_inputs*(i+1),p:n_inputs*(counter+1))= eye(n_inputs);
        i=i+1;
    end
    p=p+n_inputs;
end

Cu=[eye(dim);-eye(dim)];
CuE=Cu*E;
CC(2*dim+1:4*dim,1:dim)=CuE;
ddeltaU=zeros(2*dim,1);
t=1;
for i=1:n_inputs:size(ddeltaU,1)
    ddeltaU(i:t*n_inputs,1) = Dumax;
    t=t+1;
end
du=zeros(2*dim,1);
%upper limits
t=1;
for i=1:n_inputs:size(du,1)/2
    du(i:t*n_inputs,1)=u_max;
    t=t+1;
end
% lower limits
for i=1+size(du,1)/2 : n_inputs: size(du,1)
    du(i:t*n_inputs,1)=-u_min;
    t=t+1;
end

dd=[ddeltaU;du];
% dupast past inputs

L=zeros();
t=1;
for i=1:n_inputs:(n_inputs*nu)
    L(i:n_inputs*t,1:n_inputs)=eye(n_inputs);
    t=t+1;
end

CuL=Cu*L;
dupast=[zeros(size(ddeltaU,1),n_inputs);-CuL];






