function [deltau,uncons]=Qphild(E,F,CC,d)
[n1,m1]=size(CC);
deltau=-E\F; %global solution without constraints
uncons=deltau;
kk=0;
for i=1:n1
    if(CC(i,:)*deltau>d(i))
        kk=kk+1;
    else
        kk=kk+0;
    end
end
if (kk==0)
    return;
end
%dual quadratic programming matrices
T=CC*(E\CC');
K=(CC*(E\F)+d);
[n,m]=size(K);
lambda=zeros(n,m);
for km= 1:15
    lambda_p=lambda; %previous lambda
    for i=1:n
        w=T(i,:)*lambda-T(i,i)*lambda(i,1);
        w=w+K(i,1);
        la=-w/T(i,i);
        lambda(i,1)=max(la);
    end
    al=(lambda-lambda_p)'*(lambda-lambda_p);
    if(al<10e-5)
        break;
    end
end

deltau=-E\F - E\CC'*lambda;