function C = mcPonderes(D)

load('data4.txt')
X=data4(:,1);
Y=data4(:,2);
plot(X,Y,'g')
hold on

N=length(X);
A=zeros(N,D+1);

for i=1:N
    for j=1:D+1
        A(i,j)=X(i)^(j-1);
    end        
end

m=mean(X);
h=10;
w=zeros(N,1);
for i=1:N
    if abs(X(i)-m) < h
        w(i)=(1-((X(i)-m)/h)^2)^4;
    else
        w(i)=0
    end
end

W=diag(sqrt(w));

Ap=W*A;
Yp=W*Y;
C=Ap\Yp;

plot(X,A*C,'r')

end