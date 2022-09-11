function C = moindrescarres2(D)

load('data4.txt')
X=data4(:,1);
Y=data4(:,2);
plot(X,Y,'g')
hold on

N=length(X);
A=zeros(N,D+1);
for i=1:N
    for j=1:D+1
        A(i,j)=sin((j-1)*X(i));
    end        
end

C=A\Y;

plot(X,A*C,'r')


end