function C = moindrescarres(D)

load('data.txt')
X=data(:,1);
Y=data(:,2);
plot(X,Y,'g')
hold on

N=length(X);
A=zeros(N,D+1);
for i=1:N
    for j=1:D+1
        A(i,j)=X(i)^(j-1);
    end        
end

C=A\Y;

plot(X,A*C,'r')

end 
