load('data.txt')

plot(data(:,1), data(:,2), 'g')
hold on


B=data(:,2);
A=ones(51,2);
for i=1:51
   A(i,2)=data(i,1);
end

X=A\B;

plot(data(:,1),A*X,'r')

clear