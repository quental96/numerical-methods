mu=3e6;
f=3e8;
Dx=0.008;
Dy=Dx;

x=0:Dx:1;
Nx=length(x)-1;

y=0:Dy:0.2;
Ny=length(y)-1;

NA=Nx*Ny;
A=spalloc(NA,NA,5*NA);

b=zeros(NA,1);

for i=1:Nx
    for j=1:Ny
        
        k=(i-1)*Ny+j;
        
        A(k,k)=-4-(f/mu)*Dx^2;
        
        if i==1 || i==Nx
            A(k,k)=1; 
            b(k)=0;
        elseif j==1 || j==Ny
            A(k,k)=1;
            b(k)=sin(3*Dx*(i-1)*pi)^2;
        else
            A(k,k+1)=1;
            A(k,k-1)=1;
            A(k,k+Ny)=1;
            A(k,k-Ny)=1;
        end
      
    end
end

g=A\b;

surf(reshape(g,[Ny,Nx]))
