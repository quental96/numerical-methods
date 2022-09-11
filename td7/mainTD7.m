mu=3e6;
f=3e8;
Dx=0.008;
Dy=Dx;

x=0:Dx:1;
Nx=length(x)-1;

y=0:Dy:0.2;
Ny=length(y)-1;

Dt=Dx/mu;
t=0:Dt:2*Nx*Dt;
Nt=length(t)-1;

G=zeros(Nx,Ny,Nt);

if i==1 || i==Nx
                G(i,j,k+1)=0;
            end
            
            if j==1 || j==Ny
                G(i,j,k+1)=sin(3*Dx*(i-1)*pi)^2;
            end

for k=1:Nt
    G(:,:,k+1)=zeros(Nx,Ny);
    for i=1:Nx
        for j=1:Ny
            
            G(i,j,k+1)=G(i,j,k)+Dt*(-f*G(i,j,k)+mu/(Dx^2)*(G(i-1,j,k)+G(i+1,j,k)+G(i,j-1,k)+G(i,j+1,k))-4*mu/(Dx^2)*G(i,j,k));
        end
    end
end

NA=Nx*Ny;
A=spalloc(NA,NA,5*NA);



b=zeros(NA,1);