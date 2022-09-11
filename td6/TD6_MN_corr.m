%% Parametrisation
mu=3e6;
f=3e8;
Dx=0.008;
Dz=Dx;

x=0:Dx:1;
Nx=length(x);

z=0:Dz:0.2;
Nz=length(z);


%% Initialisation
% Gain
G=zeros(Nz,Nx);
G0=@(u) 0 ;
G1=@(u) 0 ;
G2=@(u) sin(3*pi*u).^2 ;
G3=@(u) sin(3*pi*u).^2 ;
G(:,1)  =G0(z);
G(:,end)=G1(z);
G(1,:)  =G2(x);
G(end,:)=G3(x);


%% Calculs
% Allocations
A_size=(Nx-2)*(Nz-2);
A=spalloc(A_size,A_size, 5*A_size );
b_bords=zeros(A_size,1);

% Remplissage de la matrice A et du vecteur b_bords provenant des conditions aux bords
n=0; %Numero de la ligne
for j=2:Nx-1
    for i=2:Nz-1
        n=n+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Terme diagonal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A(n,n)=-4-f/mu*Dx^2;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Autres termes + Conditions aux bords %%%%%%%%%%%%%%%%%%%%%%%%%%
        if j==2 %%% Bord gauche
            A(n,n+(Nz-2))=1;
            b_bords(n)=b_bords(n)-G(i,j-1);
        elseif j==Nx-1 %%% Bord droit
            A(n,n-(Nz-2))=1;
            b_bords(n)=b_bords(n)-G(i,j+1);
        else
            A(n,n-(Nz-2))=1;
            A(n,n+(Nz-2))=1;
        end
        
        if i==2 %%% Bord haut
            A(n,n+1)=1;
            b_bords(n)=b_bords(n)-G(i-1,j);
        elseif i==Nz-1 %%% Bord bas
            A(n,n-1)=1;
            b_bords(n)=b_bords(n)-G(i+1,j);
        else
            A(n,n-1)=1;
            A(n,n+1)=1;
        end
        
    end
end

%% Solver
G(2:(Nz-1),2:(Nx-1))= reshape(A\b_bords,Nz-2,Nx-2) ;

%% Drawing
imagesc(x,z,flipud(G));colormap(hot);set(gca,'YDir','normal');axis image;








