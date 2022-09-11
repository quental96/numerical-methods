%% Parametrisation
mu=3e6;
f=3e8;
c=3e8;
Dx=0.008;
Dz=Dx;
Dt=Dx/c;

x=0:Dx:1;
Nx=length(x);

z=0:Dz:0.2;
Nz=length(z);

t=0:Dt:2*Nx*Dt;
Nt=length(t);


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
G=repmat(G,1,1,Nt);

A_size=(Nx-2)*(Nz-2);
A=spalloc(A_size,A_size, 5*A_size );
b_bords=zeros(A_size,1);

%% Resolution instant par instant
handle=waitbar(0,'Solving...');
for k=2:Nt
    
% Remplissage de la matrice A et du vecteur b_bords provenant des conditions aux bords
n=0; %Numero de la ligne
for j=2:Nx-1
    for i=2:Nz-1
        n=n+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Terme diagonal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A(n,n)=1+f*Dt+4*mu*Dt/(Dx^2);
        b_bords(n)=G(i,j,k-1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Autres termes + Conditions aux bords %%%%%%%%%%%%%%%%%%%%%%%%%%
        if j==2 %%% Bord gauche
            A(n,n+(Nz-2))=-mu*Dt/(Dx^2);
            b_bords(n)=b_bords(n)+mu*Dt*G(i,j-1,k)/(Dx^2);
        elseif j==Nx-1 %%% Bord droit
            A(n,n-(Nz-2))=-mu*Dt/(Dx^2);
            b_bords(n)=b_bords(n)+mu*Dt*G(i,j+1,k)/(Dx^2);
        else
            A(n,n-(Nz-2))=-mu*Dt/(Dx^2);
            A(n,n+(Nz-2))=-mu*Dt/(Dx^2);
        end
        
        if i==2 %%% Bord haut
            A(n,n+1)=-mu*Dt/(Dx^2);
            b_bords(n)=b_bords(n)+mu*Dt*G(i-1,j,k)/(Dx^2);
        elseif i==Nz-1 %%% Bord bas
            A(n,n-1)=-mu*Dt/(Dx^2);
            b_bords(n)=b_bords(n)+mu*Dt*G(i+1,j,k)/(Dx^2);
        else
            A(n,n-1)=-mu*Dt/(Dx^2);
            A(n,n+1)=-mu*Dt/(Dx^2);
        end
        
    end
end

%% Solver
G(2:(Nz-1),2:(Nx-1),k)= reshape(A\b_bords,Nz-2,Nx-2) ;

    
    % Display current Gain
    imagesc(x,z,flipud(G(:,:,k)));colormap(hot);set(gca,'YDir','normal');axis image;
    pause(0.01);

    waitbar((k-1)/Nt);
end
close(handle);

%% Movie
handle=waitbar(0,'Making movie...');
g=round(2000/Nx);
FILM=zeros(Nz*g,Nx*g,Nt);
for i_t=1:Nt
    FILM(:,:,i_t)=imresize(G(:,:,i_t),g,'nearest');
    waitbar(i_t/Nt);
end
FILM=FILM/max(FILM(:));
implay(immovie(permute(max(1,FILM*256),[1 2 4 3]),hot(256)),min(100,Nt/10));
scrsz = get(0,'ScreenSize');
scrsz(end)=scrsz(end)-110;
set(findall(0,'tag','spcui_scope_framework'),'position',scrsz);
close(handle);