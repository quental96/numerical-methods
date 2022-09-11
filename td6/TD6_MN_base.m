%% Parametrisation
mu=3e6;
f=3e8;
Dx=0.008;
Dz=Dx;

x=0:Dx:1;
Nx=length(x);

z=0:Dz:0.2;
Nz=length(z);


% Initialisation
% Gain
G=zeros(Nz,Nx);
%%% TO DO %%%   G0=@(u) ? ;
%%% TO DO %%%   G1=@(u) ? ;
%%% TO DO %%%   G2=@(u) ? ;
%%% TO DO %%%   G3=@(u) ? ;
G(:,1)  =G0(z);
G(:,end)=G1(z);
G(1,:)  =G2(x);
G(end,:)=G3(x);


% Calculs
% Allocations
%%% TO DO %%%   A_size= ? ;
%%% TO DO %%%   A=spalloc(A_size,A_size, ? );
b_bords=zeros(A_size,1);

%%% TODO: remplir la matrice A et le vecteur b_bords provenant des conditions aux bords

% Solver
%%% TO DO %%%   G(2:(Nz-1),2:(Nx-1))= ? ;

% Drawing
imagesc(x,z,flipud(G));colormap(hot);set(gca,'YDir','normal');








