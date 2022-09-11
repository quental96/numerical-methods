%% Parametrisation
mu=3e6;
f=3e8;
c=3e8;
Dx=0.008;
Dz=Dx;
Dt=4e-12;

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

%% Resolution instant par instant
handle=waitbar(0,'Solving...');
for k=2:Nt-1

    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% G resolution %%%%%%%%%%%%%%%
    
    %%% TODO %%% G(2:(Nz-1),2:(Nx-1),k+1) = ??? ;
    
    for i=2:Nz-1
        for j=2:Nx-1
            G(i,j,k+1)=G(i,j,k)+Dt*mu*(G(i+1,j,k)+G(i-1,j,k)-4*G(i,j,k)+G(i,j+1,k)+G(i,j-1,k))/(Dx^2)+f*Dt*G(i,j,k);
        end
    end
    
    % Example just to check
    %G(2:(Nz-1),2:(Nx-1),k+1) = exp(-((sqrt(repmat(x(2:Nx-1),Nz-2,1).^2+repmat(flipud(z(2:Nz-1)'),1,Nx-2).^2)-k*Dx/2)/(10*Dx)).^2);
    
    % Display current Gain
    imagesc(x,z,flipud(G(:,:,k+1)));colormap(hot);set(gca,'YDir','normal');axis image;
    pause(0.01);

    waitbar(k/Nt);
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