%% TD4 correction

clear all
close all
clc

% Constants
mass=1.44e-25;
kb=1.38e-23;

load('Data_temperature1.mat');
ydata=data1;


%% Partie 1
C=zeros(length(xdata),3);
for k=1:3
    C(:,k)=xdata.^(k-1);
end

% convert polynomial coefficients to A, mu, sigma^2:
extract = @(p) [...
    exp(p(1)+1/2*p(2)*-1/(2*p(3))*p(2)), ... % A
    -1/(2*p(3))*p(2), ... % mu
    -1/(2*p(3)) ... % sigma^2
    ];

% eval Gaussian function
gauss_func = @(param,x) param(1)*exp(-(x-param(2)).^2/(2*param(3)));

% linear least-square fit:
p = C\log(ydata);

% eval, compute residual, plot
ydatafit=gauss_func(extract(p),xdata);
fprintf('residu: %d\n',norm(ydatafit-ydata)^2);
figure
plot(xdata,ydata,xdata,ydatafit)

%% Partie 2.1

% extract extremum
[v,i] = max(ydata);
% hard constraint: [1 x_i x_i^2] * p = log(v)
A = [C'*C    C(i,:)' ; ...
     C(i,:)  0       ];
b = [C'*log(ydata) ; ...
     log(v)        ];
p = A \ b;

% eval, compute residual, plot
ydatafit=gauss_func(extract(p),xdata);
fprintf('residu: %d\n',norm(ydatafit-ydata)^2);
figure
plot(xdata,ydata,xdata,ydatafit)

%% Partie 2.2

h = 73;

l = length(xdata);
A=zeros(2*l+1,3);
b = zeros(2*l+1,1);
A(1:l,    :) =  C; b(1:l)     =  log(ydata+h);           % g(x_i) < y_i+h
A(1+l:2*l,:) = -C; b(1+l:2*l) = -log(max(1e-3,ydata-h)); % y_i-h  < g(x_i)
A(end,:) = [0 0 1];                                      % p(2) < 0

p = lsqlin(C,log(ydata),A,b);

% eval, compute residual, plot
ydatafit=gauss_func(extract(p),xdata);
fprintf('residu: %d\n',norm(ydatafit-ydata)^2);
figure
plot(xdata,ydata,xdata,ydatafit)

%% Partie 2.3

figure
plot(xdata,ydata)
hold on
e=[];
for j=0:7
    wn = ydata.^j;
    ptemp = (C'*diag(wn)*C) \ (C'*diag(wn)*(log(ydata)));
    ydatafit=gauss_func(extract(ptemp),xdata);    
    plot(xdata,ydatafit)  
    hold on;
    e(j+1)=[norm(ydatafit-ydata)^2];
end
hold off;

[v j] = min(e);

best_exponent = j-1;
fprintf('best residu: %d for j=%d\n',v,best_exponent);

% recompute best solution
wn = ydata.^j;
p = (C'*diag(wn)*C) \ (C'*diag(wn)*(log(ydata)));
ydatafit=gauss_func(extract(p),xdata);


%% Comparison to non-linear fit

figure
plot(xdata,ydata,xdata,ydatafit)
hold on;

load('data_temperature_non_lin1.mat')
plot(xdata,data_non_lin1,'linewidth',2)

fprintf('non-linear residu: %d\n',norm(data_non_lin1-ydata)^2 );

%% Partie 3
load('Data_temperature_tot.mat')

C=zeros(length(xdata),3);
for k=1:3
    C(:,k)=xdata.^(k-1);
end
    
tof=[6 8 10 12 14 16 18]*1e-3;

% Getting sigma^2
A=zeros(1,3);
A(1,3)=1;
b=[0];

figure
for j=1:7
    ydata=data(:,j);
    wn=ydata.^best_exponent;
    ptemp=lsqlin(diag(sqrt(wn))*C,log((ydata)).*sqrt(wn),A,b);
    gp = extract(ptemp);
    s2(j)=gp(3);
    ydatafit=gauss_func(gp,xdata);
    plot(xdata,ydata,xdata,ydatafit)
    hold on;
end
figure
plot(tof,s2)

% Parabolic fit of sigma^2

Cb=zeros(length(s2),2);
Cb(:,1)=1;
Cb(:,2)=tof.^2;
Ab=[0 -1];
b=[0];
pb=lsqlin(Cb,s2,Ab,b);

fprintf('Temperature = %d\n',mass*pb(2)/kb)


