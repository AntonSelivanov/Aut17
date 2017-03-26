% This MATLAB program checks the feasibility of LMIs from Theorem 1 of the paper 
% A. Selivanov and E. Fridman, "Sampled-data relay control of diffusion PDEs," Automatica, 2017. 

%% Example 1: 2D Catalytic Slab 
% Plant parameters 
betaT=50;           % heat of reaction
betaU=2;            % heat transfer coefficient 

% System parameters 
D0=1/sqrt(2)/pi^2;  % diffusion
beta=0;             % convection
A=-betaU;           % reaction
B=betaU;            % control matrix
K=4;                % nominal controller gain from Assumption 5
muT=6.15;           % sector size from Assumption 4 
muB=0;
Neum=0;             % Dirichlet boundary conditions 
Ns=36;              % number of (equal) subdomains
l=sqrt(2/Ns);       % subdomain diameter
nu=1e-5;            % tuning parameter
epsilon=1e-9;       % smoothing parameter from (4) 
alpha=2.4;          % decay rate 
h=1.4e-3;           % sampling 
rho=.01;            % parameter from Assumption 3 
ai=[-.1 .1];        % dual vectors from (13) 

[P1,~,~,betau,betaw]=LMI_Aut17_th1_opt(D0,beta,A,B,K,muT,muB,Neum,l,nu,epsilon,alpha,h,rho); 

if isempty(P1) 
    disp('Example 1: not feasible'); 
else
    C0=min(diag(ai'*K*P1^(-1)*K'*ai).^(-1))*(l^2/2); 
    Cinf=h/(2*alpha)*(eigs(betau,1)+rho^2*eigs(betaw,1))*max((1./ai).^2); 
    if Cinf<(1-rho)^2*C0
        disp('Example 1: '); 
        disp(['    C0=' num2str(C0)]); 
        disp(['    Cinf=' num2str(Cinf)]); 
    else
       disp('Example 1: ');  
       disp('  Cinf >= (1-rho)^2*C0'); % (14) is violated 
    end    
end

%% Example 2: Chemical Reactor
% Plant parameters 
Le=100;                     % Lewis number
V=1.1;                      % convective velocity
D=10;                       % diffusion coefficient 
betaPar=.45;                % parameters of g(z)
d=.2; 

% System parameters 
D0=diag([1/Le, D]);         % diffusion
beta=diag([-V/Le, -V]);     % convection
A=[0 1/Le; -betaPar -d];    % reaction
B=[1; 0];                   % control matrix
K=[2 0];                    % nominal controller gain from Assumption 5
muT=[1 0; 0 0];            % sector size from Assumption 4 
muB=zeros(2);
Neum=1;                     % Dirichlet boundary conditions 
Ns=4;                       % number of (equal) subdomains
l=1/Ns;                      % subdomain diameter
nu=1e-5;                    % tuning parameter
epsilon=1e-7;               % smoothing parameter from (4) 
alpha=.14;                  % decay rate 
h=1e-3;                     % sampling 
rho=.01;                    % parameter from Assumption 3 
ai=[-.5 .5];                % dual vectors from (13) 

[P1,~,~,betau,betaw]=LMI_Aut17_th1(D0,beta,A,B,K,muT,muB,Neum,l,nu,epsilon,alpha,h); 

if isempty(P1) 
    disp('Example 2: not feasible'); 
else 
    C0=min(diag(ai'*K*P1^(-1)*K'*ai).^(-1))*l; 
    Cinf=h/(2*alpha)*(eigs(betau,1)+rho^2*eigs(betaw,1))*max((1./ai).^2); 
    if Cinf<(1-rho)^2*C0 
        disp('Example 2:'); 
        disp(['    C0=' num2str(C0)]); 
        disp(['    Cinf=' num2str(Cinf)]); 
    else 
       disp('Example 2:');  
       disp('  Cinf >= (1-rho)^2*C0'); % (14) is violated 
    end 
end