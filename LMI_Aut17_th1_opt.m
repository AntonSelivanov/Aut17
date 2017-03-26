function [P1val,P2val,P3val,betauval,betawval]=LMI_Aut17_th1_opt(D0,beta,A,B,K,muT,muB,Neum,l,nu,epsilon,alpha,h,rho)
% This MATLAB program checks the feasibility of LMIs (and optimizes C_inf) from Theorem 1  of the paper 
% A. Selivanov and E. Fridman, "Sampled-data relay control of diffusion PDEs," Automatica, 2017. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% D0       - diagonal matrix of diffusion coefficients from Assumption 1 
% beta,A,B - the parameters of the plant (3) 
% K        - nominal controller gain from Assumption 5 
% muT, muB - sector size from Assumption 4 
% Neum     =0 for Dirichlet BC, =1 for Neumann BC 
% l        - subdomain diameter 
% nu=1e-5  - tuning parameter 
% epsilon  - smoothing parameter from (4) 
% alpha    - decay rate 
% h        - sampling
% rho      - parameter from Assumption 3 (required to minimize C_inf)

% Output: 
% P1val,P2val,P3val,betauval,betawval - values of decision matrices, empty if the LMIs are not feasible

[M,L]=size(B); 
N=size(beta,2)/M; 

%% Decision variables 
P1=diag(sdpvar(M,1)); 
P2=diag(sdpvar(M,1)); 
P3=diag(sdpvar(M,1)); 
W=sdpvar(M); 
betau=sdpvar(L); 
betaw=sdpvar(L); 
Lambdaf=diag(sdpvar(M,1)); 
Lambdakappa=diag(sdpvar(M,1)); 
if Neum
    LambdaD=zeros(M); 
else
    LambdaD=diag(sdpvar(M,1)); 
end

%% LMIs
Phi=blkvar; 
Phi(1,1)=P1*(A-B*K)+(A-B*K)'*P1+2*alpha*P1-muT*muB*Lambdaf+2*N*epsilon*(1+nu^(-1))*Lambdakappa-N*pi^2*LambdaD+h*(P2*A+(P2*A)');
Phi(1,2)=(P1+h*P2)*beta; 
Phi(1,3)=P1+h*P2+1/2*(muT+muB)*Lambdaf; 
Phi(1,4)=P1*B*K; 
Phi(1,5)=h*(A'*P3-P2); 
Phi(1,6)=h*(P1*B*K)'; 
Phi(1,7)=h*P2*B; 
Phi(1,8)=h*P2*B; 
Phi(2,2)=2*kron((alpha*h*P3-P1-h*P2)*D0,eye(N))+kron(LambdaD,eye(N))+(1+nu)*l^2/pi^2*kron(Lambdakappa,eye(N)); 
Phi(2,5)=h*(P3*beta)'; 
Phi(3,3)=-Lambdaf; 
Phi(3,5)=h*P3; 
Phi(4,4)=-Lambdakappa; 
Phi(4,6)=-h*(P1*B*K)'; 
Phi(5,5)=h*(exp(2*alpha*h)*W-2*P3); 
Phi(5,7)=h*P3*B; 
Phi(5,8)=h*P3*B; 
Phi(6,6)=-pi^2*h/4*W; 
Phi(6,7)=h*P1*B; 
Phi(6,8)=h*P1*B; 
Phi(7,7)=-h*betau; 
Phi(8,8)=-h*betaw; 
Phi=sdpvar(Phi); 

%% Solution of LMIs
LMIs=[Phi<=0,W>=0,betau>=0,betaw>=0,P1>=eye(M),P3>=0,Lambdaf>=0,Lambdakappa>=0,LambdaD>=0]; % eye(M) is a scaling to prevent P1 -> 0 during optimization

options=sdpsettings('solver','lmilab','verbose',0); 
sol=optimize(LMIs,norm(betau)+norm(betaw)*rho^2,options); % minimizing C_inf 

P1val=[]; P2val=[]; P3val=[]; betauval=[]; betawval=[]; 
if sol.problem==0
    [primal,~]=check(LMIs); 
    if min(primal)>=0 
        P1val=value(P1); 
        P2val=value(P2); 
        P3val=value(P3); 
        betauval=value(betau); 
        betawval=value(betaw); 
    end
else
    yalmiperror(sol.problem) 
end