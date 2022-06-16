function [F_sm,F_sn,Z_sm,Z_sn,P_sm,P_sn,P,Es_sm,Es_sn,obj,err] = CDLRR(Tm,Tn,Ssm,Ssn,a,b,opts)

% Solve Multi-site data heterogeneity problem by CDLRR(category-discrepancy low rank representation)
% Slove the objective function of CDLRR minimization problem by ADMM
% min_{P,Pim,Zim,Eim}sigma_{i=1:K}(sigma_{m=1:N}(||Zim||_*+α||Eim||_1+sigma_{n=1:N,n≠m}((γ/2)*||(PinSin)^T(PimSim)||_F~2)))
% s.t. PimSim=PTmZim+Eim,PimPim'=I,PinPin'=I,i=1,...K
% Introduce one relaxation variables Fim, Zim=Fim
% Augmented Lagrance Multiplierfunction
%           L = sigma_{i=1:K}(sigma_{m=1:N}(||Fim||_*+α||Eim||_1
%               +sigam_{n=1:N,n≠m}((γ/2)*||(PinSin)^T(PimSim)||_F~2)          
%               +<Y1,im,Zim-Fim>+<Y2,im,PimSim-PTmZim-Eim>
%               +(mu/2)*(||Zim-Fim||_F~2+||PimSim-PTmZim-Eim||_F~2)))
% -------------------------------------------------------------------------------------------------------------
% Example: Seven sites' data heterogeneity processing; Two categories per center
% Input: 
%       Tm     --     d*nTm matrix     --     Target domain class m samples
%       Tn     --     d*nTn matrix     --     Target domain class n samples
%       Ssm    --     Structure value in Matlab. The fields are
%            Ssm.Sim    - the i-th source domain class m samples - d*nSim matrix
%       Ssn    --     Structure value in Matlab. The fields are
%            Ssn.Sin    - the i-th source domain class n samples - d*nSin matrix  
%       α      --     >0,balance parameter Replace αwith a in the code
%       γ      --     >0,balance parameter Replace γwith b in the code
%       opts   --     Structure value in Matlab. The fields are
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.DEBUG      -   0 or 1
%
% Output:
%       P     --     d*d matrix     --     Target domain common transformer matrix
%       P_sm    --     Structure value in Matlab. The fields are
%           P_sm.P1m     - the i-th source domain class m samples transformer matrix  - d*d matrix
%       P_sn    --     Structure value in Matlab. The fields are
%           P_sn.P1n     - the i-th source domain class n samples transformer matrix  - d*d matrix
%       Z_sm    --     Structure value in Matlab. The fields are
%           Z_sm.Z_im     - the i-th source domain class m samples low-rank representation matrix  - nTm*nSim matrix
%       Z_sn    --     Structure value in Matlab. The fields are
%           Z_sn.Z_in     - the i-th source domain class n samples low-rank representation matrix  - nTn*nSin matrix
%       E_sm    --     Structure value in Matlab. The fields are
%           E_sm.E1m     - the i-th source domain class m samples sparse error matrix of data representation  - d*nSim matrix
%       E_sn    --     Structure value in Matlab. The fields are
%           E_sn.E1n     - the i-th source domain class n samples sparse error matrix of data representation  - d*nSin matrix
%       obj   --     objective function value
%       err   --     residual 
% version 3.0 - 12/11/2021
% Written by Jie Ding 


%Initialize
tol = 1e-8;
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;
DEBUG = 0;

if ~exist('opts', 'var')
    opts = [];
end

if isfield(opts, 'tol');        tol = opts.tol;            end
if isfield(opts, 'rho');        rho = opts.rho;            end
if isfield(opts, 'mu');         mu = opts.mu;              end
if isfield(opts, 'max_iter');   max_iter = opts.max_iter;  end
if isfield(opts, 'max_mu');     max_mu = opts.max_mu;      end
if isfield(opts, 'DEBUG');      DEBUG = opts.DEBUG;        end

% Target
[d,nTm] = size(Tm);
[~,nTn] = size(Tn);

% Source
S1m = Ssm.S1m;
[~,nS1m] = size(S1m);
S1n = Ssn.S1n;
[~,nS1n] = size(S1n);

S2m = Ssm.S2m;
[~,nS2m] = size(S2m);
S2n = Ssn.S2n;
[~,nS2n] = size(S2n);

S3m = Ssm.S3m;
[~,nS3m] = size(S3m);
S3n = Ssn.S3n;
[~,nS3n] = size(S3n);

S4m = Ssm.S4m;
[~,nS4m] = size(S4m);
S4n = Ssn.S4n;
[~,nS4n] = size(S4n);

S5m = Ssm.S5m;
[~,nS5m] = size(S5m);
S5n = Ssn.S5n;
[~,nS5n] = size(S5n);

S6m = Ssm.S6m;
[~,nS6m] = size(S6m);
S6n = Ssn.S6n;
[~,nS6n] = size(S6n);

I = eye(d);
mibu = 1e-7*I; % Avoiding the inverse of a matrix of insufficient rank

% 变量初始化
Z1m = zeros(nTm,nS1m);Z1n = zeros(nTn,nS1n);
Z2m = zeros(nTm,nS2m);Z2n = zeros(nTn,nS2n);
Z3m = zeros(nTm,nS3m);Z3n = zeros(nTn,nS3n);
Z4m = zeros(nTm,nS4m);Z4n = zeros(nTn,nS4n);
Z5m = zeros(nTm,nS5m);Z5n = zeros(nTn,nS5n);
Z6m = zeros(nTm,nS6m);Z6n = zeros(nTn,nS6n);
Z1m_c = Z1m;Z1n_c = Z1n;
Z2m_c = Z2m;Z2n_c = Z2n;
Z3m_c = Z3m;Z3n_c = Z3n;
Z4m_c = Z4m;Z4n_c = Z4n;
Z5m_c = Z5m;Z5n_c = Z5n;
Z6m_c = Z6m;Z6n_c = Z6n;

F1m = zeros(nTm,nS1m);F1n = zeros(nTn,nS1n);
F2m = zeros(nTm,nS2m);F2n = zeros(nTn,nS2n);
F3m = zeros(nTm,nS3m);F3n = zeros(nTn,nS3n);
F4m = zeros(nTm,nS4m);F4n = zeros(nTn,nS4n);
F5m = zeros(nTm,nS5m);F5n = zeros(nTn,nS5n);
F6m = zeros(nTm,nS6m);F6n = zeros(nTn,nS6n);

% 中间变量
F1m_c = F1m;F1n_c = F1n;
F2m_c = F2m;F2n_c = F2n;
F3m_c = F3m;F3n_c = F3n;
F4m_c = F4m;F4n_c = F4n;
F5m_c = F5m;F5n_c = F5n;
F6m_c = F6m;F6n_c = F6n;

Z_I_m = eye(nTm);
Z_I_n = eye(nTn);

E1m = zeros(d,nS1m);E1n = zeros(d,nS1n);
E2m = zeros(d,nS2m);E2n = zeros(d,nS2n);
E3m = zeros(d,nS3m);E3n = zeros(d,nS3n);
E4m = zeros(d,nS4m);E4n = zeros(d,nS4n);
E5m = zeros(d,nS5m);E5n = zeros(d,nS5n);
E6m = zeros(d,nS6m);E6n = zeros(d,nS6n);
E1m_c = E1m;E1n_c = E1n;
E2m_c = E2m;E2n_c = E2n;
E3m_c = E3m;E3n_c = E3n;
E4m_c = E4m;E4n_c = E4n;
E5m_c = E5m;E5n_c = E5n;
E6m_c = E6m;E6n_c = E6n;

P1m = I;P1n = I;
P2m = I;P2n = I;
P3m = I;P3n = I;
P4m = I;P4n = I;
P5m = I;P5n = I;
P6m = I;P6n = I;
P1m_c = P1m;P1n_c = P1n;
P2m_c = P2m;P2n_c = P2n;
P3m_c = P3m;P3n_c = P3n;
P4m_c = P4m;P4n_c = P4n;
P5m_c = P5m;P5n_c = P5n;
P6m_c = P6m;P6n_c = P6n;

P = I;

% lagrangian multiplier 拉格朗日乘子
Y1_1m = zeros(nTm,nS1m);Y1_1n = zeros(nTn,nS1n);Y2_1m = zeros(d,nS1m);Y2_1n = zeros(d,nS1n);
Y1_2m = zeros(nTm,nS2m);Y1_2n = zeros(nTn,nS2n);Y2_2m = zeros(d,nS2m);Y2_2n = zeros(d,nS2n);
Y1_3m = zeros(nTm,nS3m);Y1_3n = zeros(nTn,nS3n);Y2_3m = zeros(d,nS3m);Y2_3n = zeros(d,nS3n);
Y1_4m = zeros(nTm,nS4m);Y1_4n = zeros(nTn,nS4n);Y2_4m = zeros(d,nS4m);Y2_4n = zeros(d,nS4n);
Y1_5m = zeros(nTm,nS5m);Y1_5n = zeros(nTn,nS5n);Y2_5m = zeros(d,nS5m);Y2_5n = zeros(d,nS5n);
Y1_6m = zeros(nTm,nS6m);Y1_6n = zeros(nTn,nS6n);Y2_6m = zeros(d,nS6m);Y2_6n = zeros(d,nS6n);

for iter = 1:max_iter
    % Update Pim
        % Update Fin_c
    [F1n_c,~] = prox_nuclear(Z1n_c+Y1_1n/mu,1/mu);
    [F2n_c,~] = prox_nuclear(Z2n_c+Y1_2n/mu,1/mu);
    [F3n_c,~] = prox_nuclear(Z3n_c+Y1_3n/mu,1/mu);
    [F4n_c,~] = prox_nuclear(Z4n_c+Y1_4n/mu,1/mu);
    [F5n_c,~] = prox_nuclear(Z5n_c+Y1_5n/mu,1/mu);
    [F6n_c,~] = prox_nuclear(Z6n_c+Y1_6n/mu,1/mu);
        % Update Fim
    [F1m,nuclearnormF1m] = prox_nuclear(Z1m+Y1_1m/mu,1/mu);
    [F2m,nuclearnormF2m] = prox_nuclear(Z2m+Y1_2m/mu,1/mu);
    [F3m,nuclearnormF3m] = prox_nuclear(Z3m+Y1_3m/mu,1/mu);
    [F4m,nuclearnormF4m] = prox_nuclear(Z4m+Y1_4m/mu,1/mu);
    [F5m,nuclearnormF5m] = prox_nuclear(Z5m+Y1_5m/mu,1/mu);
    [F6m,nuclearnormF6m] = prox_nuclear(Z6m+Y1_6m/mu,1/mu);
        % Update Ein_c
    E1n_c = prox_l1(P1n_c*S1n-P*Tn*Z1n_c+Y2_1n/mu,a/mu);
    E2n_c = prox_l1(P2n_c*S2n-P*Tn*Z2n_c+Y2_2n/mu,a/mu);
    E3n_c = prox_l1(P3n_c*S3n-P*Tn*Z3n_c+Y2_3n/mu,a/mu);
    E4n_c = prox_l1(P4n_c*S4n-P*Tn*Z4n_c+Y2_4n/mu,a/mu);
    E5n_c = prox_l1(P5n_c*S5n-P*Tn*Z5n_c+Y2_5n/mu,a/mu);
    E6n_c = prox_l1(P6n_c*S6n-P*Tn*Z6n_c+Y2_6n/mu,a/mu);
        % Update Eim
    E1m = prox_l1(P1m*S1m-P*Tm*Z1m+Y2_1m/mu,a/mu);
    E2m = prox_l1(P2m*S2m-P*Tm*Z2m+Y2_2m/mu,a/mu);
    E3m = prox_l1(P3m*S3m-P*Tm*Z3m+Y2_3m/mu,a/mu);
    E4m = prox_l1(P4m*S4m-P*Tm*Z4m+Y2_4m/mu,a/mu);
    E5m = prox_l1(P5m*S5m-P*Tm*Z5m+Y2_5m/mu,a/mu);
    E6m = prox_l1(P6m*S6m-P*Tm*Z6m+Y2_6m/mu,a/mu);
        % Update Zin_c
    G2_1_c = P1n_c*S1n-E1n_c+Y2_1n/mu;
    Z1n_c = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_1_c+F1n_c-Y1_1n/mu);
    G2_2_c = P2n_c*S2n-E2n_c+Y2_2n/mu;
    Z2n_c = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_2_c+F2n_c-Y1_2n/mu);
    G2_3_c = P3n_c*S3n-E3n_c+Y2_3n/mu;
    Z3n_c = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_3_c+F3n_c-Y1_3n/mu);
    G2_4_c = P4n_c*S4n-E4n_c+Y2_4n/mu;
    Z4n_c = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_4_c+F4n_c-Y1_4n/mu);
    G2_5_c = P5n_c*S5n-E5n_c+Y2_5n/mu;
    Z5n_c = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_5_c+F5n_c-Y1_5n/mu);
    G2_6_c = P6n_c*S6n-E6n_c+Y2_6n/mu;
    Z6n_c = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_6_c+F6n_c-Y1_6n/mu);
        % Update Zim
    G1_1 = P1m*S1m - E1m + Y2_1m/mu;
    Z1m = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_1+F1m-Y1_1m/mu);
    G1_2 = P2m*S2m - E2m + Y2_2m/mu;
    Z2m = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_2+F2m-Y1_2m/mu);
    G1_3 = P3m*S3m - E3m + Y2_3m/mu;
    Z3m = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_3+F3m-Y1_3m/mu);
    G1_4 = P4m*S4m - E4m + Y2_4m/mu;
    Z4m = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_4+F4m-Y1_4m/mu);
    G1_5 = P5m*S5m - E5m + Y2_5m/mu;
    Z5m = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_5+F5m-Y1_5m/mu);
    G1_6 = P6m*S6m - E6m + Y2_6m/mu;
    Z6m = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_6+F6m-Y1_6m/mu);
        % Update Pin_c
    G5_1_c = P1m*(S1m*S1m')*P1m';G6_1_c = P*Tn*Z1n_c+E1n_c-Y2_1n/mu;
    P1n_c = ((2*b*G5_1_c+mu*I)^-1)*mu*G6_1_c*S1n'*((S1n*S1n'+mibu)^-1);
    G5_2_c = P2m*(S2m*S2m')*P2m';G6_2_c = P*Tn*Z2n_c+E2n_c-Y2_2n/mu;
    P2n_c = ((2*b*G5_2_c+mu*I)^-1)*mu*G6_2_c*S2n'*((S2n*S2n'+mibu)^-1);
    G5_3_c = P3m*(S3m*S3m')*P3m';G6_3_c = P*Tn*Z3n_c+E3n_c-Y2_3n/mu;
    P3n_c = ((2*b*G5_3_c+mu*I)^-1)*mu*G6_3_c*S3n'*((S3n*S3n'+mibu)^-1);
    G5_4_c = P4m*(S4m*S4m')*P4m';G6_4_c = P*Tn*Z4n_c+E4n_c-Y2_4n/mu;
    P4n_c = ((2*b*G5_4_c+mu*I)^-1)*mu*G6_4_c*S4n'*((S4n*S4n'+mibu)^-1);
    G5_5_c = P5m*(S5m*S5m')*P5m';G6_5_c = P*Tn*Z5n_c+E5n_c-Y2_5n/mu;
    P5n_c = ((2*b*G5_5_c+mu*I)^-1)*mu*G6_5_c*S5n'*((S5n*S5n'+mibu)^-1);
    G5_6_c = P6m*(S6m*S6m')*P6m';G6_6_c = P*Tn*Z6n_c+E6n_c-Y2_6n/mu;
    P6n_c = ((2*b*G5_6_c+mu*I)^-1)*mu*G6_6_c*S6n'*((S6n*S6n'+mibu)^-1);  
    % Update Pim
    G3_1 = P1n_c*(S1n*S1n')*P1n_c';G4_1 = P*Tm*Z1m+E1m-Y2_1m/mu;
    P1m = ((2*b*G3_1+mu*I)^-1)*mu*G4_1*S1m'*((S1m*S1m'+mibu)^-1);
    G3_2 = P2n_c*(S2n*S2n')*P2n_c';G4_2 = P*Tm*Z2m+E2m-Y2_2m/mu;
    P2m = ((2*b*G3_2+mu*I)^-1)*mu*G4_2*S2m'*((S2m*S2m'+mibu)^-1);
    G3_3 = P3n_c*(S3n*S3n')*P3n_c';G4_3 = P*Tm*Z3m+E3m-Y2_3m/mu;
    P3m = ((2*b*G3_3+mu*I)^-1)*mu*G4_3*S3m'*((S3m*S3m'+mibu)^-1);
    G3_4 = P4n_c*(S4n*S4n')*P4n_c';G4_4 = P*Tm*Z4m+E4m-Y2_4m/mu;
    P4m = ((2*b*G3_4+mu*I)^-1)*mu*G4_4*S4m'*((S4m*S4m'+mibu)^-1);
    G3_5 = P5n_c*(S5n*S5n')*P5n_c';G4_5 = P*Tm*Z5m+E5m-Y2_5m/mu;
    P5m = ((2*b*G3_5+mu*I)^-1)*mu*G4_5*S5m'*((S5m*S5m'+mibu)^-1);
    G3_6 = P6n_c*(S6n*S6n')*P6n_c';G4_6 = P*Tm*Z6m+E6m-Y2_6m/mu;
    P6m = ((2*b*G3_6+mu*I)^-1)*mu*G4_6*S6m'*((S6m*S6m'+mibu)^-1);
    
    % Pim s.t. Pim'Pim=I Orithogonal constraint
    if cast(P1m'*P1m, 'int8') ~= eye(size(P1m)); P1m = orth(P1m); end
    if cast(P2m'*P2m, 'int8') ~= eye(size(P2m)); P2m = orth(P2m); end
    if cast(P3m'*P3m, 'int8') ~= eye(size(P3m)); P3m = orth(P3m); end
    if cast(P4m'*P4m, 'int8') ~= eye(size(P4m)); P4m = orth(P4m); end
    if cast(P5m'*P5m, 'int8') ~= eye(size(P5m)); P5m = orth(P5m); end
    if cast(P6m'*P6m, 'int8') ~= eye(size(P6m)); P6m = orth(P6m); end
    
    % Update Pin
        % Update Fim_c
    [F1m_c,~] = prox_nuclear(Z1m_c+Y1_1m/mu,1/mu);
    [F2m_c,~] = prox_nuclear(Z2m_c+Y1_2m/mu,1/mu);
    [F3m_c,~] = prox_nuclear(Z3m_c+Y1_3m/mu,1/mu);
    [F4m_c,~] = prox_nuclear(Z4m_c+Y1_4m/mu,1/mu);
    [F5m_c,~] = prox_nuclear(Z5m_c+Y1_5m/mu,1/mu);
    [F6m_c,~] = prox_nuclear(Z6m_c+Y1_6m/mu,1/mu);
        % Update Fin
    [F1n,nuclearnormF1n] = prox_nuclear(Z1n+Y1_1n/mu,1/mu);
    [F2n,nuclearnormF2n] = prox_nuclear(Z2n+Y1_2n/mu,1/mu);
    [F3n,nuclearnormF3n] = prox_nuclear(Z3n+Y1_3n/mu,1/mu);
    [F4n,nuclearnormF4n] = prox_nuclear(Z4n+Y1_4n/mu,1/mu);
    [F5n,nuclearnormF5n] = prox_nuclear(Z5n+Y1_5n/mu,1/mu);
    [F6n,nuclearnormF6n] = prox_nuclear(Z6n+Y1_6n/mu,1/mu);
        % Update Eim_c
    E1m_c = prox_l1(P1m_c*S1m-P*Tm*Z1m_c+Y2_1m/mu,a/mu);
    E2m_c = prox_l1(P2m_c*S2m-P*Tm*Z2m_c+Y2_2m/mu,a/mu);
    E3m_c = prox_l1(P3m_c*S3m-P*Tm*Z3m_c+Y2_3m/mu,a/mu);
    E4m_c = prox_l1(P4m_c*S4m-P*Tm*Z4m_c+Y2_4m/mu,a/mu);
    E5m_c = prox_l1(P5m_c*S5m-P*Tm*Z5m_c+Y2_5m/mu,a/mu);
    E6m_c = prox_l1(P6m_c*S6m-P*Tm*Z6m_c+Y2_6m/mu,a/mu);
        % Update Ein
    E1n = prox_l1(P1n*S1n-P*Tn*Z1n+Y2_1n/mu,a/mu);
    E2n = prox_l1(P2n*S2n-P*Tn*Z2n+Y2_2n/mu,a/mu);
    E3n = prox_l1(P3n*S3n-P*Tn*Z3n+Y2_3n/mu,a/mu);
    E4n = prox_l1(P4n*S4n-P*Tn*Z4n+Y2_4n/mu,a/mu);
    E5n = prox_l1(P5n*S5n-P*Tn*Z5n+Y2_5n/mu,a/mu);
    E6n = prox_l1(P6n*S6n-P*Tn*Z6n+Y2_6n/mu,a/mu);
        % Update Zim_c
    G1_1_c = P1m_c*S1m - E1m_c + Y2_1m/mu;
    Z1m_c = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_1_c+F1m_c-Y1_1m/mu);
    G1_2_c = P2m_c*S2m - E2m_c + Y2_2m/mu;
    Z2m_c = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_2_c+F2m_c-Y1_2m/mu);
    G1_3_c = P3m_c*S3m - E3m_c + Y2_3m/mu;
    Z3m_c = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_3_c+F3m_c-Y1_3m/mu);
    G1_4_c = P4m_c*S4m - E4m_c + Y2_4m/mu;
    Z4m_c = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_4_c+F4m_c-Y1_4m/mu);
    G1_5_c = P5m_c*S5m - E5m_c + Y2_5m/mu;
    Z5m_c = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_5_c+F5m_c-Y1_5m/mu);
    G1_6_c = P6m_c*S6m - E6m_c + Y2_6m/mu;
    Z6m_c = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_6_c+F6m_c-Y1_6m/mu);
        % Update Zin
    G2_1 = P1n*S1n-E1n+Y2_1n/mu;
    Z1n = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_1+F1n-Y1_1n/mu);
    G2_2 = P2n*S2n-E2n+Y2_2n/mu;
    Z2n = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_2+F2n-Y1_2n/mu);
    G2_3 = P3n*S3n-E3n+Y2_3n/mu;
    Z3n = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_3+F3n-Y1_3n/mu);
    G2_4 = P4n*S4n-E4n+Y2_4n/mu;
    Z4n = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_4+F4n-Y1_4n/mu);
    G2_5 = P5n*S5n-E5n+Y2_5n/mu;
    Z5n = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_5+F5n-Y1_5n/mu);
    G2_6 = P6n*S6n-E6n+Y2_6n/mu;
    Z6n = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_6+F6n-Y1_6n/mu);
        % Update Pim_c
    G3_1_c = P1n*(S1n*S1n')*P1n';G4_1_c = P*Tm*Z1m_c+E1m_c-Y2_1m/mu;
    P1m_c = ((2*b*G3_1_c+mu*I)^-1)*mu*G4_1_c*S1m'*((S1m*S1m'+mibu)^-1);
    G3_2_c = P2n*(S2n*S2n')*P2n';G4_2_c = P*Tm*Z2m_c+E2m_c-Y2_2m/mu;
    P2m_c = ((2*b*G3_2_c+mu*I)^-1)*mu*G4_2_c*S2m'*((S2m*S2m'+mibu)^-1);
    G3_3_c = P3n*(S3n*S3n')*P3n';G4_3_c = P*Tm*Z3m_c+E3m_c-Y2_3m/mu;
    P3m_c = ((2*b*G3_3_c+mu*I)^-1)*mu*G4_3_c*S3m'*((S3m*S3m'+mibu)^-1);
    G3_4_c = P4n*(S4n*S4n')*P4n';G4_4_c = P*Tm*Z4m_c+E4m_c-Y2_4m/mu;
    P4m_c = ((2*b*G3_4_c+mu*I)^-1)*mu*G4_4_c*S4m'*((S4m*S4m'+mibu)^-1);
    G3_5_c = P5n*(S5n*S5n')*P5n';G4_5_c = P*Tm*Z5m_c+E5m_c-Y2_5m/mu;
    P5m_c = ((2*b*G3_5_c+mu*I)^-1)*mu*G4_5_c*S5m'*((S5m*S5m'+mibu)^-1);
    G3_6_c = P6n*(S6n*S6n')*P6n';G4_6_c = P*Tm*Z6m_c+E6m_c-Y2_6m/mu;
    P6m_c = ((2*b*G3_6_c+mu*I)^-1)*mu*G4_6_c*S6m'*((S6m*S6m'+mibu)^-1);
         % Update Pin
    G5_1 = P1m_c*(S1m*S1m')*P1m_c';G6_1 = P*Tn*Z1n+E1n-Y2_1n/mu;
    P1n = ((2*b*G5_1+mu*I)^-1)*mu*G6_1*S1n'*((S1n*S1n'+mibu)^-1);
    G5_2 = P2m_c*(S2m*S2m')*P2m_c';G6_2 = P*Tn*Z2n+E2n-Y2_2n/mu;
    P2n = ((2*b*G5_2+mu*I)^-1)*mu*G6_2*S2n'*((S2n*S2n'+mibu)^-1);
    G5_3 = P3m_c*(S3m*S3m')*P3m_c';G6_3 = P*Tn*Z3n+E3n-Y2_3n/mu;
    P3n = ((2*b*G5_3+mu*I)^-1)*mu*G6_3*S3n'*((S3n*S3n'+mibu)^-1);
    G5_4 = P4m_c*(S4m*S4m')*P4m_c';G6_4 = P*Tn*Z4n+E4n-Y2_4n/mu;
    P4n = ((2*b*G5_4+mu*I)^-1)*mu*G6_4*S4n'*((S4n*S4n'+mibu)^-1);
    G5_5 = P5m_c*(S5m*S5m')*P5m_c';G6_5 = P*Tn*Z5n+E5n-Y2_5n/mu;
    P5n = ((2*b*G5_5+mu*I)^-1)*mu*G6_5*S5n'*((S5n*S5n'+mibu)^-1);
    G5_6 = P6m_c*(S6m*S6m')*P6m_c';G6_6 = P*Tn*Z6n+E6n-Y2_6n/mu;
    P6n = ((2*b*G5_6+mu*I)^-1)*mu*G6_6*S6n'*((S6n*S6n'+mibu)^-1);
    
    % Pin s.t. Pin'Pin=I Orithogonal constraint 正交约束
    if cast(P1n'*P1n, 'int8') ~= eye(size(P1n)); P1n = orth(P1n); end
    if cast(P2n'*P2n, 'int8') ~= eye(size(P2n)); P2n = orth(P2n); end
    if cast(P3n'*P3n, 'int8') ~= eye(size(P3n)); P3n = orth(P3n); end
    if cast(P4n'*P4n, 'int8') ~= eye(size(P4n)); P4n = orth(P4n); end
    if cast(P5n'*P5n, 'int8') ~= eye(size(P5n)); P5n = orth(P5n); end
    if cast(P6n'*P6n, 'int8') ~= eye(size(P6n)); P6n = orth(P6n); end
    
    % Update P
    G7_1 = (P1m*S1m-E1m+Y2_1m/mu)*Z1m'*Tm'+(P1n*S1n-E1n+Y2_1n/mu)*Z1n'*Tn';
    G7_2 = (P2m*S2m-E2m+Y2_2m/mu)*Z2m'*Tm'+(P2n*S2n-E2n+Y2_2n/mu)*Z2n'*Tn';
    G7_3 = (P3m*S3m-E3m+Y2_3m/mu)*Z3m'*Tm'+(P3n*S3n-E3n+Y2_3n/mu)*Z3n'*Tn';
    G7_4 = (P4m*S4m-E4m+Y2_4m/mu)*Z4m'*Tm'+(P4n*S4n-E4n+Y2_4n/mu)*Z4n'*Tn';
    G7_5 = (P5m*S5m-E5m+Y2_5m/mu)*Z5m'*Tm'+(P5n*S5n-E5n+Y2_5n/mu)*Z5n'*Tn';
    G7_6 = (P6m*S6m-E6m+Y2_6m/mu)*Z6m'*Tm'+(P6n*S6n-E6n+Y2_6n/mu)*Z6n'*Tn';
    G7 = G7_1 + G7_2 + G7_3 + G7_4 + G7_5 + G7_6;
    G8_1 = Tm*(Z1m*Z1m')*Tm'+Tn*(Z1n*Z1n')*Tn';
    G8_2 = Tm*(Z2m*Z2m')*Tm'+Tn*(Z2n*Z2n')*Tn';
    G8_3 = Tm*(Z3m*Z3m')*Tm'+Tn*(Z3n*Z3n')*Tn';
    G8_4 = Tm*(Z4m*Z4m')*Tm'+Tn*(Z4n*Z4n')*Tn';
    G8_5 = Tm*(Z5m*Z5m')*Tm'+Tn*(Z5n*Z5n')*Tn';
    G8_6 = Tm*(Z6m*Z6m')*Tm'+Tn*(Z6n*Z6n')*Tn';
    G8 = G8_1 + G8_2 + G8_3 + G8_4 + G8_5 + G8_6;
    P = G7*((G8+mibu)^-1);
    
    % P s.t. P'P=I Orithogonal constraint
    if cast(P'*P, 'int8') ~= eye(size(P)); P = orth(P); end
    
    % calculate constraint matrix
    dY1_1m = Z1m - F1m;dY1_1n = Z1n - F1n;dY2_1m = P1m*S1m-P*Tm*Z1m-E1m;dY2_1n = P1n*S1n-P*Tn*Z1n-E1n;
    dY1_2m = Z2m - F2m;dY1_2n = Z2n - F2n;dY2_2m = P2m*S2m-P*Tm*Z2m-E2m;dY2_2n = P2n*S2n-P*Tn*Z2n-E2n;
    dY1_3m = Z3m - F3m;dY1_3n = Z3n - F3n;dY2_3m = P3m*S3m-P*Tm*Z3m-E3m;dY2_3n = P3n*S3n-P*Tn*Z3n-E3n;
    dY1_4m = Z4m - F4m;dY1_4n = Z4n - F4n;dY2_4m = P4m*S4m-P*Tm*Z4m-E4m;dY2_4n = P4n*S4n-P*Tn*Z4n-E4n;
    dY1_5m = Z5m - F5m;dY1_5n = Z5n - F5n;dY2_5m = P5m*S5m-P*Tm*Z5m-E5m;dY2_5n = P5n*S5n-P*Tn*Z5n-E5n;
    dY1_6m = Z6m - F6m;dY1_6n = Z6n - F6n;dY2_6m = P6m*S6m-P*Tm*Z6m-E6m;dY2_6n = P6n*S6n-P*Tn*Z6n-E6n;
    
    % Calculate the infinite norm of matrices
    inf_dY1_1m = norm(dY1_1m,inf);inf_dY1_1n = norm(dY1_1n,inf);inf_dY2_1m = norm(dY2_1m,inf);inf_dY2_1n = norm(dY2_1n,inf);
    inf_dY1_2m = norm(dY1_2m,inf);inf_dY1_2n = norm(dY1_2n,inf);inf_dY2_2m = norm(dY2_2m,inf);inf_dY2_2n = norm(dY2_2n,inf);
    inf_dY1_3m = norm(dY1_3m,inf);inf_dY1_3n = norm(dY1_3n,inf);inf_dY2_3m = norm(dY2_3m,inf);inf_dY2_3n = norm(dY2_3n,inf);
    inf_dY1_4m = norm(dY1_4m,inf);inf_dY1_4n = norm(dY1_4n,inf);inf_dY2_4m = norm(dY2_4m,inf);inf_dY2_4n = norm(dY2_4n,inf);
    inf_dY1_5m = norm(dY1_5m,inf);inf_dY1_5n = norm(dY1_5n,inf);inf_dY2_5m = norm(dY2_5m,inf);inf_dY2_5n = norm(dY2_5n,inf);
    inf_dY1_6m = norm(dY1_6m,inf);inf_dY1_6n = norm(dY1_6n,inf);inf_dY2_6m = norm(dY2_6m,inf);inf_dY2_6n = norm(dY2_6n,inf);
    
    % Calculate iterations loss 计算迭代损失
    if DEBUG
        if iter == 1 || mod(iter,2) == 0
            F_norm_1 = nuclearnormF1m + nuclearnormF1n + nuclearnormF2m + nuclearnormF2n + nuclearnormF3m + nuclearnormF3n + nuclearnormF4m + nuclearnormF4n + nuclearnormF5m + nuclearnormF5n + nuclearnormF6m + nuclearnormF6n;
            Es_l1_1 = a*(norm(E1m(:),1)+norm(E1n(:),1)+norm(E2m(:),1)+norm(E2n(:),1)+norm(E3m(:),1)+norm(E3n(:),1)+norm(E4m(:),1)+norm(E4n(:),1)+norm(E5m(:),1)+norm(E5n(:),1)+norm(E6m(:),1)+norm(E6n(:),1));
            CI_term_1_1 = norm((P1n*S1n)'*(P1m*S1m),'fro') + norm((P1m*S1m)'*(P1n*S1n),'fro');
            CI_term_1_2 = norm((P2n*S2n)'*(P2m*S2m),'fro') + norm((P2m*S2m)'*(P2n*S2n),'fro');
            CI_term_1_3 = norm((P3n*S3n)'*(P3m*S3m),'fro') + norm((P3m*S3m)'*(P3n*S3n),'fro');
            CI_term_1_4 = norm((P4n*S4n)'*(P4m*S4m),'fro') + norm((P4m*S4m)'*(P4n*S4n),'fro');
            CI_term_1_5 = norm((P5n*S5n)'*(P5m*S5m),'fro') + norm((P5m*S5m)'*(P5n*S5n),'fro');
            CI_term_1_6 = norm((P6n*S6n)'*(P6m*S6m),'fro') + norm((P6m*S6m)'*(P6n*S6n),'fro');
            CI_term_1 = (b/2)*(CI_term_1_1 + CI_term_1_2 + CI_term_1_3 + CI_term_1_4 + CI_term_1_5 + CI_term_1_6);
            % Lagrange multiplier subterm
            LgM_term_1_1 = trace(Y1_1m'*dY1_1m) + trace(Y1_1n'*dY1_1n) + trace(Y2_1m'*dY2_1m) + trace(Y2_1n'*dY2_1n);
            LgM_term_1_2 = trace(Y1_2m'*dY1_2m) + trace(Y1_2n'*dY1_2n) + trace(Y2_2m'*dY2_2m) + trace(Y2_2n'*dY2_2n);
            LgM_term_1_3 = trace(Y1_3m'*dY1_3m) + trace(Y1_3n'*dY1_3n) + trace(Y2_3m'*dY2_3m) + trace(Y2_3n'*dY2_3n);
            LgM_term_1_4 = trace(Y1_4m'*dY1_4m) + trace(Y1_4n'*dY1_4n) + trace(Y2_4m'*dY2_4m) + trace(Y2_4n'*dY2_4n);
            LgM_term_1_5 = trace(Y1_5m'*dY1_5m) + trace(Y1_5n'*dY1_5n) + trace(Y2_5m'*dY2_5m) + trace(Y2_5n'*dY2_5n);
            LgM_term_1_6 = trace(Y1_6m'*dY1_6m) + trace(Y1_6n'*dY1_6n) + trace(Y2_6m'*dY2_6m) + trace(Y2_6n'*dY2_6n);
            LgM_term_1 = LgM_term_1_1 + LgM_term_1_2 + LgM_term_1_3 + LgM_term_1_4 + LgM_term_1_5 + LgM_term_1_6;
            % penalty term
            Penalty_term_1_1 = norm(dY1_1m,'fro')^2 + norm(dY1_1n,'fro')^2 + norm(dY2_1m,'fro')^2 + norm(dY2_1n,'fro')^2;
            Penalty_term_1_2 = norm(dY1_2m,'fro')^2 + norm(dY1_2n,'fro')^2 + norm(dY2_2m,'fro')^2 + norm(dY2_2n,'fro')^2;
            Penalty_term_1_3 = norm(dY1_3m,'fro')^2 + norm(dY1_3n,'fro')^2 + norm(dY2_3m,'fro')^2 + norm(dY2_3n,'fro')^2;
            Penalty_term_1_4 = norm(dY1_4m,'fro')^2 + norm(dY1_4n,'fro')^2 + norm(dY2_4m,'fro')^2 + norm(dY2_4n,'fro')^2;
            Penalty_term_1_5 = norm(dY1_5m,'fro')^2 + norm(dY1_5n,'fro')^2 + norm(dY2_5m,'fro')^2 + norm(dY2_5n,'fro')^2;
            Penalty_term_1_6 = norm(dY1_6m,'fro')^2 + norm(dY1_6n,'fro')^2 + norm(dY2_6m,'fro')^2 + norm(dY2_6n,'fro')^2;
            Penalty_term_1 = (mu/2)*(Penalty_term_1_1 + Penalty_term_1_2 + Penalty_term_1_3 + Penalty_term_1_4 + Penalty_term_1_5 + Penalty_term_1_6);
            obj_1 = F_norm_1 + Es_l1_1 + CI_term_1 + LgM_term_1 + Penalty_term_1;
            err_1 = sqrt(Penalty_term_1_1 + Penalty_term_1_2 + Penalty_term_1_3 + Penalty_term_1_4 + Penalty_term_1_5 + Penalty_term_1_6);
            
            disp(['iter' num2str(iter) ', mu=' num2str(mu) ...
                  ', obj=' num2str(obj_1) ', err=' num2str(err_1)]);  
        end
    end
    
    % Determining convergence or not
    if inf_dY1_1m < tol && inf_dY1_1n < tol && inf_dY2_1m < tol && inf_dY2_1n < to1 && inf_dY1_2m < tol && inf_dY1_2n < tol && inf_dY2_2m < tol && inf_dY2_2n < tol && inf_dY1_3m < tol && inf_dY1_3n < tol && inf_dY2_3m < tol && inf_dY2_3n < tol && inf_dY1_4m < tol && inf_dY1_4n < tol && inf_dY2_4m < tol && inf_dY2_4n < tol && inf_dY1_5m < tol && inf_dY1_5n < tol && inf_dY2_5m < tol && inf_dY2_5n < tol && inf_dY1_6m < tol && inf_dY1_6n < tol && inf_dY2_6m < tol && inf_dY2_6n < tol
        break;
    end
    
     % Update the multipliers and parameters
    Y1_1m = Y1_1m + mu*dY1_1m;Y1_1n = Y1_1n + mu*dY1_1n;Y2_1m = Y2_1m + mu*dY2_1m;Y2_1n = Y2_1n + mu*dY2_1n;
    Y1_2m = Y1_2m + mu*dY1_2m;Y1_2n = Y1_2n + mu*dY1_2n;Y2_2m = Y2_2m + mu*dY2_2m;Y2_2n = Y2_2n + mu*dY2_2n;
    Y1_3m = Y1_3m + mu*dY1_3m;Y1_3n = Y1_3n + mu*dY1_3n;Y2_3m = Y2_3m + mu*dY2_3m;Y2_3n = Y2_3n + mu*dY2_3n;
    Y1_4m = Y1_4m + mu*dY1_4m;Y1_4n = Y1_4n + mu*dY1_4n;Y2_4m = Y2_4m + mu*dY2_4m;Y2_4n = Y2_4n + mu*dY2_4n;
    Y1_5m = Y1_5m + mu*dY1_5m;Y1_5n = Y1_5n + mu*dY1_5n;Y2_5m = Y2_5m + mu*dY2_5m;Y2_5n = Y2_5n + mu*dY2_5n;
    Y1_6m = Y1_6m + mu*dY1_6m;Y1_6n = Y1_6n + mu*dY1_6n;Y2_6m = Y2_6m + mu*dY2_6m;Y2_6n = Y2_6n + mu*dY2_6n;
    
    mu = min(mu*rho, max_mu);
end
F_norm = nuclearnormF1m + nuclearnormF1n + nuclearnormF2m + nuclearnormF2n + nuclearnormF3m + nuclearnormF3n + nuclearnormF4m + nuclearnormF4n + nuclearnormF5m + nuclearnormF5n + nuclearnormF6m + nuclearnormF6n;
Es_l1 = a*(norm(E1m(:),1)+norm(E1n(:),1)+norm(E2m(:),1)+norm(E2n(:),1)+norm(E3m(:),1)+norm(E3n(:),1)+norm(E4m(:),1)+norm(E4n(:),1)+norm(E5m(:),1)+norm(E5n(:),1)+norm(E6m(:),1)+norm(E6n(:),1));
CI_term_1 = norm((P1n*S1n)'*(P1m*S1m),'fro') + norm((P1m*S1m)'*(P1n*S1n),'fro');
CI_term_2 = norm((P2n*S2n)'*(P2m*S2m),'fro') + norm((P2m*S2m)'*(P2n*S2n),'fro');
CI_term_3 = norm((P3n*S3n)'*(P3m*S3m),'fro') + norm((P3m*S3m)'*(P3n*S3n),'fro');
CI_term_4 = norm((P4n*S4n)'*(P4m*S4m),'fro') + norm((P4m*S4m)'*(P4n*S4n),'fro');
CI_term_5 = norm((P5n*S5n)'*(P5m*S5m),'fro') + norm((P5m*S5m)'*(P5n*S5n),'fro');
CI_term_6 = norm((P6n*S6n)'*(P6m*S6m),'fro') + norm((P6m*S6m)'*(P6n*S6n),'fro');
CI_term = (b/2)*(CI_term_1 + CI_term_2 + CI_term_3 + CI_term_4 + CI_term_5 + CI_term_6);
% Lagrange multiplier subterm 拉格朗日乘子项
LgM_term_1 = trace(Y1_1m'*dY1_1m) + trace(Y1_1n'*dY1_1n) + trace(Y2_1m'*dY2_1m) + trace(Y2_1n'*dY2_1n);
LgM_term_2 = trace(Y1_2m'*dY1_2m) + trace(Y1_2n'*dY1_2n) + trace(Y2_2m'*dY2_2m) + trace(Y2_2n'*dY2_2n);
LgM_term_3 = trace(Y1_3m'*dY1_3m) + trace(Y1_3n'*dY1_3n) + trace(Y2_3m'*dY2_3m) + trace(Y2_3n'*dY2_3n);
LgM_term_4 = trace(Y1_4m'*dY1_4m) + trace(Y1_4n'*dY1_4n) + trace(Y2_4m'*dY2_4m) + trace(Y2_4n'*dY2_4n);
LgM_term_5 = trace(Y1_5m'*dY1_5m) + trace(Y1_5n'*dY1_5n) + trace(Y2_5m'*dY2_5m) + trace(Y2_5n'*dY2_5n);
LgM_term_6 = trace(Y1_6m'*dY1_6m) + trace(Y1_6n'*dY1_6n) + trace(Y2_6m'*dY2_6m) + trace(Y2_6n'*dY2_6n);
LgM_term = LgM_term_1 + LgM_term_2 + LgM_term_3 + LgM_term_4 + LgM_term_5 + LgM_term_6;
% penalty term 惩罚项
Penalty_term_1 = norm(dY1_1m,'fro')^2 + norm(dY1_1n,'fro')^2 + norm(dY2_1m,'fro')^2 + norm(dY2_1n,'fro')^2;
Penalty_term_2 = norm(dY1_2m,'fro')^2 + norm(dY1_2n,'fro')^2 + norm(dY2_2m,'fro')^2 + norm(dY2_2n,'fro')^2;
Penalty_term_3 = norm(dY1_3m,'fro')^2 + norm(dY1_3n,'fro')^2 + norm(dY2_3m,'fro')^2 + norm(dY2_3n,'fro')^2;
Penalty_term_4 = norm(dY1_4m,'fro')^2 + norm(dY1_4n,'fro')^2 + norm(dY2_4m,'fro')^2 + norm(dY2_4n,'fro')^2;
Penalty_term_5 = norm(dY1_5m,'fro')^2 + norm(dY1_5n,'fro')^2 + norm(dY2_5m,'fro')^2 + norm(dY2_5n,'fro')^2;
Penalty_term_6 = norm(dY1_6m,'fro')^2 + norm(dY1_6n,'fro')^2 + norm(dY2_6m,'fro')^2 + norm(dY2_6n,'fro')^2;
Penalty_term = (mu/2)*(Penalty_term_1 + Penalty_term_2 + Penalty_term_3 + Penalty_term_4 + Penalty_term_5 + Penalty_term_6);
% 目标函数值
obj = F_norm + Es_l1 + CI_term + LgM_term + Penalty_term;
err = sqrt(Penalty_term_1 + Penalty_term_2 + Penalty_term_3 + Penalty_term_4 + Penalty_term_5 + Penalty_term_6);
% 保存变量
F_sm.F1m = F1m;F_sm.F2m = F2m;F_sm.F3m = F3m;F_sm.F4m = F4m;F_sm.F5m = F5m;F_sm.F6m = F6m;
F_sn.F1n = F1n;F_sn.F2n = F2n;F_sn.F3n = F3n;F_sn.F4n = F4n;F_sn.F5n = F5n;F_sn.F6n = F6n;
Z_sm.Z1m = Z1m;Z_sm.Z2m = Z2m;Z_sm.Z3m = Z3m;Z_sm.Z4m = Z4m;Z_sm.Z5m = Z5m;Z_sm.Z6m = Z6m;
Z_sn.Z1n = Z1n;Z_sn.Z2n = Z2n;Z_sn.Z3n = Z3n;Z_sn.Z4n = Z4n;Z_sn.Z5n = Z5n;Z_sn.Z6n = Z6n;
P_sm.P1m = P1m;P_sm.P2m = P2m;P_sm.P3m = P3m;P_sm.P4m = P4m;P_sm.P5m = P5m;P_sm.P6m = P6m;
P_sn.P1n = P1n;P_sn.P2n = P2n;P_sn.P3n = P3n;P_sn.P4n = P4n;P_sn.P5n = P5n;P_sn.P6n = P6n;
Es_sm.E1m = E1m;Es_sm.E2m = E2m;Es_sm.E3m = E3m;Es_sm.E4m = E4m;Es_sm.E5m = E5m;Es_sm.E6m = E6m;
Es_sn.E1n = E1n;Es_sn.E2n = E2n;Es_sn.E3n = E3n;Es_sn.E4n = E4n;Es_sn.E5n = E5n;Es_sn.E6n = E6n;
end