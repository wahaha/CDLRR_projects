function [F_sm,F_sn,Z_sm,Z_sn,W_sm,W_sn,Es_sm,Es_sn,obj,err] = CDLRR_3(Tm,Tn,Ssm,Ssn,a,b,opts)
%Initialize
to1 = 1e-8;
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;
DEBUG = 0;

if ~exist('opts', 'var')
    opts = [];
end 
if isfield(opts, 'tol');        tol = opts.tol;            end
if isfield(opts, 'max_iter');   max_iter = opts.max_iter;  end
if isfield(opts, 'max_mu');     max_mu = opts.max_mu;      end
if isfield(opts, 'rho');        rho = opts.rho;            end
if isfield(opts, 'mu');         mu = opts.mu;              end
if isfield(opts, 'DEBUG');      DEBUG = opts.DEBUG;        end

[~,nTm] = size(Tm);
[d,nTn] = size(Tn);

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
mibu = 1e-7*I;

Z1m = zeros(nTm,nS1m);Z1n = zeros(nTn,nS1n);
Z2m = zeros(nTm,nS2m);Z2n = zeros(nTn,nS2n);
Z3m = zeros(nTm,nS3m);Z3n = zeros(nTn,nS3n);
Z4m = zeros(nTm,nS4m);Z4n = zeros(nTn,nS4n);
Z5m = zeros(nTm,nS5m);Z5n = zeros(nTn,nS5n);
Z6m = zeros(nTm,nS6m);Z6n = zeros(nTn,nS6n);
Z1m_c = Z1m;Z2m_c = Z2m;Z3m_c = Z3m;Z4m_c = Z4m;Z5m_c = Z5m;Z6m_c = Z6m;
Z1n_c = Z1n;Z2n_c = Z2n;Z3n_c = Z3n;Z4n_c = Z4n;Z5n_c = Z5n;Z6n_c = Z6n;

F1m = zeros(nTm,nS1m);F1n = zeros(nTn,nS1n);
F2m = zeros(nTm,nS2m);F2n = zeros(nTn,nS2n);
F3m = zeros(nTm,nS3m);F3n = zeros(nTn,nS3n);
F4m = zeros(nTm,nS4m);F4n = zeros(nTn,nS4n);
F5m = zeros(nTm,nS5m);F5n = zeros(nTn,nS5n);
F6m = zeros(nTm,nS6m);F6n = zeros(nTn,nS6n);
F1m_c = F1m;F1n_c = F1n;
F2m_c = F2m;F2n_c = F2n;
F3m_c = F3m;F3n_c = F3n;
F4m_c = F4m;F4n_c = F4n;
F5m_c = F5m;F5n_c = F5n;
F6m_c = F6m;F6n_c = F6n;

Z_mI = eye(nTm);
Z_nI = eye(nTn);

Es1m = zeros(d,nS1m);Es1n = zeros(d,nS1n);
Es2m = zeros(d,nS2m);Es2n = zeros(d,nS2n);
Es3m = zeros(d,nS3m);Es3n = zeros(d,nS3n);
Es4m = zeros(d,nS4m);Es4n = zeros(d,nS4n);
Es5m = zeros(d,nS5m);Es5n = zeros(d,nS5n);
Es6m = zeros(d,nS6m);Es6n = zeros(d,nS6n);
Es1m_c = Es1m;Es1n_c = Es1n;
Es2m_c = Es2m;Es2n_c = Es2n;
Es3m_c = Es3m;Es3n_c = Es3n;
Es4m_c = Es4m;Es4n_c = Es4n;
Es5m_c = Es5m;Es5n_c = Es5n;
Es6m_c = Es6m;Es6n_c = Es6n;

W1m = I;W1n = I;
W2m = I;W2n = I;
W3m = I;W3n = I;
W4m = I;W4n = I;
W5m = I;W5n = I;
W6m = I;W6n = I;
W1m_c = W1m;W1n_c = W1n;
W2m_c = W2m;W2n_c = W2n;
W3m_c = W3m;W3n_c = W3n;
W4m_c = W4m;W4n_c = W4n;
W5m_c = W5m;W5n_c = W5n;
W6m_c = W6m;W6n_c = W6n;

% lagrangian multiplier
Y1_1m = zeros(nTm,nS1m);Y1_2m = zeros(nTm,nS2m);Y1_3m = zeros(nTm,nS3m);Y1_4m = zeros(nTm,nS4m);Y1_5m = zeros(nTm,nS5m);Y1_6m = zeros(nTm,nS6m);
Y1_1n = zeros(nTn,nS1n);Y1_2n = zeros(nTn,nS2n);Y1_3n = zeros(nTn,nS3n);Y1_4n = zeros(nTn,nS4n);Y1_5n = zeros(nTn,nS5n);Y1_6n = zeros(nTn,nS6n);
Y2_1m = zeros(d,nS1m);Y2_2m = zeros(d,nS2m);Y2_3m = zeros(d,nS3m);Y2_4m = zeros(d,nS4m);Y2_5m = zeros(d,nS5m);Y2_6m = zeros(d,nS6m);
Y2_1n = zeros(d,nS1n);Y2_2n = zeros(d,nS2n);Y2_3n = zeros(d,nS3n);Y2_4n = zeros(d,nS4n);Y2_5n = zeros(d,nS5n);Y2_6n = zeros(d,nS6n);

for iter = 1:max_iter
   % Update Win
   % first super block {Fm1,Fm2}
   % Fmi = argmin_{Fmi} ||Fmi||_*+<Y1_mi,Zmi-Fmi>+(mu/2)*||Zmi-Fmi||_F~2
   % proximal operator
   [F1m,~] = prox_nuclear(Z1m+Y1_1m/mu,1/mu);
   [F2m,~] = prox_nuclear(Z2m+Y1_2m/mu,1/mu);
   [F3m,~] = prox_nuclear(Z3m+Y1_3m/mu,1/mu);
   [F4m,~] = prox_nuclear(Z4m+Y1_4m/mu,1/mu);
   [F5m,~] = prox_nuclear(Z5m+Y1_5m/mu,1/mu);
   [F6m,~] = prox_nuclear(Z6m+Y1_6m/mu,1/mu);
   [F1n_c,nuclearnormF1n_c] = prox_nuclear(Z1n_c+Y1_1n/mu,1/mu);
   [F2n_c,nuclearnormF2n_c] = prox_nuclear(Z2n_c+Y1_2n/mu,1/mu);
   [F3n_c,nuclearnormF3n_c] = prox_nuclear(Z3n_c+Y1_3n/mu,1/mu);
   [F4n_c,nuclearnormF4n_c] = prox_nuclear(Z4n_c+Y1_4n/mu,1/mu);
   [F5n_c,nuclearnormF5n_c] = prox_nuclear(Z5n_c+Y1_5n/mu,1/mu);
   [F6n_c,nuclearnormF6n_c] = prox_nuclear(Z6n_c+Y1_6n/mu,1/mu);

   
   % proximal operator
   Es1m = prox_l1(W1m*S1m-Tm*Z1m+Y2_1m/mu,a/mu);
   Es2m = prox_l1(W2m*S2m-Tm*Z2m+Y2_2m/mu,a/mu);
   Es3m = prox_l1(W3m*S3m-Tm*Z3m+Y2_3m/mu,a/mu);
   Es4m = prox_l1(W4m*S4m-Tm*Z4m+Y2_4m/mu,a/mu);
   Es5m = prox_l1(W5m*S5m-Tm*Z5m+Y2_5m/mu,a/mu);
   Es6m = prox_l1(W6m*S6m-Tm*Z6m+Y2_6m/mu,a/mu);

   Es1n_c = prox_l1(W1n_c*S1n-Tn*Z1n_c+Y2_1n/mu,a/mu);
   Es2n_c = prox_l1(W2n_c*S2n-Tn*Z2n_c+Y2_2n/mu,a/mu);
   Es3n_c = prox_l1(W3n_c*S3n-Tn*Z3n_c+Y2_3n/mu,a/mu);
   Es4n_c = prox_l1(W4n_c*S4n-Tn*Z4n_c+Y2_4n/mu,a/mu);
   Es5n_c = prox_l1(W5n_c*S5n-Tn*Z5n_c+Y2_5n/mu,a/mu);
   Es6n_c = prox_l1(W6n_c*S6n-Tn*Z6n_c+Y2_6n/mu,a/mu);

   
   G1_1 = W1m*S1m-Es1m+Y2_1m/mu;
   Z1m = ((Z_mI+Tm'*Tm)^-1)*(F1m-Y1_1m/mu+Tm'*G1_1);
   G1_2 = W2m*S2m-Es2m+Y2_2m/mu;
   Z2m = ((Z_mI+Tm'*Tm)^-1)*(F2m-Y1_2m/mu+Tm'*G1_2);
   G1_3 = W3m*S3m-Es3m+Y2_3m/mu;
   Z3m = ((Z_mI+Tm'*Tm)^-1)*(F3m-Y1_3m/mu+Tm'*G1_3);
   G1_4 = W4m*S4m-Es4m+Y2_4m/mu;
   Z4m = ((Z_mI+Tm'*Tm)^-1)*(F4m-Y1_4m/mu+Tm'*G1_4);
   G1_5 = W5m*S5m-Es5m+Y2_5m/mu;
   Z5m = ((Z_mI+Tm'*Tm)^-1)*(F5m-Y1_5m/mu+Tm'*G1_5);
   G1_6 = W6m*S6m-Es6m+Y2_6m/mu;
   Z6m = ((Z_mI+Tm'*Tm)^-1)*(F6m-Y1_6m/mu+Tm'*G1_6); 
       
   G2_1_c = W1n_c*S1n-Es1n_c+Y2_1n/mu;
   Z1n_c = ((Z_nI+Tn'*Tn)^-1)*(F1n_c-Y1_1n/mu+Tn'*G2_1_c);
   G2_2_c = W2n_c*S2n-Es2n_c+Y2_2n/mu;
   Z2n_c = ((Z_nI+Tn'*Tn)^-1)*(F2n_c-Y1_2n/mu+Tn'*G2_2_c);
   G2_3_c = W3n_c*S3n-Es3n_c+Y2_3n/mu;
   Z3n_c = ((Z_nI+Tn'*Tn)^-1)*(F3n_c-Y1_3n/mu+Tn'*G2_3_c);
   G2_4_c = W4n_c*S4n-Es4n_c+Y2_4n/mu;
   Z4n_c = ((Z_nI+Tn'*Tn)^-1)*(F4n_c-Y1_4n/mu+Tn'*G2_4_c);
   G2_5_c = W5n_c*S5n-Es5n_c+Y2_5n/mu;
   Z5n_c = ((Z_nI+Tn'*Tn)^-1)*(F5n_c-Y1_5n/mu+Tn'*G2_5_c);
   G2_6_c = W6n_c*S6n-Es6n_c+Y2_6n/mu;
   Z6n_c = ((Z_nI+Tn'*Tn)^-1)*(F6n_c-Y1_6n/mu+Tn'*G2_6_c);
   
   G3_1 = W1n_c*(S1n*S1n')*W1n_c';
   W1m = mu*((2*b*G3_1+mu*I)^-1)*((Tm*Z1m+Es1m-Y2_1m/mu)*S1m')*((S1m*S1m'+mibu)^-1);
   G3_2 = W2n_c*(S2n*S2n')*W2n_c';
   W2m = mu*((2*b*G3_2+mu*I)^-1)*((Tm*Z2m+Es2m-Y2_2m/mu)*S2m')*((S2m*S2m'+mibu)^-1);
   G3_3 = W3n_c*(S3n*S3n')*W3n_c';
   W3m = mu*((2*b*G3_3+mu*I)^-1)*((Tm*Z3m+Es3m-Y2_3m/mu)*S3m')*((S3m*S3m'+mibu)^-1);
   G3_4 = W4n_c*(S4n*S4n')*W4n_c';
   W4m = mu*((2*b*G3_4+mu*I)^-1)*((Tm*Z4m+Es4m-Y2_4m/mu)*S4m')*((S4m*S4m'+mibu)^-1);
   G3_5 = W5n_c*(S5n*S5n')*W5n_c';
   W5m = mu*((2*b*G3_5+mu*I)^-1)*((Tm*Z5m+Es5m-Y2_5m/mu)*S5m')*((S5m*S5m'+mibu)^-1);
   G3_6 = W6n_c*(S6n*S6n')*W6n_c';
   W6m = mu*((2*b*G3_6+mu*I)^-1)*((Tm*Z6m+Es6m-Y2_6m/mu)*S6m')*((S6m*S6m'+mibu)^-1);
   
   G4_1_c = W1m*(S1m*S1m')*W1m';
   W1n_c = mu*((2*b*G4_1_c+mu*I)^-1)*((Tn*Z1n_c+Es1n_c-Y2_1n/mu)*S1n')*((S1n*S1n'+mibu)^-1);
   G4_2_c = W2m*(S2m*S2m')*W2m';
   W2n_c = mu*((2*b*G4_2_c+mu*I)^-1)*((Tn*Z2n_c+Es2n_c-Y2_2n/mu)*S2n')*((S2n*S2n'+mibu)^-1);
   G4_3_c = W3m*(S3m*S3m')*W3m';
   W3n_c = mu*((2*b*G4_3_c+mu*I)^-1)*((Tn*Z3n_c+Es3n_c-Y2_3n/mu)*S3n')*((S3n*S3n'+mibu)^-1);
   G4_4_c = W4m*(S4m*S4m')*W4m';
   W4n_c = mu*((2*b*G4_4_c+mu*I)^-1)*((Tn*Z4n_c+Es4n_c-Y2_4n/mu)*S4n')*((S4n*S4n'+mibu)^-1);
   G4_5_c = W5m*(S5m*S5m')*W5m';
   W5n_c = mu*((2*b*G4_5_c+mu*I)^-1)*((Tn*Z5n_c+Es5n_c-Y2_5n/mu)*S5n')*((S5n*S5n'+mibu)^-1);
   G4_6_c = W6m*(S6m*S6m')*W6m';
   W6n_c = mu*((2*b*G4_6_c+mu*I)^-1)*((Tn*Z6n_c+Es6n_c-Y2_6n/mu)*S6n')*((S6n*S6n'+mibu)^-1);
   
    if ~isequal(cast(W1n_c'*W1n_c,'int8'),eye(size(W1n_c))); W1n_c = orth(W1n_c); end
    if ~isequal(cast(W2n_c'*W2n_c,'int8'),eye(size(W2n_c))); W2n_c = orth(W2n_c); end
    if ~isequal(cast(W3n_c'*W3n_c,'int8'),eye(size(W3n_c))); W3n_c = orth(W3n_c); end
    if ~isequal(cast(W4n_c'*W4n_c,'int8'),eye(size(W4n_c))); W4n_c = orth(W4n_c); end
    if ~isequal(cast(W5n_c'*W5n_c,'int8'),eye(size(W5n_c))); W5n_c = orth(W5n_c); end
    if ~isequal(cast(W6n_c'*W6n_c,'int8'),eye(size(W6n_c))); W6n_c = orth(W6n_c); end
    if ~isequal(cast(W7n_c'*W7n_c,'int8'),eye(size(W7n_c))); W7n_c = orth(W7n_c); end

   % Update Wim
   [F1n,~] = prox_nuclear(Z1n+Y1_1n/mu,1/mu);
   [F2n,~] = prox_nuclear(Z2n+Y1_2n/mu,1/mu);
   [F3n,~] = prox_nuclear(Z3n+Y1_3n/mu,1/mu);
   [F4n,~] = prox_nuclear(Z4n+Y1_4n/mu,1/mu);
   [F5n,~] = prox_nuclear(Z5n+Y1_5n/mu,1/mu);
   [F6n,~] = prox_nuclear(Z6n+Y1_6n/mu,1/mu);
   [F1m_c,nuclearnormF1m_c] = prox_nuclear(Z1m_c+Y1_1m/mu,1/mu);
   [F2m_c,nuclearnormF2m_c] = prox_nuclear(Z2m_c+Y1_2m/mu,1/mu);
   [F3m_c,nuclearnormF3m_c] = prox_nuclear(Z3m_c+Y1_3m/mu,1/mu);
   [F4m_c,nuclearnormF4m_c] = prox_nuclear(Z4m_c+Y1_4m/mu,1/mu);
   [F5m_c,nuclearnormF5m_c] = prox_nuclear(Z5m_c+Y1_5m/mu,1/mu);
   [F6m_c,nuclearnormF6m_c] = prox_nuclear(Z6m_c+Y1_6m/mu,1/mu);
   
   Es1n = prox_l1(W1n*S1n-Tn*Z1n+Y2_1n/mu,a/mu);
   Es2n = prox_l1(W2n*S2n-Tn*Z2n+Y2_2n/mu,a/mu); 
   Es3n = prox_l1(W3n*S3n-Tn*Z3n+Y2_3n/mu,a/mu);
   Es4n = prox_l1(W4n*S4n-Tn*Z4n+Y2_4n/mu,a/mu);
   Es5n = prox_l1(W5n*S5n-Tn*Z5n+Y2_5n/mu,a/mu);
   Es6n = prox_l1(W6n*S6n-Tn*Z6n+Y2_6n/mu,a/mu);
   Es1m_c = prox_l1(W1m_c*S1m-Tm*Z1m_c+Y2_1m/mu,a/mu);
   Es2m_c = prox_l1(W2m_c*S2m-Tm*Z2m_c+Y2_2m/mu,a/mu);
   Es3m_c = prox_l1(W3m_c*S3m-Tm*Z3m_c+Y2_3m/mu,a/mu);
   Es4m_c = prox_l1(W4m_c*S4m-Tm*Z4m_c+Y2_4m/mu,a/mu);
   Es5m_c = prox_l1(W5m_c*S5m-Tm*Z5m_c+Y2_5m/mu,a/mu);
   Es6m_c = prox_l1(W6m_c*S6m-Tm*Z6m_c+Y2_6m/mu,a/mu);
   
   G2_1 = W1n*S1n-Es1n+Y2_1n/mu;
   Z1n = ((Z_nI+Tn'*Tn)^-1)*(F1n-Y1_1n/mu+Tn'*G2_1);
   G2_2 = W2n*S2n-Es2n+Y2_2n/mu;
   Z2n = ((Z_nI+Tn'*Tn)^-1)*(F2n-Y1_2n/mu+Tn'*G2_2);
   G2_3 = W3n*S3n-Es3n+Y2_3n/mu;
   Z3n = ((Z_nI+Tn'*Tn)^-1)*(F3n-Y1_3n/mu+Tn'*G2_3);
   G2_4 = W4n*S4n-Es4n+Y2_4n/mu;
   Z4n = ((Z_nI+Tn'*Tn)^-1)*(F4n-Y1_4n/mu+Tn'*G2_4);
   G2_5 = W5n*S5n-Es5n+Y2_5n/mu;
   Z5n = ((Z_nI+Tn'*Tn)^-1)*(F5n-Y1_5n/mu+Tn'*G2_5);
   G2_6 = W6n*S6n-Es6n+Y2_6n/mu;
   Z6n = ((Z_nI+Tn'*Tn)^-1)*(F6n-Y1_6n/mu+Tn'*G2_6);
   
   G1_1_c = W1m_c*S1m-Es1m_c+Y2_1m/mu;
   Z1m_c = ((Z_mI+Tm'*Tm)^-1)*(F1m_c-Y1_1m/mu+Tm'*G1_1_c);
   G1_2_c = W2m_c*S2m-Es2m_c+Y2_2m/mu;
   Z2m_c = ((Z_mI+Tm'*Tm)^-1)*(F2m_c-Y1_2m/mu+Tm'*G1_2_c);
   G1_3_c = W3m_c*S3m-Es3m_c+Y2_3m/mu;
   Z3m_c = ((Z_mI+Tm'*Tm)^-1)*(F3m_c-Y1_3m/mu+Tm'*G1_3_c);
   G1_4_c = W4m_c*S4m-Es4m_c+Y2_4m/mu;
   Z4m_c = ((Z_mI+Tm'*Tm)^-1)*(F4m_c-Y1_4m/mu+Tm'*G1_4_c);
   G1_5_c = W5m_c*S5m-Es5m_c+Y2_5m/mu;
   Z5m_c = ((Z_mI+Tm'*Tm)^-1)*(F5m_c-Y1_5m/mu+Tm'*G1_5_c);
   G1_6_c = W6m_c*S6m-Es6m_c+Y2_6m/mu;
   Z6m_c = ((Z_mI+Tm'*Tm)^-1)*(F6m_c-Y1_6m/mu+Tm'*G1_6_c);
   
   G4_1 = W1m_c*(S1m*S1m')*W1m_c';
   W1n = mu*((2*b*G4_1+mu*I)^-1)*((Tn*Z1n+Es1n-Y2_1n/mu)*S1n')*((S1n*S1n'+mibu)^-1);
   G4_2 = W2m_c*(S2m*S2m')*W2m_c';
   W2n = mu*((2*b*G4_2+mu*I)^-1)*((Tn*Z2n+Es2n-Y2_2n/mu)*S2n')*((S2n*S2n'+mibu)^-1);
   G4_3 = W3m_c*(S3m*S3m')*W3m_c';
   W3n = mu*((2*b*G4_3+mu*I)^-1)*((Tn*Z3n+Es3n-Y2_3n/mu)*S3n')*((S3n*S3n'+mibu)^-1);
   G4_4 = W4m_c*(S4m*S4m')*W4m_c';
   W4n = mu*((2*b*G4_4+mu*I)^-1)*((Tn*Z4n+Es4n-Y2_4n/mu)*S4n')*((S4n*S4n'+mibu)^-1);
   G4_5 = W5m_c*(S5m*S5m')*W5m_c';
   W5n = mu*((2*b*G4_5+mu*I)^-1)*((Tn*Z5n+Es5n-Y2_5n/mu)*S5n')*((S5n*S5n'+mibu)^-1);
   G4_6 = W6m_c*(S6m*S6m')*W6m_c';
   W6n = mu*((2*b*G4_6+mu*I)^-1)*((Tn*Z6n+Es6n-Y2_6n/mu)*S6n')*((S6n*S6n'+mibu)^-1);

   
   G3_1_c = W1n*(S1n*S1n')*W1n';
   W1m_c = mu*((2*b*G3_1_c+mu*I)^-1)*((Tm*Z1m_c+Es1m_c-Y2_1m/mu)*S1m')*((S1m*S1m'+mibu)^-1);
   G3_2_c = W2n*(S2n*S2n')*W2n';
   W2m_c = mu*((2*b*G3_2_c+mu*I)^-1)*((Tm*Z2m_c+Es2m_c-Y2_2m/mu)*S2m')*((S2m*S2m'+mibu)^-1);
   G3_3_c = W3n*(S3n*S3n')*W3n';
   W3m_c = mu*((2*b*G3_3_c+mu*I)^-1)*((Tm*Z3m_c+Es3m_c-Y2_3m/mu)*S3m')*((S3m*S3m'+mibu)^-1);
   G3_4_c = W4n*(S4n*S4n')*W4n';
   W4m_c = mu*((2*b*G3_4_c+mu*I)^-1)*((Tm*Z4m_c+Es4m_c-Y2_4m/mu)*S4m')*((S4m*S4m'+mibu)^-1);
   G3_5_c = W5n*(S5n*S5n')*W5n';
   W5m_c = mu*((2*b*G3_5_c+mu*I)^-1)*((Tm*Z5m_c+Es5m_c-Y2_5m/mu)*S5m')*((S5m*S5m'+mibu)^-1);
   G3_6_c = W6n*(S6n*S6n')*W6n';
   W6m_c = mu*((2*b*G3_6_c+mu*I)^-1)*((Tm*Z6m_c+Es6m_c-Y2_6m/mu)*S6m')*((S6m*S6m'+mibu)^-1);
   
    if ~isequal(cast(W1m_c'*W1m_c,'int8'),eye(size(W1m_c))); W1m_c = orth(W1m_c); end
    if ~isequal(cast(W2m_c'*W2m_c,'int8'),eye(size(W2m_c))); W2m_c = orth(W2m_c); end
    if ~isequal(cast(W3m_c'*W3m_c,'int8'),eye(size(W3m_c))); W3m_c = orth(W3m_c); end
    if ~isequal(cast(W4m_c'*W4m_c,'int8'),eye(size(W4m_c))); W4m_c = orth(W4m_c); end
    if ~isequal(cast(W5m_c'*W5m_c,'int8'),eye(size(W5m_c))); W5m_c = orth(W5m_c); end
    if ~isequal(cast(W6m_c'*W6m_c,'int8'),eye(size(W6m_c))); W6m_c = orth(W6m_c); end
    if ~isequal(cast(W7m_c'*W7m_c,'int8'),eye(size(W7m_c))); W7m_c = orth(W7m_c); end
    % calculate constraint matrix
   dY1_1m = Z1m_c - F1m_c;dY1_2m = Z2m_c - F2m_c;dY1_3m = Z3m_c - F3m_c;dY1_4m = Z4m_c - F4m_c;dY1_5m = Z5m_c - F5m_c;dY1_6m = Z6m_c - F6m_c;
   dY1_1n = Z1n_c - F1n_c;dY1_2n = Z2n_c - F2n_c;dY1_3n = Z3n_c - F3n_c;dY1_4n = Z4n_c - F4n_c;dY1_5n = Z5n_c - F5n_c;dY1_6n = Z6n_c - F6n_c;
   dY2_1m = W1m_c*S1m-Tm*Z1m_c-Es1m_c;dY2_2m = W2m_c*S2m-Tm*Z2m_c-Es2m_c;dY2_3m = W3m_c*S3m-Tm*Z3m_c-Es3m_c;dY2_4m = W4m_c*S4m-Tm*Z4m_c-Es4m_c;dY2_5m = W5m_c*S5m-Tm*Z5m_c-Es5m_c;dY2_6m = W6m_c*S6m-Tm*Z6m_c-Es6m_c;
   dY2_1n = W1n_c*S1n-Tn*Z1n_c-Es1n_c;dY2_2n = W2n_c*S2n-Tn*Z2n_c-Es2n_c;dY2_3n = W3n_c*S3n-Tn*Z3n_c-Es3n_c;dY2_4n = W4n_c*S4n-Tn*Z4n_c-Es4n_c;dY2_5n = W5n_c*S5n-Tn*Z5n_c-Es5n_c;dY2_6n = W6n_c*S6n-Tn*Z6n_c-Es6n_c;
   
   % Calculate the infinite norm of matrices
   inf_dY1_1m = norm(dY1_1m, inf);inf_dY1_2m = norm(dY1_2m, inf);inf_dY1_3m = norm(dY1_3m, inf);inf_dY1_4m = norm(dY1_4m, inf);inf_dY1_5m = norm(dY1_5m, inf);inf_dY1_6m = norm(dY1_6m, inf);
   inf_dY1_1n = norm(dY1_1n, inf);inf_dY1_2n = norm(dY1_2n, inf);inf_dY1_3n = norm(dY1_3n, inf);inf_dY1_4n = norm(dY1_4n, inf);inf_dY1_5n = norm(dY1_5n, inf);inf_dY1_6n = norm(dY1_6n, inf);
   inf_dY2_1m = norm(dY2_1m, inf);inf_dY2_2m = norm(dY2_2m, inf);inf_dY2_3m = norm(dY2_3m, inf);inf_dY2_4m = norm(dY2_4m, inf);inf_dY2_5m = norm(dY2_5m, inf);inf_dY2_6m = norm(dY2_6m, inf);
   inf_dY2_1n = norm(dY2_1n, inf);inf_dY2_2n = norm(dY2_2n, inf);inf_dY2_3n = norm(dY2_3n, inf);inf_dY2_4n = norm(dY2_4n, inf);inf_dY2_5n = norm(dY2_5n, inf);inf_dY2_6n = norm(dY2_6n, inf);
   
   % Calculate iterations loss
   if DEBUG
       if iter == 1 || mod(iter,2) == 0
           F_norm_1 = nuclearnormF1m_c + nuclearnormF1n_c + nuclearnormF2m_c + nuclearnormF2n_c + nuclearnormF3m_c + nuclearnormF3n_c + nuclearnormF4m_c + nuclearnormF4n_c + nuclearnormF5m_c + nuclearnormF5n_c + nuclearnormF6m_c + nuclearnormF6n_c;
           Es_l1_1 = a*(norm(Es1m_c(:),1)+norm(Es2m_c(:),1)+norm(Es3m_c(:),1)+norm(Es4m_c(:),1)+norm(Es5m_c(:),1)+norm(Es6m_c(:),1)++norm(Es1n_c(:),1)+norm(Es2n_c(:),1)+norm(Es3n_c(:),1)+norm(Es4n_c(:),1)+norm(Es5n_c(:),1)+norm(Es6n_c(:),1));
           CI_term_1 = (b/2)*(norm((W1n_c*S1n)'*(W1m_c*S1m),'fro')+norm((W2n_c*S2n)'*(W2m_c*S2m),'fro')+norm((W3n_c*S3n)'*(W3m_c*S3m),'fro')+norm((W4n_c*S4n)'*(W4m_c*S4m),'fro')+norm((W5n_c*S5n)'*(W5m_c*S5m),'fro')+norm((W6n_c*S6n)'*(W6m_c*S6m),'fro'))+(b/2)*(norm((W1m_c*S1m)'*(W1n_c*S1n),'fro')+norm((W2m_c*S2m)'*(W2n_c*S2n),'fro')+norm((W3m_c*S3m)'*(W3n_c*S3n),'fro')+norm((W4m_c*S4m)'*(W4n_c*S4n),'fro')+norm((W5m_c*S5m)'*(W5n_c*S5n),'fro')+norm((W6m_c*S6m)'*(W6n_c*S6n),'fro'));
           % Lagrange multiplier subterm
           LgM_term_1 = trace(Y1_1m'*dY1_1m)+trace(Y1_2m'*dY1_2m)+trace(Y1_3m'*dY1_3m)+trace(Y1_4m'*dY1_4m)+trace(Y1_5m'*dY1_5m)+trace(Y1_6m'*dY1_6m)+trace(Y1_1n'*dY1_1n)+trace(Y1_2n'*dY1_2n)+trace(Y1_3n'*dY1_3n)+trace(Y1_4n'*dY1_4n)+trace(Y1_5n'*dY1_5n)+trace(Y1_6n'*dY1_6n)+trace(Y2_1m'*dY2_1m)+trace(Y2_2m'*dY2_2m)+trace(Y2_3m'*dY2_3m)+trace(Y2_4m'*dY2_4m)+trace(Y2_5m'*dY2_5m)+trace(Y2_6m'*dY2_6m)+trace(Y2_1n'*dY2_1n)+trace(Y2_2n'*dY2_2n)+trace(Y2_3n'*dY2_3n)+trace(Y2_4n'*dY2_4n)+trace(Y2_5n'*dY2_5n)+trace(Y2_6n'*dY2_6n);
           % penalty term
           Penalty_term_1 = (mu/2)*(norm(dY1_1m,'fro')+norm(dY1_2m,'fro')+norm(dY1_3m,'fro')+norm(dY1_4m,'fro')+norm(dY1_5m,'fro')+norm(dY1_6m,'fro')+norm(dY1_1n,'fro')+norm(dY1_2n,'fro')+norm(dY1_3n,'fro')+norm(dY1_4n,'fro')+norm(dY1_5n,'fro')+norm(dY1_6n,'fro')+norm(dY2_1m,'fro')+norm(dY2_2m,'fro')+norm(dY2_3m,'fro')+norm(dY2_4m,'fro')+norm(dY2_5m,'fro')+norm(dY2_6m,'fro')+norm(dY2_1n,'fro')+norm(dY2_2n,'fro')+norm(dY2_3n,'fro')+norm(dY2_4n,'fro')+norm(dY2_5n,'fro')+norm(dY2_6n,'fro'));
           obj_1 = F_norm_1 + Es_l1_1 + CI_term_1 + LgM_term_1 + Penalty_term_1;
           err_1 = sqrt(norm(dY1_1m,'fro')^2+norm(dY1_2m,'fro')^2+norm(dY1_3m,'fro')^2+norm(dY1_4m,'fro')^2+norm(dY1_5m,'fro')^2+norm(dY1_6m,'fro')^2+norm(dY1_1n,'fro')^2+norm(dY1_2n,'fro')^2+norm(dY1_3n,'fro')^2+norm(dY1_4n,'fro')^2+norm(dY1_5n,'fro')^2+norm(dY1_6n,'fro')^2+norm(dY2_1m,'fro')^2+norm(dY2_2m,'fro')^2+norm(dY2_3m,'fro')^2+norm(dY2_4m,'fro')^2+norm(dY2_5m,'fro')^2+norm(dY2_6m,'fro')^2+norm(dY2_1n,'fro')^2+norm(dY2_2n,'fro')^2+norm(dY2_3n,'fro')^2+norm(dY2_4n,'fro')^2+norm(dY2_5n,'fro')^2+norm(dY2_6n,'fro')^2);
             
           disp(['iter' num2str(iter) ', mu=' num2str(mu) ...
                  ', obj=' num2str(obj_1) ', err=' num2str(err_1)]);        
       end
   end
   if inf_dY1_1m < tol && inf_dY1_2m < tol && inf_dY1_3m < tol && inf_dY1_4m < tol && inf_dY1_5m < tol && inf_dY1_6m < tol && inf_dY1_7m < tol && inf_dY1_1n < tol && inf_dY1_2n < tol && inf_dY1_3n < tol && inf_dY1_4n < tol && inf_dY1_5n < tol && inf_dY1_6n < tol && inf_dY1_7n < tol && inf_dY2_1m < tol && inf_dY2_2m < tol && inf_dY2_3m < tol && inf_dY2_4m < tol && inf_dY2_5m < tol && inf_dY2_6m < tol && inf_dY2_7m < tol && inf_dY2_1n < tol && inf_dY2_2n < tol && inf_dY2_3n < tol && inf_dY2_4n < tol && inf_dY2_5n < tol && inf_dY2_6n < tol && inf_dY2_7n < tol
        break
   end
   
   Y1_1m = Y1_1m + mu*dY1_1m;Y1_2m = Y1_2m + mu*dY1_2m;Y1_3m = Y1_3m + mu*dY1_3m;Y1_4m = Y1_4m + mu*dY1_4m;Y1_5m = Y1_5m + mu*dY1_5m;Y1_6m = Y1_6m + mu*dY1_6m;
   Y1_1n = Y1_1n + mu*dY1_1n;Y1_2n = Y1_2n + mu*dY1_2n;Y1_3n = Y1_3n + mu*dY1_3n;Y1_4n = Y1_4n + mu*dY1_4n;Y1_5n = Y1_5n + mu*dY1_5n;Y1_6n = Y1_6n + mu*dY1_6n;
   Y2_1m = Y2_1m + mu*dY2_1m;Y2_2m = Y2_2m + mu*dY2_2m;Y2_3m = Y2_3m + mu*dY2_3m;Y2_4m = Y2_4m + mu*dY2_4m;Y2_5m = Y2_5m + mu*dY2_5m;Y2_6m = Y2_6m + mu*dY2_6m;
   Y2_1n = Y2_1n + mu*dY2_1n;Y2_2n = Y2_2n + mu*dY2_2n;Y2_3n = Y2_3n + mu*dY2_3n;Y2_4n = Y2_4n + mu*dY2_4n;Y2_5n = Y2_5n + mu*dY2_5n;Y2_6n = Y2_6n + mu*dY2_6n;
   mu = min(mu*rho, max_mu);
   
end
F_norm = nuclearnormF1m_c + nuclearnormF1n_c + nuclearnormF2m_c + nuclearnormF2n_c + nuclearnormF3m_c + nuclearnormF3n_c + nuclearnormF4m_c + nuclearnormF4n_c + nuclearnormF5m_c + nuclearnormF5n_c + nuclearnormF6m_c + nuclearnormF6n_c;
Es_l1 = a*(norm(Es1m_c(:),1)+norm(Es2m_c(:),1)+norm(Es3m_c(:),1)+norm(Es4m_c(:),1)+norm(Es5m_c(:),1)+norm(Es6m_c(:),1)++norm(Es1n_c(:),1)+norm(Es2n_c(:),1)+norm(Es3n_c(:),1)+norm(Es4n_c(:),1)+norm(Es5n_c(:),1)+norm(Es6n_c(:),1));
CI_term = (b/2)*(norm((W1n_c*S1n)'*(W1m_c*S1m),'fro')+norm((W2n_c*S2n)'*(W2m_c*S2m),'fro')+norm((W3n_c*S3n)'*(W3m_c*S3m),'fro')+norm((W4n_c*S4n)'*(W4m_c*S4m),'fro')+norm((W5n_c*S5n)'*(W5m_c*S5m),'fro')+norm((W6n_c*S6n)'*(W6m_c*S6m),'fro'))+(b/2)*(norm((W1m_c*S1m)'*(W1n_c*S1n),'fro')+norm((W2m_c*S2m)'*(W2n_c*S2n),'fro')+norm((W3m_c*S3m)'*(W3n_c*S3n),'fro')+norm((W4m_c*S4m)'*(W4n_c*S4n),'fro')+norm((W5m_c*S5m)'*(W5n_c*S5n),'fro')+norm((W6m_c*S6m)'*(W6n_c*S6n),'fro'));
% Lagrange multiplier subterm
LgM_term = trace(Y1_1m'*dY1_1m)+trace(Y1_2m'*dY1_2m)+trace(Y1_3m'*dY1_3m)+trace(Y1_4m'*dY1_4m)+trace(Y1_5m'*dY1_5m)+trace(Y1_6m'*dY1_6m)+trace(Y1_1n'*dY1_1n)+trace(Y1_2n'*dY1_2n)+trace(Y1_3n'*dY1_3n)+trace(Y1_4n'*dY1_4n)+trace(Y1_5n'*dY1_5n)+trace(Y1_6n'*dY1_6n)+trace(Y2_1m'*dY2_1m)+trace(Y2_2m'*dY2_2m)+trace(Y2_3m'*dY2_3m)+trace(Y2_4m'*dY2_4m)+trace(Y2_5m'*dY2_5m)+trace(Y2_6m'*dY2_6m)+trace(Y2_1n'*dY2_1n)+trace(Y2_2n'*dY2_2n)+trace(Y2_3n'*dY2_3n)+trace(Y2_4n'*dY2_4n)+trace(Y2_5n'*dY2_5n)+trace(Y2_6n'*dY2_6n);
% penalty term
Penalty_term = (mu/2)*(norm(dY1_1m,'fro')+norm(dY1_2m,'fro')+norm(dY1_3m,'fro')+norm(dY1_4m,'fro')+norm(dY1_5m,'fro')+norm(dY1_6m,'fro')+norm(dY1_1n,'fro')+norm(dY1_2n,'fro')+norm(dY1_3n,'fro')+norm(dY1_4n,'fro')+norm(dY1_5n,'fro')+norm(dY1_6n,'fro')+norm(dY2_1m,'fro')+norm(dY2_2m,'fro')+norm(dY2_3m,'fro')+norm(dY2_4m,'fro')+norm(dY2_5m,'fro')+norm(dY2_6m,'fro')+norm(dY2_1n,'fro')+norm(dY2_2n,'fro')+norm(dY2_3n,'fro')+norm(dY2_4n,'fro')+norm(dY2_5n,'fro')+norm(dY2_6n,'fro'));
obj= F_norm + Es_l1 + CI_term + LgM_term + Penalty_term;
err = sqrt(norm(dY1_1m,'fro')^2+norm(dY1_2m,'fro')^2+norm(dY1_3m,'fro')^2+norm(dY1_4m,'fro')^2+norm(dY1_5m,'fro')^2+norm(dY1_6m,'fro')^2+norm(dY1_1n,'fro')^2+norm(dY1_2n,'fro')^2+norm(dY1_3n,'fro')^2+norm(dY1_4n,'fro')^2+norm(dY1_5n,'fro')^2+norm(dY1_6n,'fro')^2+norm(dY2_1m,'fro')^2+norm(dY2_2m,'fro')^2+norm(dY2_3m,'fro')^2+norm(dY2_4m,'fro')^2+norm(dY2_5m,'fro')^2+norm(dY2_6m,'fro')^2+norm(dY2_1n,'fro')^2+norm(dY2_2n,'fro')^2+norm(dY2_3n,'fro')^2+norm(dY2_4n,'fro')^2+norm(dY2_5n,'fro')^2+norm(dY2_6n,'fro')^2);
F_sm.F1m = F1m_c;F_sm.F2m = F2m_c;F_sm.F3m = F3m_c;F_sm.F4m = F4m_c;F_sm.F5m = F5m_c;F_sm.F6m = F6m_c;
F_sn.F1n = F1n_c;F_sn.F2n = F2n_c;F_sn.F3n = F3n_c;F_sn.F4n = F4n_c;F_sn.F5n = F5n_c;F_sn.F6n = F6n_c;
Z_sm.Z1m = Z1m_c;Z_sm.Z2m = Z2m_c;Z_sm.Z3m = Z3m_c;Z_sm.Z4m = Z4m_c;Z_sm.Z5m = Z5m_c;Z_sm.Z6m = Z6m_c;
Z_sn.Z1n = Z1n_c;Z_sn.Z2n = Z2n_c;Z_sn.Z3n = Z3m_c;Z_sn.Z4n = Z4n_c;Z_sn.Z5n = Z5n_c;Z_sn.Z6n = Z6n_c;
W_sm.W1m = W1m_c;W_sm.W2m = W2m_c;W_sm.W3m = W3m_c;W_sm.W4m = W4m_c;W_sm.W5m = W5m_c;W_sm.W6m = W6m_c;
W_sn.W1n = W1n_c;W_sn.W2n = W2n_c;W_sn.W3n = W3n_c;W_sn.W4n = W4n_c;W_sn.W5n = W5n_c;W_sn.W6n = W6n_c;
Es_sm.Es1m = Es1m_c;Es_sm.Es2m = Es2m_c;Es_sm.Es3m = Es3m_c;Es_sm.Es4m = Es4m_c;Es_sm.Es5m = Es5m_c;Es_sm.Es6m = Es6m_c;
Es_sn.Es1n = Es1n_c;Es_sn.Es2n = Es2n_c;Es_sn.Es3n = Es3n_c;Es_sn.Es4n = Es4n_c;Es_sn.Es5n = Es5n_c;Es_sn.Es6n = Es6n_c;               
end