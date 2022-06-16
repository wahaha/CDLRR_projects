clc
load Synthetic_Data.txt;    % 加载人工数据
data = Synthetic_Data;

% Obtain Target domain datas
%  class 1
Tm1 = data(1,1:200);Tm2 = data(2,1:200);
Tm = [Tm1;Tm2];
%  class 2
Tn1 = data(1,201:400);Tn2 = data(2,201:400);
Tn = [Tn1;Tn2];
% T = [Tm Tn];

% Obtain Source1 domain datas
%  class 1
Sm1_1 = data(3,1:200);Sm1_2 = data(4,1:200);
Sm1 = [Sm1_1;Sm1_2];
S1m = Sm1;
%  class 2
Sn1_1 = data(3,201:400);Sn1_2 = data(4,201:400);
Sn1 = [Sn1_1;Sn1_2];
S1n = Sn1;

% Obtain Source2 domain datas
%  class 1
Sm2_1 = data(5,1:200);Sm2_2 = data(6,1:200);
Sm2 = [Sm2_1;Sm2_2];
S2m = Sm2;
%  class 2
Sn2_1 = data(5,201:400);Sn2_2 = data(6,201:400);
Sn2 = [Sn2_1;Sn2_2];
S2n = Sn2;

Ssm.S1m = S1m;Ssm.S2m = S2m;
Ssn.S1n = S1n;Ssn.S2n = S2n;

opts.tol = 1e-7;
opts.max_iter = 120;
opts.max_mu = 1e7;
opts.DEBUG = 1;

% 给定寻参范围
a = optimizableVariable('a', [0 0.1], 'Type', 'real');
b = optimizableVariable('b', [0 0.1], 'Type', 'real');
rho = optimizableVariable('rho', [1 1.2], 'Type', 'real');
mu = optimizableVariable('mu', [1e-7 1e-5], 'Type', 'real');


parameter = [a, b, rho, mu];
objFun = @(parameter) getObjValue(parameter, Tm, Tn, Ssm, Ssn, opts);

iters = 60;
points = 10;
ratio  = 0.5;
% acqfn = 'expected-improvement-plus';
acqfn = 'expected-improvement-per-second-plus';
% results = bayesopt(objFun, parameter,'Verbose',1,'ExplorationRatio',ratio,...
%                    'NumCoupledConstraints',1,'MaxObjectiveEvaluations',iters,...
%                    'AcquisitionFunctionName',acqfn,'NumSeedPoints',points);
results = bayesopt(objFun, parameter,'Verbose',1,'ExplorationRatio',ratio,...
                   'MaxObjectiveEvaluations',iters,...
                   'AcquisitionFunctionName',acqfn,'NumSeedPoints',points);
               
function [obj] = getObjValue(parameter, Tm, Tn, Ssm, Ssn, opts)
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
if isfield(opts, 'DEBUG');      DEBUG = opts.DEBUG;        end

a = parameter.a;
b = parameter.b;
rho = parameter.rho;
mu = parameter.mu;

[d,nTm] = size(Tm);
[~,nTn] = size(Tn);

S1m = Ssm.S1m;
[~,nS1m] = size(S1m);
S1n = Ssn.S1n;
[~,nS1n] = size(S1n);

S2m = Ssm.S2m;
[~,nS2m] = size(S2m);
S2n = Ssn.S2n;
[~,nS2n] = size(S2n);

I = eye(d);

mibu = 1e-7*I;

Z1m = zeros(nTm,nS1m);Z1n = zeros(nTn,nS1n);
Z2m = zeros(nTm,nS2m);Z2n = zeros(nTn,nS2n);
Z1m_c = Z1m;Z1n_c = Z1n;
Z2m_c = Z2m;Z2n_c = Z2n;

F1m = zeros(nTm,nS1m);F1n = zeros(nTn,nS1n);
F2m = zeros(nTm,nS2m);F2n = zeros(nTn,nS2n);
F1m_c = F1m;F1n_c = F1n;
F2m_c = F2m;F2n_c = F2n;

Z_I_m = eye(nTm);
Z_I_n = eye(nTn);

E1m = zeros(d,nS1m);E1n = zeros(d,nS1n);
E2m = zeros(d,nS2m);E2n = zeros(d,nS2n);
E1m_c = E1m;E1n_c = E1n;
E2m_c = E2m;E2n_c = E2n;

P1m = I;P1n = I;
P2m = I;P2n = I;
P1m_c = P1m;P1n_c = P1n;
P2m_c = P2m;P2n_c = P2n;

P = I;

% lagrangian multiplier
Y1_1m = zeros(nTm,nS1m);Y1_1n = zeros(nTn,nS1n);Y2_1m = zeros(d,nS1m);Y2_1n = zeros(d,nS1n);
Y1_2m = zeros(nTm,nS2m);Y1_2n = zeros(nTn,nS2n);Y2_2m = zeros(d,nS2m);Y2_2n = zeros(d,nS2n);

for iter = 1:max_iter
    % Update Pim
        % Update Fin_c
    [F1n_c,~] = prox_nuclear(Z1n_c+Y1_1n/mu,1/mu);
    [F2n_c,~] = prox_nuclear(Z2n_c+Y1_2n/mu,1/mu);
        % Update Fim
    [F1m,nuclearnormF1m] = prox_nuclear(Z1m+Y1_1m/mu,1/mu);
    [F2m,nuclearnormF2m] = prox_nuclear(Z2m+Y1_2m/mu,1/mu);
        % Update Ein_c
    E1n_c = prox_l1(P1n_c*S1n-P*Tn*Z1n_c+Y2_1n/mu,a/mu);
    E2n_c = prox_l1(P2n_c*S2n-P*Tn*Z2n_c+Y2_2n/mu,a/mu);
        % Update Eim
    E1m = prox_l1(P1m*S1m-P*Tm*Z1m+Y2_1m/mu,a/mu);
    E2m = prox_l1(P2m*S2m-P*Tm*Z2m+Y2_2m/mu,a/mu);
        % Update Zin_c
    G2_1_c = P1n_c*S1n-E1n_c+Y2_1n/mu;
    Z1n_c = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_1_c+F1n_c-Y1_1n/mu);
    G2_2_c = P2n_c*S2n-E2n_c+Y2_2n/mu;
    Z2n_c = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_2_c+F2n_c-Y1_2n/mu);
        % Update Zim
    G1_1 = P1m*S1m - E1m + Y2_1m/mu;
    Z1m = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_1+F1m-Y1_1m/mu);
    G1_2 = P2m*S2m - E2m + Y2_2m/mu;
    Z2m = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_2+F2m-Y1_2m/mu);
        % Update Pin_c
    G5_1_c = P1m*(S1m*S1m')*P1m';G6_1_c = P*Tn*Z1n_c+E1n_c-Y2_1n/mu;
    P1n_c = ((2*b*G5_1_c+mu*I)^-1)*mu*G6_1_c*S1n'*((S1n*S1n')^-1);
    G5_2_c = P2m*(S2m*S2m')*P2m';G6_2_c = P*Tn*Z2n_c+E2n_c-Y2_2n/mu;
    P2n_c = ((2*b*G5_2_c+mu*I)^-1)*mu*G6_2_c*S2n'*((S2n*S2n')^-1);
        % Update Pim
    G3_1 = P1n_c*(S1n*S1n')*P1n_c';G4_1 = P*Tm*Z1m+E1m-Y2_1m/mu;
    P1m = ((2*b*G3_1+mu*I)^-1)*mu*G4_1*S1m'*((S1m*S1m')^-1);
    G3_2 = P2n_c*(S2n*S2n')*P2n_c';G4_2 = P*Tm*Z2m+E2m-Y2_2m/mu;
    P2m = ((2*b*G3_2+mu*I)^-1)*mu*G4_2*S2m'*((S2m*S2m')^-1);
    
    % 正交约束 Pim
    if cast(P1m'*P1m, 'int8') ~= eye(size(P1m)); P1m = orth(P1m); end
    if cast(P2m'*P2m, 'int8') ~= eye(size(P2m)); P2m = orth(P2m); end
    
    % Update Pin
        % Update Fim_c
    [F1m_c,~] = prox_nuclear(Z1m_c+Y1_1m/mu,1/mu);
    [F2m_c,~] = prox_nuclear(Z2m_c+Y1_2m/mu,1/mu);
        % Update Fin
    [F1n,nuclearnormF1n] = prox_nuclear(Z1n+Y1_1n/mu,1/mu);
    [F2n,nuclearnormF2n] = prox_nuclear(Z2n+Y1_2n/mu,1/mu);
        % Update Eim_c
    E1m_c = prox_l1(P1m_c*S1m-P*Tm*Z1m_c+Y2_1m/mu,a/mu);
    E2m_c = prox_l1(P2m_c*S2m-P*Tm*Z2m_c+Y2_2m/mu,a/mu);
        % Update Ein
    E1n = prox_l1(P1n*S1n-P*Tn*Z1n+Y2_1n/mu,a/mu);
    E2n = prox_l1(P2n*S2n-P*Tn*Z2n+Y2_2n/mu,a/mu);
        % Update Zim_c
    G1_1_c = P1m_c*S1m - E1m_c + Y2_1m/mu;
    Z1m_c = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_1_c+F1m_c-Y1_1m/mu);
    G1_2_c = P2m_c*S2m - E2m_c + Y2_2m/mu;
    Z2m_c = ((Tm'*(P'*P)*Tm+Z_I_m)^-1)*(Tm'*P'*G1_2_c+F2m_c-Y1_2m/mu);
        % Update Zin
    G2_1 = P1n*S1n-E1n+Y2_1n/mu;
    Z1n = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_1+F1n-Y1_1n/mu);
    G2_2 = P2n*S2n-E2n+Y2_2n/mu;
    Z2n = ((Tn'*(P'*P)*Tn+Z_I_n)^-1)*(Tn'*P'*G2_2+F2n-Y1_2n/mu);
        % Update Pim_c
    G3_1_c = P1n*(S1n*S1n')*P1n';G4_1_c = P*Tm*Z1m_c+E1m_c-Y2_1m/mu;
    P1m_c = ((2*b*G3_1_c+mu*I)^-1)*mu*G4_1_c*S1m'*((S1m*S1m')^-1);
    G3_2_c = P2n*(S2n*S2n')*P2n';G4_2_c = P*Tm*Z2m_c+E2m_c-Y2_2m/mu;
    P2m_c = ((2*b*G3_2_c+mu*I)^-1)*mu*G4_2_c*S2m'*((S2m*S2m')^-1);
        % Update Pin
    G5_1 = P1m_c*(S1m*S1m')*P1m_c';G6_1 = P*Tn*Z1n+E1n-Y2_1n/mu;
    P1n = ((2*b*G5_1+mu*I)^-1)*mu*G6_1*S1n'*((S1n*S1n')^-1);
    G5_2 = P2m_c*(S2m*S2m')*P2m_c';G6_2 = P*Tn*Z2n+E2n-Y2_2n/mu;
    P2n = ((2*b*G5_2+mu*I)^-1)*mu*G6_2*S2n'*((S2n*S2n')^-1);
    
    % 正交约束 Pin
    if cast(P1n'*P1n, 'int8') ~= eye(size(P1n)); P1n = orth(P1n); end
    if cast(P2n'*P2n, 'int8') ~= eye(size(P2n)); P2n = orth(P2n); end
    
    % P s.t. P'P=I Orithogonal constraint
    if cast(P'*P, 'int8') ~= eye(size(P)); P = orth(P); end
    
    % calculate constraint matrix
    dY1_1m = Z1m - F1m;dY1_1n = Z1n - F1n;dY2_1m = P1m*S1m-P*Tm*Z1m-E1m;dY2_1n = P1n*S1n-P*Tn*Z1n-E1n;
    dY1_2m = Z2m - F2m;dY1_2n = Z2n - F2n;dY2_2m = P2m*S2m-P*Tm*Z2m-E2m;dY2_2n = P2n*S2n-P*Tn*Z2n-E2n;
    
    % Calculate the infinite norm of matrices
    inf_dY1_1m = norm(dY1_1m,inf);inf_dY1_1n = norm(dY1_1n,inf);inf_dY2_1m = norm(dY2_1m,inf);inf_dY2_1n = norm(dY2_1n,inf);
    inf_dY1_2m = norm(dY1_2m,inf);inf_dY1_2n = norm(dY1_2n,inf);inf_dY2_2m = norm(dY2_2m,inf);inf_dY2_2n = norm(dY2_2n,inf);
    
    % Calculate iterations loss
    if DEBUG
        if iter == 1 || mod(iter,2) == 0
            F_norm_1 = nuclearnormF1m + nuclearnormF1n + nuclearnormF2m + nuclearnormF2n;
            Es_l1_1 = a*(norm(E1m(:),1)+norm(E1n(:),1)+norm(E2m(:),1)+norm(E2n(:),1));
            CI_term_1 = (b/2)*(norm((P1n*S1n)'*(P1m*S1m),'fro')+norm((P1m*S1m)'*(P1n*S1n),'fro')+norm((P2n*S2n)'*(P2m*S2m),'fro')+norm((P2m*S2m)'*(P2n*S2n),'fro'));
            % Lagrange multiplier subterm
            LgM_term_1 = trace(Y1_1m'*dY1_1m)+trace(Y1_1n'*dY1_1n)+trace(Y2_1m'*dY2_1m)+trace(Y2_1n'*dY2_1n)+trace(Y1_2m'*dY1_2m)+trace(Y1_2n'*dY1_2n)+trace(Y2_2m'*dY2_2m)+trace(Y2_2n'*dY2_2n);
            % penalty term
            Penalty_term_1 = (mu/2)*(norm(dY1_1m,'fro')^2+norm(dY1_1n,'fro')^2+norm(dY2_1m,'fro')^2+norm(dY2_1n,'fro')^2+norm(dY1_2m,'fro')^2+norm(dY1_2n,'fro')^2+norm(dY2_2m,'fro')^2+norm(dY2_2n,'fro')^2);
            obj_1 = F_norm_1 + Es_l1_1 + CI_term_1 + LgM_term_1 + Penalty_term_1;
            err_1 = sqrt(norm(dY1_1m,'fro')^2+norm(dY1_1n,'fro')^2+norm(dY2_1m,'fro')^2+norm(dY2_1n,'fro')^2+norm(dY1_2m,'fro')^2+norm(dY1_2n,'fro')^2+norm(dY2_2m,'fro')^2+norm(dY2_2n,'fro')^2);
            
            disp(['iter' num2str(iter) ', mu=' num2str(mu) ...
                  ', obj=' num2str(obj_1) ', err=' num2str(err_1)]);  
        end
    end
    
    % 判断是否收敛
    if inf_dY1_1m < tol && inf_dY1_1n < tol && inf_dY2_1m < tol && inf_dY2_1n < to1 && inf_dY1_2m < tol && inf_dY1_2n < tol && inf_dY2_2m < tol && inf_dY2_2n < tol
        break;
    end
    
    % Update the multipliers and parameters
    Y1_1m = Y1_1m + mu*dY1_1m;Y1_1n = Y1_1n + mu*dY1_1n;Y2_1m = Y2_1m + mu*dY2_1m;Y2_1n = Y2_1n + mu*dY2_1n;
    Y1_2m = Y1_2m + mu*dY1_2m;Y1_2n = Y1_2n + mu*dY1_2n;Y2_2m = Y2_2m + mu*dY2_2m;Y2_2n = Y2_2n + mu*dY2_2n;
    
    mu = min(mu*rho, max_mu);
end
F_norm = nuclearnormF1m + nuclearnormF1n + nuclearnormF2m + nuclearnormF2n;
Es_l1 = a*(norm(E1m(:),1)+norm(E1n(:),1)+norm(E2m(:),1)+norm(E2n(:),1));
CI_term = (b/2)*(norm((P1n*S1n)'*(P1m*S1m),'fro')+norm((P1m*S1m)'*(P1n*S1n),'fro')+norm((P2n*S2n)'*(P2m*S2m),'fro')+norm((P2m*S2m)'*(P2n*S2n),'fro'));
% Lagrange multiplier subterm
LgM_term = trace(Y1_1m'*dY1_1m)+trace(Y1_1n'*dY1_1n)+trace(Y2_1m'*dY2_1m)+trace(Y2_1n'*dY2_1n)+trace(Y1_2m'*dY1_2m)+trace(Y1_2n'*dY1_2n)+trace(Y2_2m'*dY2_2m)+trace(Y2_2n'*dY2_2n);
% penalty term
Penalty_term = (mu/2)*(norm(dY1_1m,'fro')^2+norm(dY1_1n,'fro')^2+norm(dY2_1m,'fro')^2+norm(dY2_1n,'fro')^2+norm(dY1_2m,'fro')^2+norm(dY1_2n,'fro')^2+norm(dY2_2m,'fro')^2+norm(dY2_2n,'fro')^2);
obj = F_norm + Es_l1 + CI_term + LgM_term + Penalty_term;
err = sqrt(norm(dY1_1m,'fro')^2+norm(dY1_1n,'fro')^2+norm(dY2_1m,'fro')^2+norm(dY2_1n,'fro')^2+norm(dY1_2m,'fro')^2+norm(dY1_2n,'fro')^2+norm(dY2_2m,'fro')^2+norm(dY2_2n,'fro')^2);
            
end

