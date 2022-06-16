function [FS,ZS,WS,ES,obj,err] = CDLRR_1(XT,XS,a,opts)
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

[d,nT] = size(XT);

Xs1 = XS.Xs1;
[~,nS1] = size(Xs1);
Xs2 = XS.Xs2;
[~,nS2] = size(Xs2);
Xs3 = XS.Xs3;
[~,nS3] = size(Xs3);
Xs4 = XS.Xs4;
[~,nS4] = size(Xs4);
Xs5 = XS.Xs5;
[~,nS5] = size(Xs5);
Xs6 = XS.Xs6;
[~,nS6] = size(Xs6);

I = eye(d);
Z_I = eye(nT);
mibu = 1e-7*I;

W1 = I;W2 = I;W3 = I;
W4 = I;W5 = I;W6 = I;

Z1 = zeros(nT,nS1);
Z2 = zeros(nT,nS2);
Z3 = zeros(nT,nS3);
Z4 = zeros(nT,nS4);
Z5 = zeros(nT,nS5);
Z6 = zeros(nT,nS6);

Es1 = zeros(d,nS1);Es2 = zeros(d,nS2);
Es3 = zeros(d,nS3);Es4 = zeros(d,nS4);
Es5 = zeros(d,nS5);Es6 = zeros(d,nS6);

% lagrangian multiplier
Y1_1 = zeros(nT,nS1);Y2_1 = zeros(d,nS1);
Y1_2 = zeros(nT,nS2);Y2_2 = zeros(d,nS2);
Y1_3 = zeros(nT,nS3);Y2_3 = zeros(d,nS3);
Y1_4 = zeros(nT,nS4);Y2_4 = zeros(d,nS4);
Y1_5 = zeros(nT,nS5);Y2_5 = zeros(d,nS5);
Y1_6 = zeros(nT,nS6);Y2_6 = zeros(d,nS6);

for iter = 1:max_iter
    % update Fi
    [F1,nuclearnormF1] = prox_nuclear(Z1+Y1_1/mu, 1/mu);
    [F2,nuclearnormF2] = prox_nuclear(Z2+Y1_2/mu, 1/mu);
    [F3,nuclearnormF3] = prox_nuclear(Z3+Y1_3/mu, 1/mu);
    [F4,nuclearnormF4] = prox_nuclear(Z4+Y1_4/mu, 1/mu);
    [F5,nuclearnormF5] = prox_nuclear(Z5+Y1_5/mu, 1/mu);
    [F6,nuclearnormF6] = prox_nuclear(Z6+Y1_6/mu, 1/mu);
    
    % update Esi
    Es1 = prox_l1(W1*Xs1-XT*Z1+Y2_1/mu,a/mu);
    Es2 = prox_l1(W2*Xs2-XT*Z2+Y2_2/mu,a/mu);
    Es3 = prox_l1(W3*Xs3-XT*Z3+Y2_3/mu,a/mu);
    Es4 = prox_l1(W4*Xs4-XT*Z4+Y2_4/mu,a/mu);
    Es5 = prox_l1(W5*Xs5-XT*Z5+Y2_5/mu,a/mu);
    Es6 = prox_l1(W6*Xs6-XT*Z6+Y2_6/mu,a/mu);
    
    % update Zi
    G1_1 = W1*Xs1-Es1+Y2_1/mu;
    Z1 = ((XT'*XT+Z_I)^-1)*(XT'*G1_1+F1-Y1_1/mu);
    G1_2 = W2*Xs2-Es2+Y2_2/mu;
    Z2 = ((XT'*XT+Z_I)^-1)*(XT'*G1_2+F2-Y1_2/mu);
    G1_3 = W3*Xs3-Es3+Y2_3/mu;
    Z3 = ((XT'*XT+Z_I)^-1)*(XT'*G1_3+F3-Y1_3/mu);
    G1_4 = W4*Xs4-Es4+Y2_4/mu;
    Z4 = ((XT'*XT+Z_I)^-1)*(XT'*G1_4+F4-Y1_4/mu);
    G1_5 = W5*Xs5-Es5+Y2_5/mu;
    Z5 = ((XT'*XT+Z_I)^-1)*(XT'*G1_5+F5-Y1_5/mu);
    G1_6 = W6*Xs6-Es6+Y2_6/mu;
    Z6 = ((XT'*XT+Z_I)^-1)*(XT'*G1_6+F6-Y1_6/mu);
    
    % update Wi
    W1 = (XT*Z1+Es1-Y2_1/mu)*Xs1'*((Xs1*Xs1'+mibu)^-1);
    W2 = (XT*Z2+Es2-Y2_2/mu)*Xs2'*((Xs2*Xs2'+mibu)^-1);
    W3 = (XT*Z3+Es3-Y2_3/mu)*Xs3'*((Xs3*Xs3'+mibu)^-1);
    W4 = (XT*Z4+Es4-Y2_4/mu)*Xs4'*((Xs4*Xs4'+mibu)^-1);
    W5 = (XT*Z5+Es5-Y2_5/mu)*Xs5'*((Xs5*Xs5'+mibu)^-1);
    W6 = (XT*Z6+Es6-Y2_6/mu)*Xs6'*((Xs6*Xs6'+mibu)^-1);
    
    % Wi s.t. Wi'Wi=I Orthogonal constraint
    if cast(W1'*W1, 'int8') ~= eye(size(W1)); W1 = orth(W1); end
    if cast(W2'*W2, 'int8') ~= eye(size(W2)); W2 = orth(W2); end
    if cast(W3'*W3, 'int8') ~= eye(size(W3)); W3 = orth(W3); end
    if cast(W4'*W4, 'int8') ~= eye(size(W4)); W4 = orth(W4); end
    if cast(W5'*W5, 'int8') ~= eye(size(W5)); W5 = orth(W5); end
    if cast(W6'*W6, 'int8') ~= eye(size(W6)); W6 = orth(W6); end
    
    % calculate constraint matrix
    dY1_1 = Z1-F1;dY2_1 = W1*Xs1-XT*Z1-Es1;
    dY1_2 = Z2-F2;dY2_2 = W2*Xs2-XT*Z2-Es2;
    dY1_3 = Z3-F3;dY2_3 = W3*Xs3-XT*Z3-Es3;
    dY1_4 = Z4-F4;dY2_4 = W4*Xs4-XT*Z4-Es4;
    dY1_5 = Z5-F5;dY2_5 = W5*Xs5-XT*Z5-Es5;
    dY1_6 = Z6-F6;dY2_6 = W6*Xs6-XT*Z6-Es6;
    
    % Calculate the infinite norm of matrices
    inf_dY1_1 = norm(dY1_1, inf);inf_dY2_1 = norm(dY2_1, inf);
    inf_dY1_2 = norm(dY1_2, inf);inf_dY2_2 = norm(dY2_2, inf);
    inf_dY1_3 = norm(dY1_3, inf);inf_dY2_3 = norm(dY2_3, inf);
    inf_dY1_4 = norm(dY1_4, inf);inf_dY2_4 = norm(dY2_4, inf);
    inf_dY1_5 = norm(dY1_5, inf);inf_dY2_5 = norm(dY2_5, inf);
    inf_dY1_6 = norm(dY1_6, inf);inf_dY2_6 = norm(dY2_6, inf);
    
    % Calculate the Fro2 norm of matrices
    fro2_dY1_1 = norm(dY1_1,'fro')^2;fro2_dY2_1 = norm(dY2_1,'fro')^2;
    fro2_dY1_2 = norm(dY1_2,'fro')^2;fro2_dY2_2 = norm(dY2_2,'fro')^2;
    fro2_dY1_3 = norm(dY1_3,'fro')^2;fro2_dY2_3 = norm(dY2_3,'fro')^2;
    fro2_dY1_4 = norm(dY1_4,'fro')^2;fro2_dY2_4 = norm(dY2_4,'fro')^2;
    fro2_dY1_5 = norm(dY1_5,'fro')^2;fro2_dY2_5 = norm(dY2_5,'fro')^2;
    fro2_dY1_6 = norm(dY1_6,'fro')^2;fro2_dY2_6 = norm(dY2_6,'fro')^2;
    
    if DEBUG
        if iter == 1 || mod(iter,2) == 0
            Fi_norm1 = nuclearnormF1 + nuclearnormF2 + nuclearnormF3 + nuclearnormF4 + nuclearnormF5 + nuclearnormF6;
            Esi_l1_1 = a*(norm(Es1(:),1) + norm(Es2(:),1) + norm(Es3(:),1) + norm(Es4(:),1) + norm(Es5(:),1) + norm(Es6(:),1));
            % Lagrange multiplier subterm
            lgm_term_1 = trace(Y1_1'*dY1_1)+trace(Y1_2'*dY1_2)+trace(Y1_3'*dY1_3)+trace(Y1_4'*dY1_4)+trace(Y1_5'*dY1_5)+trace(Y1_6'*dY1_6);
            lgm_term_2 = trace(Y2_1'*dY2_1)+trace(Y2_2'*dY2_2)+trace(Y2_3'*dY2_3)+trace(Y2_4'*dY2_4)+trace(Y2_5'*dY2_5)+trace(Y2_6'*dY2_6);
            LgM_term_sum1 = lgm_term_1 + lgm_term_2;
            % penalty term
            penalty_term1 = fro2_dY1_1+fro2_dY1_2+fro2_dY1_3+fro2_dY1_4+fro2_dY1_5+fro2_dY1_6;
            penalty_term2 = fro2_dY2_1+fro2_dY2_2+fro2_dY2_3+fro2_dY2_4+fro2_dY2_5+fro2_dY2_6;
            Penalty_term_1 = (mu/2)*(penalty_term1+penalty_term2);
            obj_1 = Fi_norm1 + Esi_l1_1 + LgM_term_sum1 + Penalty_term_1;
            err_1 = sqrt(penalty_term1 + penalty_term2);
            disp(['iter' num2str(iter) ', mu=' num2str(mu) ...
                   ', obj=' num2str(obj_1) ', err=' num2str(err_1)]);  
        end
    end
     % Check the convergence conditions
    if inf_dY1_1 < tol && inf_dY1_2 < tol && inf_dY1_3 < tol && inf_dY1_4 < tol && inf_dY1_5 < tol && inf_dY1_6 < tol && inf_dY2_1 < tol && inf_dY2_2 < tol && inf_dY2_3 < tol && inf_dY2_4 < tol && inf_dY2_5 < tol && inf_dY2_6 < tol
       break;
    end    
    % Update the multipliers and parameters
    Y1_1 = Y1_1 + mu*dY1_1;Y2_1 = Y2_1 + mu*dY2_1;
    Y1_2 = Y1_2 + mu*dY1_2;Y2_2 = Y2_2 + mu*dY2_2;
    Y1_3 = Y1_3 + mu*dY1_3;Y2_3 = Y2_3 + mu*dY2_3;
    Y1_4 = Y1_4 + mu*dY1_4;Y2_4 = Y2_4 + mu*dY2_4;
    Y1_5 = Y1_5 + mu*dY1_5;Y2_5 = Y2_5 + mu*dY2_5;
    Y1_6 = Y1_6 + mu*dY1_6;Y2_6 = Y2_6 + mu*dY2_6;
    mu = min(mu*rho, max_mu); 
end
Fi_norm = nuclearnormF1 + nuclearnormF2 + nuclearnormF3 + nuclearnormF4 + nuclearnormF5 + nuclearnormF6 ;
Esi_l1 = a*(norm(Es1(:),1) + norm(Es2(:),1) + norm(Es3(:),1) + norm(Es4(:),1) + norm(Es5(:),1) + norm(Es6(:),1));
% Lagrange multiplier subterm
lgm_term_1 = trace(Y1_1'*dY1_1)+trace(Y1_2'*dY1_2)+trace(Y1_3'*dY1_3)+trace(Y1_4'*dY1_4)+trace(Y1_5'*dY1_5)+trace(Y1_6'*dY1_6);
lgm_term_2 = trace(Y2_1'*dY2_1)+trace(Y2_2'*dY2_2)+trace(Y2_3'*dY2_3)+trace(Y2_4'*dY2_4)+trace(Y2_5'*dY2_5)+trace(Y2_6'*dY2_6);
LgM_term_sum = lgm_term_1 + lgm_term_2;
% penalty term
penalty_term1 = fro2_dY1_1+fro2_dY1_2+fro2_dY1_3+fro2_dY1_4+fro2_dY1_5+fro2_dY1_6;
penalty_term2 = fro2_dY2_1+fro2_dY2_2+fro2_dY2_3+fro2_dY2_4+fro2_dY2_5+fro2_dY2_6;
Penalty_term = (mu/2)*(penalty_term1+penalty_term2);
obj = Fi_norm + Esi_l1 + LgM_term_sum + Penalty_term;
err = sqrt(penalty_term1 + penalty_term2);
FS.F1 = F1;FS.F2 = F2;FS.F3 = F3;FS.F4 = F4;FS.F5 = F5;FS.F6 = F6;
ZS.Z1 = Z1;ZS.Z2 = Z2;ZS.Z3 = Z3;ZS.Z4 = Z4;ZS.Z5 = Z5;ZS.Z6 = Z6;
WS.W1 = W1;WS.W2 = W2;WS.W3 = W3;WS.W4 = W4;WS.W5 = W5;WS.W6 = W6;
ES.Es1 = Es1;ES.Es2 = Es2;ES.Es3 = Es3;ES.Es4 = Es4;ES.Es5 = Es5;ES.Es6 = Es6;
end