function [FS,ZS,PS,P,ES,obj,err] = CDLRR_2(T,XS,a,opts)
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

[d,nT] = size(T);

S1 = XS.S1;
[~,nS1] = size(S1);
S2 = XS.S2;
[~,nS2] = size(S2);
S3 = XS.S3;
[~,nS3] = size(S3);
S4 = XS.S4;
[~,nS4] = size(S4);
S5 = XS.S5;
[~,nS5] = size(S5);
S6 = XS.S6;
[~,nS6] = size(S6);

I = eye(d);
Z_I = eye(nT);
mibu = 1e-7*I;

P1 = I;P2 = I;P3 = I;
P4 = I;P5 = I;P6 = I;
P = I;

Z1 = zeros(nT,nS1);
Z2 = zeros(nT,nS2);
Z3 = zeros(nT,nS3);
Z4 = zeros(nT,nS4);
Z5 = zeros(nT,nS5);
Z6 = zeros(nT,nS6);

E1 = zeros(d,nS1);E2 = zeros(d,nS2);
E3 = zeros(d,nS3);E4 = zeros(d,nS4);
E5 = zeros(d,nS5);E6 = zeros(d,nS6);

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
    
    % update Ei
    E1 = prox_l1(P1*S1-P*T*Z1+Y2_1/mu,a/mu);
    E2 = prox_l1(P2*S2-P*T*Z2+Y2_2/mu,a/mu);
    E3 = prox_l1(P3*S3-P*T*Z3+Y2_3/mu,a/mu);
    E4 = prox_l1(P4*S4-P*T*Z4+Y2_4/mu,a/mu);
    E5 = prox_l1(P5*S5-P*T*Z5+Y2_5/mu,a/mu);
    E6 = prox_l1(P6*S6-P*T*Z6+Y2_6/mu,a/mu);
    
    % update Zi
    G1_1 = P1*S1-E1+Y2_1/mu;
    Z1 = ((T'*(P'*P)*T+Z_I)^-1)*(T'*P'*G1_1+F1-Y1_1/mu);
    G1_2 = P2*S2-E2+Y2_2/mu;
    Z2 = ((T'*(P'*P)*T+Z_I)^-1)*(T'*P'*G1_2+F2-Y1_2/mu);
    G1_3 = P3*S3-E3+Y2_3/mu;
    Z3 = ((T'*(P'*P)*T+Z_I)^-1)*(T'*P'*G1_3+F3-Y1_3/mu);
    G1_4 = P4*S4-E4+Y2_4/mu;
    Z4 = ((T'*(P'*P)*T+Z_I)^-1)*(T'*P'*G1_4+F4-Y1_4/mu);
    G1_5 = P5*S5-E5+Y2_5/mu;
    Z5 = ((T'*(P'*P)*T+Z_I)^-1)*(T'*P'*G1_5+F5-Y1_5/mu);
    G1_6 = P6*S6-E6+Y2_6/mu;
    Z6 = ((T'*(P'*P)*T+Z_I)^-1)*(T'*P'*G1_6+F6-Y1_6/mu);
    
    % update Pi
    P1 = (P*T*Z1+E1-Y2_1/mu)*S1'*((S1*S1'+mibu)^-1);
    P2 = (P*T*Z2+E2-Y2_2/mu)*S2'*((S2*S2'+mibu)^-1);
    P3 = (P*T*Z3+E3-Y2_3/mu)*S3'*((S3*S3'+mibu)^-1);
    P4 = (P*T*Z4+E4-Y2_4/mu)*S4'*((S4*S4'+mibu)^-1);
    P5 = (P*T*Z5+E5-Y2_5/mu)*S5'*((S5*S5'+mibu)^-1);
    P6 = (P*T*Z6+E6-Y2_6/mu)*S6'*((S6*S6'+mibu)^-1);
    
    % Pi s.t. Pi'Pi=I Orthogonal constraint
    if cast(P1'*P1,'int8') ~= eye(size(P1)); P1 = orth(P1); end
    if cast(P2'*P2,'int8') ~= eye(size(P2)); P2 = orth(P2); end
    if cast(P3'*P3,'int8') ~= eye(size(P3)); P3 = orth(P3); end
    if cast(P4'*P4,'int8') ~= eye(size(P4)); P4 = orth(P4); end
    if cast(P5'*P5,'int8') ~= eye(size(P5)); P5 = orth(P5); end
    if cast(P6'*P6,'int8') ~= eye(size(P6)); P6 = orth(P6); end
    
    % update P
    G2_1 = (P1*S1-E1+Y2_1/mu)*Z1'*T';
    G2_2 = (P2*S2-E2+Y2_2/mu)*Z2'*T';
    G2_3 = (P3*S3-E3+Y2_3/mu)*Z3'*T';
    G2_4 = (P4*S4-E4+Y2_4/mu)*Z4'*T';
    G2_5 = (P5*S5-E5+Y2_5/mu)*Z5'*T';
    G2_6 = (P6*S6-E6+Y2_6/mu)*Z6'*T';
    G2 = G2_1 + G2_2 + G2_3 + G2_4 + G2_5 + G2_6;
    
    G3_1 = T*(Z1*Z1')*T';
    G3_2 = T*(Z2*Z2')*T';
    G3_3 = T*(Z3*Z3')*T';
    G3_4 = T*(Z4*Z4')*T';
    G3_5 = T*(Z5*Z5')*T';
    G3_6 = T*(Z6*Z6')*T';
    G3 = G3_1 + G3_2 + G3_3 + G3_4 + G3_5 + G3_6;
    P = G2*((G3+mibu)^-1);
    
    % calculate constraint matrix
    dY1_1 = Z1-F1;dY2_1 = P1*S1-P*T*Z1-E1;
    dY1_2 = Z2-F2;dY2_2 = P2*S2-P*T*Z2-E2;
    dY1_3 = Z3-F3;dY2_3 = P3*S3-P*T*Z3-E3;
    dY1_4 = Z4-F4;dY2_4 = P4*S4-P*T*Z4-E4;
    dY1_5 = Z5-F5;dY2_5 = P5*S5-P*T*Z5-E5;
    dY1_6 = Z6-F6;dY2_6 = P6*S6-P*T*Z6-E6;
    
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
            Esi_l1_1 = a*(norm(E1(:),1) + norm(E2(:),1) + norm(E3(:),1) + norm(E4(:),1) + norm(E5(:),1) + norm(E6(:),1));
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
Esi_l1 = a*(norm(E1(:),1) + norm(E2(:),1) + norm(E3(:),1) + norm(E4(:),1) + norm(E5(:),1) + norm(E6(:),1));
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
PS.P1 = P1;PS.P2 = P2;PS.P3 = P3;PS.P4 = P4;PS.P5 = P5;PS.P6 = P6;
ES.E1 = E1;ES.E2 = E2;ES.E3 = E3;ES.E4 = E4;ES.E5 = E5;ES.E6 = E6;
end