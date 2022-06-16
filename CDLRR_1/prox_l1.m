function X = prox_l1(B,lambda)

% The proximal operator of the l1 norm
% 
% min_x lambda*||X||_1+0.5*||X-B||_F^2
%
% version 1.0 - 19/07/2021
%
% Written by wahaha (dingjie@njtech.edu.cn)
% 

X = max(0,B-lambda)+min(0,B+lambda);
