function X = shrink_l1(B,lambda)

% The shrinkage operator of the l1 norm
% 
% min_X lambda*||X||_1+0.5*||X-B||_F^2
%
% version 1.0 - 19/07/2021
%
% Written by wahaha (dingjie@njtech.edu.cn)
% 

X = sign(B) .* max(abs(B)-lambda,0);