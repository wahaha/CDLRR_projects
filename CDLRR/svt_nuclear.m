function [X,nuclearnorm] = svt_nuclear(B,lambda)

% The Sigular Value Thresholding of the nuclear norm of a matrix
% 
% min_X lambda*||X||_*+0.5*||X-B||_F^2
%
% version 1.0 - 19/07/2021
%
% Written by wahaha (dingjie@njtech.edu.cn)
% 

[U,S,V] = svd(B,'econ');
S = diag(S);
S = sign(S) .* max(abs(S)-lambda, 0);
X = U*diag(S)*V';   
nuclearnorm = sum(S);

end