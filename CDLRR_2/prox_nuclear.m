function [X,nuclearnorm] = prox_nuclear(B,lambda)

% The proximal operator of the nuclear norm of a matrix
% 
% min_X lambda*||X||_*+0.5*||X-B||_F^2
%
% version 1.0 - 18/06/2016
%
% Written by Canyi Lu (canyilu@gmail.com)
% 

[U,S,V] = svd(B, 'econ');   %对矩阵B进行精简分解
S = diag(S);    %去除S中的所有0，并将非0值置为1列
svp = length(find(S>lambda));   %求出S中的值大于lambda个数
if svp>=1
    S = S(1:svp) - lambda;  %S中前svp个数均减去lambda,此步将S中小于lambda的值舍去，实现降秩
    X = U(:,1:svp)*diag(S)*V(:,1:svp)'; %X的秩小于等于B的秩 S的维度发生变化，变为[svp,svp]
    nuclearnorm = sum(S);
else
    X = zeros(size(B));
    nuclearnorm = 0;
end
