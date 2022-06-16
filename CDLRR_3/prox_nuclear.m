function [X,nuclearnorm] = prox_nuclear(B,lambda)

% The proximal operator of the nuclear norm of a matrix
% 
% min_X lambda*||X||_*+0.5*||X-B||_F^2
%
% version 1.0 - 18/06/2016
%
% Written by Canyi Lu (canyilu@gmail.com)
% 

[U,S,V] = svd(B, 'econ');   %�Ծ���B���о���ֽ�
S = diag(S);    %ȥ��S�е�����0��������0ֵ��Ϊ1��
svp = length(find(S>lambda));   %���S�е�ֵ����lambda����
if svp>=1
    S = S(1:svp) - lambda;  %S��ǰsvp��������ȥlambda,�˲���S��С��lambda��ֵ��ȥ��ʵ�ֽ���
    X = U(:,1:svp)*diag(S)*V(:,1:svp)'; %X����С�ڵ���B���� S��ά�ȷ����仯����Ϊ[svp,svp]
    nuclearnorm = sum(S);
else
    X = zeros(size(B));
    nuclearnorm = 0;
end
