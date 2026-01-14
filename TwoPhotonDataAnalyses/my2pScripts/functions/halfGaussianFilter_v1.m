function [B,w] = halfGaussianFilter_v1(A,windowSize,alpha)

if(~exist('alpha','var'))
    alpha = 2.5;%2.5-->4-->6
end

L = 1+2*(windowSize-1);
w_raw = gausswin(L,alpha);
w = w_raw(windowSize:end);
w = w ./ sum(w);
w = w';
% B = conv(A,w);
B = conv2(A,w);
% B = conv2(A,w,'same');
a = 1;
temp1 = 3;%4-->3
% B = B(temp1+(1:length(A)));
B = B(:,temp1+(1:length(A)));
