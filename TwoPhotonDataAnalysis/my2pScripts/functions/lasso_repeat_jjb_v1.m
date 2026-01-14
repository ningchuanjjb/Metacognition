function [B,r2] = lasso_repeat_jjb_v1(x,y,repeatNum,fun_lasso_jjb)

% repeatNum = 10;
numFrames = size(x,2);
temp_B = zeros(numFrames,repeatNum);
temp_r2 = zeros(1,repeatNum);
temp_B0 = zeros(1,repeatNum);
for tempj=1:repeatNum
    [temp_B(:,tempj),temp_r2(tempj),temp_B0(tempj)] = fun_lasso_jjb(x,y); %#ok<*PFBNS>
end
% B = mean(temp_B,2);
% B0 = mean(temp_B0);

% B = median(temp_B,2);
% B0 = median(temp_B0);

[~,tempI] = max(temp_r2);
B = temp_B(:,tempI);
B0 = temp_B0(tempI);


y_hat = double(x) * B + B0;
R = y-y_hat;
r2 = 1 - sum(R.^2)/sum(((y-mean(y)).^2));

a = 1;