function [B,r2,B0] = lasso_jjb_v6(x,y,weight)


% [B_raw,FitInfo] = lasso(double(x),y,'CV',3,'Standardize',false);%3-->5
%[B_raw,FitInfo] = lasso(double(x),y,'Standardize',false);
% [B_raw,FitInfo] = lasso(double(x),y,'Standardize',false,'NumLambda',50,'RelTol',1e-3,'Alpha',0.5);%NumLambda=100, RelTol=1e-4, Alpha=0.5
% [B_raw,FitInfo] = lasso(double(x),y,'Standardize',false,'NumLambda',50,'RelTol',1e-3);%NumLambda=100, RelTol=1e-4, Alpha=0.5
% [B_raw,FitInfo] = lasso(double(x),y,'CV',3,'Standardize',false,'MCReps',3);%3-->5


% [B_raw,FitInfo] = lasso(double(x),y,'Standardize',false,'Weights',weight,'NumLambda',50,'RelTol',1e-3);%NumLambda=100, RelTol=1e-4, Alpha=0.5
% [B_raw,FitInfo] = lasso(double(x),y,'CV',3,'Standardize',false,'MCReps',3,'Weights',weight);%3-->5
% [B_raw,FitInfo] = lasso(double(x),y,'CV',3,'Standardize',false,'MCReps',7,'Weights',weight);
% [B_raw,FitInfo] = lasso(double(x),y,'CV',3,'Standardize',false,'MCReps',8,'Weights',weight);% actual time  is *(3/8)*0.94
% [B_raw,FitInfo] = lasso(double(x),y,'CV',3,'Standardize',false,'MCReps',100,'Weights',weight);%3-->5
% [B_raw,FitInfo] = lasso(double(x),y,'CV',3,'Standardize',false,'MCReps',1,'Weights',weight);%3-->5

CV_num = 3;
if size(x,1) >= CV_num
    [B_raw,FitInfo] = lasso(double(x),y,'CV',CV_num,'Standardize',false,'MCReps',1,'Weights',weight);%3-->5
else
    [B_raw,FitInfo] = lasso(double(x),y,'Standardize',false,'MCReps',1,'Weights',weight);%3-->5
end

[MinMSE,IndexMinMSE] = min(FitInfo.MSE);
B = B_raw(:,IndexMinMSE);
B0 = FitInfo.Intercept(IndexMinMSE);
r2 = 1-length(y)*MinMSE/sum((y-mean(y)).^2);

