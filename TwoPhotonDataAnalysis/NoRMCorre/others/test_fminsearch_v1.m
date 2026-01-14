
range_contrast_ratio = [1,6];
% range_grid_size = [24,64];
% range_max_dev = [2,15];
% params_range = [range_contrast_ratio,range_grid_size,range_max_dev];
params_range = range_contrast_ratio;

initial_contrast_ratio = 2;
% initial_grid_size = 48;
% initial_max_dev = 8;
contrast_ratio_raw = ...
    -log(1/((initial_contrast_ratio-range_contrast_ratio(1))/(range_contrast_ratio(2)-range_contrast_ratio(1)))-1);
% grid_size_raw = ...
%     -log(1/((initial_grid_size-range_grid_size(1))/(range_grid_size(2)-range_grid_size(1)))-1);
% max_dev_raw = ...
%     -log(1/((initial_max_dev-range_max_dev(1))/(range_max_dev(2)-range_max_dev(1)))-1);
% initialParams = [contrast_ratio_raw,grid_size_raw,max_dev_raw];
initialParams = contrast_ratio_raw;

ModelName = 'nonRigidParamFit_v2';
Model = str2func(ModelName);
func = @(params) Model(params,gather_template_in,options_nonrigid,gather_template_n11n,params_range);

% Use fminsearch
options = optimset('fminsearch');
options.TolFun = 1e-5;%1e-7
options.TolX = 1e-1;%1e-1
options.MaxFunEvals = 40;
options.Display = 'iter'; % off, iter, final
options.FunValCheck = 'on';
options.PlotFcns = [];%@optimplotfval


t0 = tic;
modelParams = fminsearch(func, initialParams, options);
[loss,M_output,rho] = Model(modelParams,gather_template_in,options_nonrigid,gather_template_n11n,params_range);
fprintf('fimsearch time=%.1f seconds, rho=%.5f\n',toc(t0),rho);

contrast_ratio_raw = modelParams(1);
% grid_size_raw = modelParams(2);
% max_dev_raw = modelParams(3);
contrast_ratio = (1/(1+exp(-contrast_ratio_raw)))*(range_contrast_ratio(2)-range_contrast_ratio(1))+range_contrast_ratio(1);% constrain in (1, 6)
% grid_size = round((1/(1+exp(-grid_size_raw)))*(range_grid_size(2)-range_grid_size(2))+range_grid_size(1));% constrain in (24, 64), int
% max_dev = round((1/(1+exp(-max_dev_raw)))*(range_max_dev(2)-range_max_dev(1))+range_max_dev(1));% constrain in (2, 15), int


a = 1;

