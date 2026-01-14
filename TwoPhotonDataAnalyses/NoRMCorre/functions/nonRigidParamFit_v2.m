function [loss,M_output,rho] = nonRigidParamFit_v2(params,gather_template_in,options_nonrigid,gather_template_n11n,params_range)
%initialParams = [contrast_ratio,grid_size,max_dev];
% params_range = [range_contrast_ratio,range_grid_size,range_max_dev];

temp_options_nonrigid = options_nonrigid;

contrast_ratio_raw = params(1);
% grid_size_raw = params(2);
% max_dev_raw = params(3);

range_contrast_ratio = params_range(1:2);
% range_grid_size = params_range(3:4);
% range_max_dev = params_range(5:6);



% a = 1/(1+exp(-a_raw)); % Use sigmoid to constrain 'a' in (0, 1)
contrast_ratio = (1/(1+exp(-contrast_ratio_raw)))*(range_contrast_ratio(2)-range_contrast_ratio(1))+range_contrast_ratio(1);% constrain in (1, 6)
% grid_size = round((1/(1+exp(-grid_size_raw)))*(range_grid_size(2)-range_grid_size(2))+range_grid_size(1));% constrain in (24, 64), int
% max_dev = round((1/(1+exp(-max_dev_raw)))*(range_max_dev(2)-range_max_dev(1))+range_max_dev(1));% constrain in (2, 15), int

% To inhibit high brightness so that low brightness region could be more accurate
gather_template_in = gather_template_in.^(1/contrast_ratio);
gather_template_n11n = gather_template_n11n.^(1/contrast_ratio);


% temp_options_nonrigid.grid_size = [grid_size,grid_size,1];
% temp_options_nonrigid.overlap_pre = [grid_size,grid_size,1];
% temp_options_nonrigid.overlap_post = [grid_size,grid_size,1];
% temp_options_nonrigid.max_dev = [max_dev,max_dev,1];

M_output = normcorre_jjb_v2(gather_template_in,temp_options_nonrigid,gather_template_n11n);

rho = corr2(M_output,gather_template_n11n);
loss = (1 - rho)^2;



