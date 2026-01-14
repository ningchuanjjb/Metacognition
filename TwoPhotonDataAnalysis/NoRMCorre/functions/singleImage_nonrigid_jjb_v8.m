function [D,temp_reg] = singleImage_nonrigid_jjb_v8(template_in, template)
%template_in, template_in_raw; % This is current template.
%template, template_n11n; % This is existing template.

a = 1;
if_profile = 0;
if_showTime = 0;
if if_profile == 1
    profile on
end
t0 = tic;

positive_index = (template_in > 0);
template_in = template_in.*positive_index;
template_in_raw = template_in;

[min_template_in_raw,max_template_in_raw] = bounds(gather(template_in_raw),'all');
[min_template,max_template] = bounds(gather(template),'all'); %#ok<*ASGLU>
median_template = median(template,'all');

template_n11n = rescale(template,min_template_in_raw,max_template_in_raw);
a = 1;

options_nonrigid = NoRMCorreSetParms(...
    'd1',512,...
    'd2',512,...
    'grid_size',[32,32],...%48
    'overlap_pre',[32,32],...%48
    'overlap_post',[32,32],...%48
    'us_fac',50,...%50
    'max_dev',[12 12],...%8
    'max_shift',[48 48],...  %[48 48]
    'phase_flag',false,... %true-->false
    'upd_template',false,... %true
    'iter',3,... %2
    'boundary','zero',...
    'output_type','mat',...
    'correct_bidir',true);

options_nonrigid.shifts_method = 'cubic';
options_nonrigid.iter = 1;


% %% To fit nBins_histeq, which is for histogram equlization
% optimal_NumTiles_histeq = 8;
% 
% gather_template_n11n = gather(template_n11n);
% gather_template_in = gather(template_in);
% initial_template_n11n = rescale(gather_template_n11n,0,1);
% initial_template_in = rescale(gather_template_in,0,1);
% 
% % temp_start = 8;
% % temp_end = 10;%13-->10
% % temp_step = 1;%0.05-->0.5-->0.7-->1
% % tempN = floor((temp_end-temp_start)/temp_step)+1;
% % temp_nBins_histeq_power = temp_start:temp_step:temp_end;
% % rho = zeros(1,tempN);
% % 
% % temp_options_nonrigid = options_nonrigid;
% % t0_temp = tic;
% % % parfor tempi=1:tempN
% % for tempi=1:tempN
% %     temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',floor(2^(temp_nBins_histeq_power(tempi))),'Range','original','Distribution','rayleigh','Alpha',0.4);
% %     temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',floor(2^(temp_nBins_histeq_power(tempi))),'Range','original','Distribution','rayleigh','Alpha',0.4);
% %     temp_output = normcorre_jjb_v2(temp_template_in,temp_options_nonrigid,temp_template_n11n);
% %     rho(tempi) = corr2(temp_output,temp_template_n11n);
% % end
% % [~,I]=max(rho);
% % optimal_nBins_histeq_power = temp_nBins_histeq_power(I);
% % optimal_nBins_histeq = floor(2^optimal_nBins_histeq_power);
% % if if_showTime == 1
% %     fprintf(' optimal_nBins_histeq=%d, time=%.1f secs.',optimal_nBins_histeq,toc(t0_temp));
% % end
% 
% optimal_nBins_histeq = 256;
% 
% 
% temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',0.4);
% temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',0.4);
% gather_template_n11n = rescale(temp_template_n11n,min_template_in_raw,max_template_in_raw); %#ok<*NASGU>
% gather_template_in = rescale(temp_template_in,min_template_in_raw,max_template_in_raw);
% 
% 
% %% To fit alpha_histeq, which is for histogram equlization
% gather_template_n11n = gather(template_n11n);
% gather_template_in = gather(template_in);
% initial_template_n11n = rescale(gather_template_n11n,0,1);
% initial_template_in = rescale(gather_template_in,0,1);
% 
% % temp_start = 0.1;%0.1
% % temp_end = 1;%10-->5-->3-->1
% % temp_step = 0.4;%0.1-->0.025-->0.01-->0.06-->0.09
% % tempN = floor((temp_end-temp_start)/temp_step)+1;
% % temp_alpha_histeq = temp_start:temp_step:temp_end;
% % rho = zeros(1,tempN);
% % 
% % temp_options_nonrigid = options_nonrigid;
% % t0_temp = tic;
% % % parfor tempi=1:tempN
% % for tempi=1:tempN
% %     temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',temp_alpha_histeq(tempi));
% %     temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',temp_alpha_histeq(tempi));
% %     temp_output = normcorre_jjb_v2(temp_template_in,temp_options_nonrigid,temp_template_n11n);
% %     rho(tempi) = corr2(temp_output,temp_template_n11n);
% % end
% % [~,I]=max(rho);
% % optimal_alpha_histeq = temp_alpha_histeq(I);
% % if if_showTime == 1
% %     fprintf(' optimal_alpha_histeq=%.3f, time=%.1f secs.',optimal_alpha_histeq,toc(t0_temp));
% % end
% 
% optimal_alpha_histeq = 0.177777;
% 
% temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);
% temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);
% gather_template_n11n = rescale(temp_template_n11n,min_template_in_raw,max_template_in_raw);
% gather_template_in = rescale(temp_template_in,min_template_in_raw,max_template_in_raw);
% 
% 
% %% To fit NumTiles_histeq, which is for histogram equlization
% gather_template_n11n = gather(template_n11n);
% gather_template_in = gather(template_in);
% initial_template_n11n = rescale(gather_template_n11n,0,1);
% initial_template_in = rescale(gather_template_in,0,1);
% 
% % temp_start = 8;%2
% % temp_end = 64;%128
% % temp_step = 18;%1-->4-->7
% % tempN = floor((temp_end-temp_start)/temp_step)+1;
% % temp_NumTiles_histeq = temp_start:temp_step:temp_end;
% % rho = zeros(1,tempN);
% % 
% % temp_options_nonrigid = options_nonrigid;
% % t0_temp = tic;
% % % parfor tempi=1:tempN
% % for tempi=1:tempN
% %     temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*temp_NumTiles_histeq(tempi),'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);
% %     temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*temp_NumTiles_histeq(tempi),'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);
% %     temp_output = normcorre_jjb_v2(temp_template_in,temp_options_nonrigid,temp_template_n11n);
% %     rho(tempi) = corr2(temp_output,temp_template_n11n);
% % end
% % [~,I]=max(rho);
% % optimal_NumTiles_histeq = temp_NumTiles_histeq(I);
% % if if_showTime == 1
% %     fprintf(' optimal_NumTiles_histeq=%d, time=%.1f secs.',optimal_NumTiles_histeq,toc(t0_temp));
% % end
% 
% optimal_NumTiles_histeq = 8;
% 
% temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);
% temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);
% gather_template_n11n = rescale(temp_template_n11n,min_template_in_raw,max_template_in_raw); %#ok<*NASGU>
% gather_template_in = rescale(temp_template_in,min_template_in_raw,max_template_in_raw);


%% adapthisteq
gather_template_n11n = gather(template_n11n);
gather_template_in = gather(template_in);
initial_template_n11n = rescale(gather_template_n11n,0,1);
initial_template_in = rescale(gather_template_in,0,1);

optimal_NumTiles_histeq = 32;%8
optimal_nBins_histeq = 4096*2;%256
optimal_alpha_histeq = 0.11777;%0.4,0.11777
ClipLimit = 0.1;%0.1

temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,...
    'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq,'ClipLimit',ClipLimit);
gather_template_n11n = rescale(temp_template_n11n,min_template_in_raw,max_template_in_raw); %#ok<*NASGU>

temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,...
    'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq,'ClipLimit',ClipLimit);
gather_template_in = rescale(temp_template_in,min_template_in_raw,max_template_in_raw);

max_gather_template_n11n = max(gather_template_n11n,[],'all');
max_gather_template_in = max(gather_template_in,[],'all');

%% To fit optimal grid size
temp_start = 24;%32-->24-->16
temp_end = 64;%100
temp_step = 10;%1-->3-->5-->10
tempN = floor((temp_end-temp_start)/temp_step)+1;
temp_grid_size = temp_start:temp_step:temp_end;
rho = zeros(1,tempN);

t0_temp = tic;
% parfor tempi=1:tempN
for tempi=1:tempN
    temp_options_nonrigid = options_nonrigid;
    temp_options_nonrigid.grid_size = [temp_grid_size(tempi),temp_grid_size(tempi),1];
    temp_options_nonrigid.overlap_pre = [temp_grid_size(tempi),temp_grid_size(tempi),1];
    temp_options_nonrigid.overlap_post = [temp_grid_size(tempi),temp_grid_size(tempi),1];
    temp_output = normcorre_jjb_v2(gather_template_in,temp_options_nonrigid,gather_template_n11n);
    rho(tempi) = corr2(temp_output,gather_template_n11n);
end
[~,I]=max(rho);
optimal_grid_size = temp_grid_size(I);
optimal_overlap = optimal_grid_size;
options_nonrigid.grid_size = [optimal_grid_size,optimal_grid_size,1];
options_nonrigid.overlap_pre = [optimal_overlap,optimal_overlap,1];
options_nonrigid.overlap_post = [optimal_overlap,optimal_overlap,1];
if if_showTime == 1
    fprintf(' optimal_grid_size=%d, time=%.1f secs.',optimal_grid_size,toc(t0_temp));
end


%% To fit optimal max_dev
temp_start = 2;
temp_end = 15;
temp_step = 6;%1-->2-->4
tempN = floor((temp_end-temp_start)/temp_step)+1;
temp_max_dev = temp_start:temp_step:temp_end;
rho = zeros(1,tempN);

t0_temp = tic;
% parfor tempi=1:tempN
for tempi=1:tempN
    temp_options_nonrigid = options_nonrigid;
    temp_options_nonrigid.max_dev = [temp_max_dev(tempi),temp_max_dev(tempi),1];
    temp_output = normcorre_jjb_v2(gather_template_in,temp_options_nonrigid,gather_template_n11n);
    rho(tempi) = corr2(temp_output,gather_template_n11n);
end
[~,I]=max(rho);
optimal_max_dev = temp_max_dev(I);
options_nonrigid.max_dev = [optimal_max_dev,optimal_max_dev,1];
if if_showTime == 1
    fprintf(' optimal_max_dev=%d, time=%.1f secs.\n',optimal_max_dev,toc(t0_temp));
end


%% To fit optimal contrast_ratio by fminsearch
% range_contrast_ratio = [0.3,3.01];%[0.99,4]-->[0.99,3.01]-->[0.1,3.01]-->[0.3,3.01]
% params_range = range_contrast_ratio;
% 
% ModelName = 'nonRigidParamFit_v2';
% Model = str2func(ModelName);
% func = @(params) Model(params,gather_template_in,options_nonrigid,gather_template_n11n,params_range);
% 
% % Use fminsearch
% optimset_options = optimset('fminsearch');
% optimset_options.TolFun = 1e-7;%1e-6
% optimset_options.TolX = 1e-3;%1e-1
% optimset_options.MaxFunEvals = 2;%40-->3
% optimset_options.MaxIter = 2;%40-->25-->10-->5-->3
% optimset_options.Display = 'off'; % off, iter, final
% optimset_options.FunValCheck = 'on';
% optimset_options.PlotFcns = [];%@optimplotfval
% 
% 
% initial_contrast_ratio = 0.5:0.3:1.8;%1:0.1428:2 --> 1:0.28570:3 --> 0.5:0.1:1.8 --> 0.5:0.2:1.8
% contrast_ratio_raw = ...
%     -log(1./((initial_contrast_ratio-range_contrast_ratio(1))./(range_contrast_ratio(2)-range_contrast_ratio(1)))-1);
% initialParams = contrast_ratio_raw;
% tempN = length(initialParams);
% rho = zeros(1,tempN);
% temp_modelParams = zeros(1,tempN);
% 
% t0_contrast = tic;
% % parfor tempi=1:tempN
% for tempi=1:tempN
%     temp_optimset_options = optimset_options;
%     if tempi==1
%         temp_optimset_options.Display = 'off'; % off, iter, final
%     end
%     temp_modelParams(tempi) = fminsearch(func, initialParams(tempi), temp_optimset_options);
%     [~,~,rho(tempi)] = Model(temp_modelParams(tempi),gather_template_in,options_nonrigid,gather_template_n11n,params_range); %#ok<*PFBNS>
%     a = 1;
% end
% a = 1;
% [rho2,I]=max(rho);
% modelParams = temp_modelParams(I);
% 
% 
% contrast_ratio_raw = modelParams(1);
% contrast_ratio = (1/(1+exp(-contrast_ratio_raw)))*(range_contrast_ratio(2)-range_contrast_ratio(1))+range_contrast_ratio(1);% constrain in range_contrast_ratio
% fprintf('fimsearch time=%.1f seconds, rho=%.5f, contrast_ratio=%.3f(from%.3f)\n',...
%     toc(t0_contrast),rho2,contrast_ratio,initial_contrast_ratio(I));

contrast_ratio = 1;%0.8

gather_template_in = gather_template_in.^(1/contrast_ratio);
gather_template_n11n = gather_template_n11n.^(1/contrast_ratio);
[min_gather_template_in,max_gather_template_in] = bounds(gather_template_in,'all');

% To do coarse correction firstly
options_nonrigid.iter = 3;%3
% options_nonrigid.shifts_method = 'FFT';% slow but good
options_nonrigid.shifts_method = 'cubic';% fast

t0_temp = tic;
M1_final = normcorre_jjb_v2(gather_template_in,options_nonrigid,gather_template_n11n);
if if_showTime == 1
    fprintf('time (normcorre_jjb_v2)=%.1f secs.\n',toc(t0_temp));
end

overMax_index = M1_final > max_gather_template_in;
underMin_index = M1_final < min_gather_template_in;
M1_final = M1_final.*(~overMax_index) + max_gather_template_in.*overMax_index;
M1_final = M1_final.*(~underMin_index) + min_gather_template_in.*underMin_index;

if_nonRigidRefine = 1;
if if_nonRigidRefine == 1
    % To do refined correction secondly
    t0_temp = tic;
    %options_nonrigid.shifts_method = 'cubic';
    options_nonrigid.overlap_pre = [16,16,1];%16
    options_nonrigid.overlap_post = [16,16,1];
    
    temp_start = 64;
    temp_end = 32;
    temp_step = -4;%-1-->-2-->-4
    rho_0 = corr2(M1_final,gather_template_n11n);
    for tempi=temp_start:temp_step:temp_end %64 to 32
        options_nonrigid.grid_size = [tempi,tempi,1];
        temp_output = normcorre_jjb_v2(M1_final,options_nonrigid,gather_template_n11n);
        temp_rho = corr2(temp_output,gather_template_n11n);
        if temp_rho > rho_0
            M1_final = temp_output;
            rho_0 = temp_rho;
            
            overMax_index = M1_final > max_gather_template_in;
            underMin_index = M1_final < min_gather_template_in;
            M1_final = M1_final.*(~overMax_index) + max_gather_template_in.*overMax_index;
            M1_final = M1_final.*(~underMin_index) + min_gather_template_in.*underMin_index;
        end
    end
    if if_showTime == 1
        fprintf('Nonrigid refine time=%.2f seconds.\n',toc(t0_temp));
    end
end


a = 1;

% t0_temp = tic;

% when params is [9,1000], it leads huge fake alignment in 0509 FOV, so here adjust to [5,100]
PyramidLevels = 5;%9
iterNum = 100;%1000
temp_iterRange = ceil(iterNum:-iterNum/PyramidLevels:1);

[D,temp_reg] = imregdemons(gather_template_in,gather(M1_final),...
    temp_iterRange,...
    'AccumulatedFieldSmoothing',1,...
    'PyramidLevels',PyramidLevels,...
    'DisplayWaitbar',false);

% fprintf('time (imregdemons)=%.1f secs.\n',toc(t0_temp));

fprintf('Single image non-rigid done, time=%.1f secs.\n',toc(t0));

clear M1_final 
clear max_template_in_raw min_template_in_raw template_n11n


if if_profile == 1
    profile viewer
end

a = 1;

%% End
