clear all
close all

template_0105A = read_file('template_20230105A_Ding_Site09_20230108T184435.tif');
template_0106A_rotate = read_file('template_20230106A_Ding_Site09_sameFOV0105A_justTest_20230109T000802.tif');

template_1 = gpuArray(single(template_0105A));
template_2_raw = gpuArray(single(template_0106A_rotate));
clear template_0*


median_template_1 = median(template_1,'all');
positive_index = template_2_raw>0;
template_2_raw = template_2_raw.*positive_index + median_template_1.*(~positive_index);
median_template_2_raw = median(template_2_raw,'all');
% template_2 = template_2_raw + median_template_1 - median_template_2_raw;

[min_template_1,max_template_1] = bounds(template_1,'all');
[min_template_2_raw,max_template_2_raw] = bounds(template_2_raw,'all');
template_2_n11n = (template_2_raw-min_template_2_raw)./(max_template_2_raw-min_template_2_raw);
template_2_n11n = template_2_n11n.*(max_template_1-min_template_1)+min_template_1;
template_2 = template_2_n11n;

median_template_2 = median(template_2,'all');
[min_template_2,max_template_2] = bounds(template_2,'all');
clear template_2_n11n template_2_raw positive_index min* max*

template_1;
template_2;
% threshold = 300;
% bw_template_1 = template_1.*(template_1>threshold);
% bw_template_2 = template_2.*(template_2>threshold);
Sensitivity = 0.4;
bw_template_1 = imbinarize(gather(template_1./max(template_1,[],'all')),'adaptive','Sensitivity',Sensitivity);
bw_template_2 = imbinarize(gather(template_2./max(template_2,[],'all')),'adaptive','Sensitivity',Sensitivity);

% I_1 = template_1;
% I_2 = template_2;
I_1 = bw_template_1;
% I_1 = template_1;
% I_2 = bw_template_1;



% template_xcorr = xcorr2(template_1,template_2);
% [~,max_index] = max(template_xcorr,[],'all','linear');
% [max_row,max_col] = ind2sub(size(template_xcorr),max_index);
% max_row = max_row - 512;
% max_col = max_col - 512;
% I_4 = template_xcorr;

[~,Greg] = dftregistration_min_max_gpu_jjb_v10(...
    fftn(template_1),fftn(template_2),50,[-32 -32],[32 32],1,1,1);
template_2_shift = real(ifftn(Greg));
positive_index = template_2_shift > 0;
template_2_shift = template_2_shift.*positive_index + median_template_2.*(~positive_index);
overMax_index = template_2_shift > max(template_2,[],'all');
underMin_index = template_2_shift < min(template_2,[],'all');

template_2_shift = template_2_shift.*(~overMax_index) + max(template_2,[],'all').*overMax_index;
template_2_shift = template_2_shift.*(~underMin_index) + min(template_2,[],'all').*underMin_index;

bw_template_2_shift = imbinarize(gather(template_2_shift./max(template_2_shift,[],'all')),'adaptive','Sensitivity',Sensitivity);
I_4 = template_2_shift;
I_2 = bw_template_2_shift;    

% template_corr = corr(template_1,template_2_shift);
% template_corr = abs(corr(template_1,template_2_shift));
% template_corr = corr(bw_template_1,bw_template_2_shift);
% template_corr = abs(corr(bw_template_1,bw_template_2_shift));
% median_template_corr_raw  = median(template_corr,'all','omitnan');
% nan_index = isnan(template_corr);
% template_corr(nan_index) = median_template_corr_raw;
% negative_index = template_corr<0;
% template_corr = template_corr.*(~negative_index) + median_template_corr_raw.*negative_index;
% 
% I_3 = template_corr;

% template_intersect = bw_template_1 & bw_template_2;
template_intersect = bw_template_1 & bw_template_2_shift;
% I_3 = template_intersect;

options = NoRMCorreSetParms(...
    'd1',512,...
    'd2',512,...
    'grid_size',[32,32],...
    'overlap_pre',[32,32],...
    'overlap_post',[32,32],...
    'us_fac',50,...%50
    'max_dev',[8 8],...
    'max_shift',[32 32],...  %[64 64]
    'phase_flag',true,... %true
    'upd_template',false,... %true  
    'iter',1,... %2
    'boundary','zero',...
    'output_type','mat');
[M_final,shifts] = normcorre(bw_template_2_shift,options,bw_template_1);

overMax_index = M_final > 1;
underMin_index = M_final < 0;
M_final = M_final.*(~overMax_index) + 1.*overMax_index;
M_final = M_final.*(~underMin_index) + 0.*underMin_index;
M_final_raw = M_final;
bw_index = M_final > 0.5;
% M_final = M_final.*bw_index + 0.*(~bw_index);
M_final = bw_index;
I_3 = M_final;


% I_1 = template_1;
% I_2 = template_2;

fig1 = figure('Name','Fig1, I_1','NumberTitle','off');
set(gcf,'Position',[0 0 1920 1080]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(I_1.^(1/6));
colormap(gray);
axis off equal

fig2 = figure('Name','Fig2, I_2','NumberTitle','off');
set(gcf,'Position',[0 0 1920 1080]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(I_2.^(1/6));
colormap(gray);
axis off equal

fig3 = figure('Name','Fig3, I_3','NumberTitle','off');
set(gcf,'Position',[0 0 1920 1080]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(I_3.^(1/6));
colormap(gray);
axis off equal

% fig4 = figure('Name','Fig4, I_4','NumberTitle','off');
% set(gcf,'Position',[0 0 1920 1080]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
% set (gca,'position',[0.01,0.01,0.98,0.98] )
% % imagesc(I_4);
% imagesc(I_4.^(1/6));
% colormap(gray);
% axis off equal

% fig5 = figure('Name','Fig4, I_5','NumberTitle','off');
% set(gcf,'Position',[0 0 1920 1080]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
% set (gca,'position',[0.01,0.01,0.98,0.98] )
% % imagesc(I_4);
% imagesc(I_5.^(1/6));
% colormap(gray);
% axis off equal