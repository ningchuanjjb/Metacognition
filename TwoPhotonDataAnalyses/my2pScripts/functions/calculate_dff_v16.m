function calculate_dff_v16(path,options_dff)
%% Initialization
%clear all
% home;

if(~exist('path','var'))
    %path = 'C:\ASDROOT\STUDY\temp';
    %path = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230511A_Ding_Site05\Result20230513T132357\plane0';
    path = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230528A_Ding_Site02B_sameFOV0527\Result20230529T105723\plane0';    
    %path = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230615A_Ding_Site05D_sameFOV0614\Result20230616T000339\plane0';                    
end
if(~exist('options_dff','var'))
    df_window = 1000;
    if_save = 0;
    if_plot = 1;
    roi_ID_forPlot = 523+1;
    %roi_ID_forPlot = 1170+1;
else
    df_window = options_dff.df_window;
    if_save = options_dff.if_save_dff;
    if_plot = options_dff.if_plot_dff;
    roi_ID_forPlot = options_dff.roi_ID_forPlot+1;
end

if(~exist('if_compute','var'))
    if_compute = 1;
else
    if_compute = 0;%1
end

% if(~exist('df_window','var'))
%     df_window = 1000;
% end
% if(~exist('if_save','var'))
%     if_save = 0;
% end
% if(~exist('if_plot','var'))
%     if_plot = 1;
% end

path_my2pScripts = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts'; %#ok<*UNRCH>
% cd(path_my2pScripts)

% if_compute = 1;%1
% df_window = 1000;%500-->1000
% roi_ID_forPlot = 15+1;
% T = length(temp_F);
T = 380000;%6200-->20000
T0 =1;%1

fileName_Fall = 'Fall.mat';
fileName_iscell = 'iscell.npy';

fullFileName_Fall = [path,'\',fileName_Fall];
fullFileName_iscell = [path,'\',fileName_iscell];


%% Main block for dff
if if_compute == 1
    load(fullFileName_Fall,'F');
    iscell = readNPY(fullFileName_iscell);

    F_raw = F + 1;
    %F = F_raw(iscell(:,1)==1,:);
    %fprintf('Get %d valid F.\n',size(F,1));
    
    dffName = 'dff_jjb';
    dff_v = autoGetFunName_myScripts(dffName, [path_my2pScripts,'\','functions']);
    fun_dff = str2func(dff_v);
    %fprintf('Choose to run %s.\n',dff_v);
    
    t0 = tic;
    %profile on
    [F_dff_raw,F0_raw] = fun_dff(F_raw,df_window,path_my2pScripts);
    %F_dff = F_dff_raw(iscell(:,1)==1,:);
    %F0 = F0_raw(iscell(:,1)==1,:);
    %profile viewer
    %fprintf('Time of dff is %.2f secs.\n',toc(t0));
    
    if if_save == 1
        %save(fullFileName_Fall,'F_dff','-append','-nocompression');
        %save(fullFileName_Fall,'F0','-append','-nocompression');
        save(fullFileName_Fall,'F_dff_raw','-append','-nocompression');
        save(fullFileName_Fall,'F0_raw','-append','-nocompression');        
        %fprintf('Time of save is %.2f secs.\n',toc(t0));
    end
    fprintf('Time of dff is %.2f secs.\n',toc(t0));
end


%% Plot
if if_plot == 1
    %temp_id2 = sum(iscell(1:roi_ID_forPlot,1));
    %temp_F = F(temp_id2,1:T);
    %temp_F0 = F0(temp_id2,1:T);
    %temp_F_dff = F_dff(temp_id2,1:T);
    
    temp_id2 = roi_ID_forPlot;
    temp_F = F_raw(temp_id2,1:T);
    temp_F0 = F0_raw(temp_id2,1:T);
    temp_F_dff = F_dff_raw(temp_id2,1:T);
        
    close all
    fig1 = figure('Name','Fig1','NumberTitle','off');
    set(gcf,'Position',[35+118 35+225 1600 550]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point    
    t = tiledlayout(4,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
    
    nexttile
    plot(temp_F, '-', 'LineWidth', 1, 'Color', [0.3010 0.7450 0.9330]);
    hold on
    plot(temp_F0, '-', 'LineWidth', 2, 'Color', [0.40 0.40 0.40]);
    hold on
    xlim([T0 T]);
    set(gca,'yticklabel',[])
    set(gca, 'FontSize', 11)
    set(gca,'box','off');% 取消右、上边框
    
    nexttile
    plot(temp_F0, '-', 'LineWidth', 2, 'Color', [0.40 0.40 0.40]);
    hold on    
    xlim([T0 T]);
    set(gca,'yticklabel',[])
    set(gca, 'FontSize', 11)
    set(gca,'box','off');% 取消右、上边框
    
    nexttile
    plot(temp_F_dff, '-', 'LineWidth', 1, 'Color', [0.4660 0.6740 0.1880]);
    xlim([T0 T]);
    set(gca,'yticklabel',[])
    set(gca, 'FontSize', 11)
    set(gca,'box','off');% 取消右、上边框    
end

a = 1;
%% END
