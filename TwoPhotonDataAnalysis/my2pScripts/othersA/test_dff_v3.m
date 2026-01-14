%% Initialization
%clear all
home;

path_my2pScripts = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts'; %#ok<*UNRCH>
% cd(path_my2pScripts)

if_compute = 1;
df_window = 1000;%500-->1000
temp_id = 15+1;
% T = length(temp_F);
T = 16200;%6200

%% Main block for dff
if if_compute == 1      
    load('C:\ASDROOT\STUDY\temp\F.mat')
    load iscell
    F_raw = F;
    F = F_raw(iscell(:,1)==1,:);
    
    dffName = 'dff_jjb';
    dff_v = autoGetFunName_myScripts(dffName, [path_my2pScripts,'\','functions']);
    fun_dff = str2func(dff_v);
    fprintf('Choose to run %s.\n',dff_v);
        
    t0 = tic;
    %profile on
    [F_dff,F0] = fun_dff(F,df_window,path_my2pScripts);
    %profile viewer
    fprintf('Time of dff = %.2f secs.\n',toc(t0))
end

temp_id2 = sum(iscell(1:temp_id,1));
temp_F = F(temp_id2,1:T);
temp_F0 = F0(temp_id2,1:T);
temp_F_dff = F_dff(temp_id2,1:T);

%% Plot
close all
fig1 = figure('Name','Fig1','NumberTitle','off');
set(gcf,'Position',[35+118 35+225 1600 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

nexttile
plot(temp_F, '-', 'LineWidth', 1, 'Color', [0.3010 0.7450 0.9330]);
hold on
plot(temp_F0, '-', 'LineWidth', 1, 'Color', [0.40 0.40 0.40]);
hold on
xlim([1 T]);
set(gca,'yticklabel',[])
set(gca, 'FontSize', 11)
set(gca,'box','off');% 取消右、上边框

nexttile
plot(temp_F0, '-', 'LineWidth', 1, 'Color', [0.40 0.40 0.40]);
xlim([1 T]);
set(gca,'yticklabel',[])
set(gca, 'FontSize', 11)
set(gca,'box','off');% 取消右、上边框

nexttile
plot(temp_F_dff, '-', 'LineWidth', 1, 'Color', [0.6350 0.0780 0.1840]);
xlim([1 T]);
set(gca,'yticklabel',[])
set(gca, 'FontSize', 11)
set(gca,'box','off');% 取消右、上边框

%% END
