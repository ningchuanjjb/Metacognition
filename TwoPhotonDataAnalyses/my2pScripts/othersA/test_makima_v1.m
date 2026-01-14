a = 1;
dx = diff(Xpcs,1,2);
tempA = kron(Xpcs(:,1:end-1),ones(1,shift));
tempB = kron(dx,1/shift*(0:shift-1));
Xf = tempA + tempB;

temp_length = (length(Xpcs)-1)*shift; %#ok<*NASGU>

% x = [0 1 2.5 3.6 5 7 8.1 10];
% y = cos(x);
% xq = 0:.25:10;
% yq = makima(x,y,xq);


x = 1:(length(Xpcs)-1);
y = Xpcs(:,1:end-1);
xq = linspace(min(x),max(x),temp_length);
yq_raw = makima(x,y,xq);
yq = [repmat(yq_raw(1),1,shift/2) yq_raw(1:end-shift/2)];

close all
%% Plot
fig1 = figure('Name','Fig1','NumberTitle','off');
set(gcf,'Position',[35+350 35+225 1000 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

nexttile
plot(Xf, '-', 'LineWidth', 1, 'Color', [0.3010 0.7450 0.9330]);
hold on
plot(yq, '-', 'LineWidth', 1, 'Color', [0.40 0.40 0.40]);
hold on
% xlim([1 T]);
% ylim([0 temp_max]);
% set(gca,'yticklabel',[])
set(gca, 'FontSize', 11)
set(gca,'box','off');% 取消右、上边框



