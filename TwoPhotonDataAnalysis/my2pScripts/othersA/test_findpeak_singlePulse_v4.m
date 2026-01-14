close all

%% part 1

% temp_diff = diff(temp1D);
% temp_diff = [0 temp_diff];
% tempBoolIndex_nagative = temp_diff < 0;
% tempBoolIndex_positive = temp_diff >= 0;
% temp_diff2 = diff(tempBoolIndex_positive);
% temp_valleyBoolIndex = temp_diff2 == 1;
% locs = find(temp_valleyBoolIndex==true);
% valleys = temp_max - temp1D(temp_valleyBoolIndex);
% 
% temp1D_raw = tempID_F_smooth(tempj+(0:threshold_spks_start2));
% if floor(length(temp1D_raw)/2)*2 ~= length(temp1D_raw)
%     temp1D_raw(end) = [];
% end

temp_min = min(temp1D);

fig1 = figure('Name','Fig1','NumberTitle','off');
set(gcf,'Position',[35+0 35+0 900 700]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

nexttile
plot(temp1D,'r');
hold on

% [temp1D_envelope,~] = envelope(temp1D,10,'analytic');%10
% [temp1D_envelope,~] = envelope(temp1D,1,'peak');%10
% temp1D_envelope(1:10) = 0;
% temp1D_envelope(end-9:end) = 0;
% temp1D_envelope = abs(hilbert(temp1D));

[temp_peaks,temp_locs] = findpeaks(temp1D); % Time-consuming!!!
temp_peaks = [temp1D(1) temp_peaks temp1D(end)]; %#ok<*AGROW>
temp_locs = [1 temp_locs size(temp1D,2)];
temp1D_envelope = interp1(temp_locs,temp_peaks,temp_locs(1):temp_locs(end));

plot(temp1D_envelope,'g');
hold on

for temptempi=1:length(locs)
   plot(locs(temptempi)*[1 1],[temp_min, temp_max],'k'); 
   hold on
end

plot(temp_spkEndIndex*[1 1],[temp_min, temp_max],'k','lineWidth',2);
hold on

nexttile
plot(temp1D_raw,'r');
hold on



%% part 2

% tempCalciumPulse = temp1D(1:temp_spkEndIndex);
tempCalciumPulse = temp1D(temp_peakIndex:temp_spkEndIndex);
x = 1:length(tempCalciumPulse);
y = tempCalciumPulse;
% % temp_mdl = fitglm(x,y,'Link','log');
% temp_mdl = fitglm(x,y,'Distribution','poisson');
% % temp_mdl = fitglm(x,y,'Distribution','gamma');
% r2 = temp_mdl.Rsquared.Adjusted;
% 
% 
% % lambdahat = poissfit(tempCalciumPulse);
% 
% x_fit = 1:0.01:length(tempCalciumPulse);
% y_fit = predict(temp_mdl,x_fit')';


% [xData, yData] = prepareCurveData(x,double(y));
% ft = fittype( 'exp1' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = [1.55650486517981 -0.0909553835273515];
% [fitresult, gof] = fit(xData,yData,ft,opts);

% cftool(x,y);

% fig2 = figure('Name','Fig2','NumberTitle','off');
% set(gcf,'Position',[35+0 35+0 900 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
% t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

nexttile
plot(tempCalciumPulse,'r');
hold on
% plot(x_fit,y_fit,'b');
% hold on









