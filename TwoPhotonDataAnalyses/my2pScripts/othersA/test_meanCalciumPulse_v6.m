close all




%% part 2


temp_sum = sum(tempCalciumPulse,2);
temptempBoolIndex = temp_sum~=0;
temptempIndex = find(temptempBoolIndex==true);

tempCalciumPulse_valid = tempCalciumPulse(temptempIndex,:);
tempCalciumPulse_length = zeros(length(temptempIndex),1);
for tempkk=1:size(tempCalciumPulse_valid,1)
    tempCalciumPulse_length(tempkk) = find(tempCalciumPulse_valid(tempkk,:)==0,1);
end
% temp_validLength = round(median(tempCalciumPulse_length));



% temp_validLength = min(tempCalciumPulse_length);





% tempCalciumPulse_mean = mean(tempCalciumPulse(temptempIndex,1:temp_validLength));
tempCalciumPulse_mean = mean(tempCalciumPulse(temptempIndex,:));
% tempCalciumPulse_mean = median(tempCalciumPulse(temptempIndex,:));

% tempCalciumPulse_mean = median(tempCalciumPulse(temptempIndex,1:temp_validLength));
% tempCalciumPulse_mean = mean(tempCalciumPulse(temptempIndex,:));
% tempCalciumPulse_mean = tempCalciumPulse_mean(1:find(tempCalciumPulse_mean~=0,1,'last'));


tempCalciumPulse_mean_smooth = smoothdata(tempCalciumPulse_mean,2,'gaussian',windowSize);
% tempCalciumPulse_mean_smooth = tempCalciumPulse_mean;
% tempCalciumPulse_mean_smooth = smoothdata(tempCalciumPulse_mean,2,'movmean',windowSize);

% tempCalciumPulse_mean = smoothdata(tempCalciumPulse_mean,2,'gaussian',5);

[temp_min2, temp_max2] = bounds(tempCalciumPulse_mean_smooth);
[valleys2,locs2] = findpeaks(temp_max2-tempCalciumPulse_mean_smooth); % Time-consuming!!!

tempCalciumPulse_mean_raw = tempCalciumPulse_mean;
if isempty(locs2) == false
% if false
    temp_validLength2 = locs2(1);
    %tempCalciumPulse_mean = tempCalciumPulse_mean(1:locs2(1));
else
    temp_validLength2 = 0;
end

temp_validLength_final = max(temp_validLength,temp_validLength2);
tempCalciumPulse_mean = tempCalciumPulse_mean(1:temp_validLength_final);
tempCalciumPulse_mean_smooth = tempCalciumPulse_mean_smooth(1:temp_validLength_final);


x = 1:length(tempCalciumPulse_mean);
y = tempCalciumPulse_mean;
% cftool(x,y);

[xData, yData] = prepareCurveData(x,double(y));
% ft = fittype( 'exp1' );
ft = fittype('a*(b^x)+c',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a','b','c'});
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 min(y)-100];
opts.Upper = [5000 1 min(y)+100];
% opts.StartPoint = [1.55650486517981 -0.0909553835273515 0];
opts.StartPoint = [100 0 0];
[fitresult, gof] = fit(xData,yData,ft,opts);

temp_r2_calciumPulse = gof.rsquare;
% temp_calciumPulse_decayCoeff = exp(fitresult.b);
temp_calciumPulse_decayCoeff = fitresult.b;

% temp_calciumPulse_decayCoeff = 0.7;


x_fit = 1:0.01:length(tempCalciumPulse_mean);
% y_fit = fitresult.a.*exp(fitresult.b.*x_fit);
y_fit = fitresult.a.*(temp_calciumPulse_decayCoeff.^x_fit) + fitresult.c;
% y_fit = fitresult.a.*(temp_calciumPulse_decayCoeff.^x_fit);



fig1 = figure('Name','Fig1','NumberTitle','off');
set(gcf,'Position',[35+0 35+0 900 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

nexttile
plot(tempCalciumPulse_mean,'r');
hold on
plot(x_fit,y_fit,'b');
hold on

nexttile
plot(tempCalciumPulse_mean_smooth,'r');
hold on
for temptempi=1:length(locs2)
   plot(locs2(temptempi)*[1 1],[temp_min2, temp_max2],'k'); 
   hold on
end

