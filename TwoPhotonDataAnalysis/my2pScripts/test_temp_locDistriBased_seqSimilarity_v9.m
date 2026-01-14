%% Initialization
% close all

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)


spatialInterferenceWeight = 0.5;%0.5

valid_length = 3;

if_precision_meanProb0_sumProb1 = 0;

Posterior_2d_model;
Posterior_2d_n11n_mean;



%% score_stimuli_to_response
fun_seqSimilarity_Name_v = autoGetFunName_myScripts('fun_seqSimilarity', [targetPATH '\functions']);
fun_seqSimilarity = str2func(fun_seqSimilarity_Name_v);

% score_stimuli_to_response = nan(sum(numSeq(1:valid_length)),sum(numSeq(1:valid_length)));
score_stimuli_to_response = nan(sum(numSeq(1:4)),sum(numSeq(1:4)));
for tempi=1:sum(numSeq(1:4))
    boolIndex_location_A = boolIndex_location_seq_T(tempi,:);
    
    for tempj=1:sum(numSeq(1:4))
        boolIndex_location_B = boolIndex_location_seq_T(tempj,:);
        
        score_stimuli_to_response(tempi,tempj) = ...
            fun_seqSimilarity(boolIndex_location_A,boolIndex_location_B,spatialInterferenceWeight);
        
        if tempi == 24 && tempj == 23
            a = 1; %#ok<*NASGU>
        end
    end
end

% temp_score = score_stimuli_to_response(1,:);
%
% temp_score_unique = sort(unique(temp_score),'descend');

%% seqProb_stimuli_to_response_n11n
Posterior_2d_n11n_mean;

seqProb_stimuli_to_response = nan(sum(numSeq(1:valid_length)),sum(numSeq(1:valid_length)));
seqProb_stimuli_to_response_n11n = nan(sum(numSeq(1:valid_length)),sum(numSeq(1:valid_length)));
for tempi=1:sum(numSeq(1:valid_length))
    temp_locDistri = Posterior_2d_n11n_mean(tempi,:);
    
    if isnan(temp_locDistri(1))
        continue
    end
    for tempj=1:sum(numSeq(1:valid_length))
        temp_boolIndex_location_seq = boolIndex_location_seq_T(tempj,:);
        
        temp_p = temp_locDistri;
        temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
        temp_p_seq = prod(temp_p,2);
        
        seqProb_stimuli_to_response(tempi,tempj) = temp_p_seq;
    end
    seqProb_stimuli_to_response_n11n(tempi,:) = seqProb_stimuli_to_response(tempi,:)./sum(seqProb_stimuli_to_response(tempi,:));
end

%% seqPrecision_neuron
score_stimuli_to_response;
seqProb_stimuli_to_response_n11n;

%     temp_score = score_stimuli_to_response(1,:);
%     % temp_score_unique = sort(unique(temp_score),'descend');
%     temp_seqProb = seqProb_stimuli_to_response_n11n(1,:);

seqSigma = nan(sum(numSeq(1:valid_length)),1);
for tempi=1:sum(numSeq(1:valid_length))
    temp_score = score_stimuli_to_response(tempi,1:sum(numSeq(1:valid_length)));
    temp_seqProb = seqProb_stimuli_to_response_n11n(tempi,:);
    
    if isnan(temp_seqProb(1))
        continue
    end
    
    %     x = temp_score';
    %     y = temp_seqProb';
    %     % cftool(x,y);
    %     f = fit(x,y,'gauss1');
    %     temp_c1 = f.c1;
    %     temp_sigma = temp_c1 / sqrt(2);
    
    if if_precision_meanProb0_sumProb1 == 0
        x = temp_score';
        y = temp_seqProb';
        % cftool(x,y);
        
    elseif if_precision_meanProb0_sumProb1 == 1
        temp_score_unique = unique(temp_score);
        temp_seqProb_uniqueSum = nan(size(temp_score_unique));
        for tempj=1:length(temp_score_unique)
            temptempBoolIndex = temp_score == temp_score_unique(tempj);
            
            temp_seqProb_uniqueSum(tempj) = sum(temp_seqProb(temptempBoolIndex));
        end
        x = temp_score_unique';
        y = temp_seqProb_uniqueSum';
        
    end
    
    
    mu = max(x);
    
    [xData, yData] = prepareCurveData(x,double(y));
    
    temp_str = sprintf('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
    
    %ft = fittype('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-6)/sigma)^2))',...
    %   'dependent',{'y'},'independent',{'x'},...
    %   'coefficients',{'sigma'});
    ft = fittype(temp_str,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'sigma'});
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = 0.5;
    [fitresult, gof] = fit(xData,yData,ft,opts);
    
    temp_r2 = gof.rsquare;
    temp_sigma = fitresult.sigma;
    
    
    
    seqSigma(tempi) = temp_sigma;
    
    if tempi == 25
        a = 1;
    end
end

seqPrecision_neuron = 1./seqSigma;


%% seqPrecision_behavior
score_stimuli_to_response;

responseMatrix_noChoice;
responseMatrix_noChoice_valid = responseMatrix_noChoice(1:sum(numSeq(1:valid_length)),1:sum(numSeq(1:valid_length)));
responseMatrix_noChoice_valid_n11n = responseMatrix_noChoice_valid ./ sum(responseMatrix_noChoice_valid,2);

seqSigma = nan(sum(numSeq(1:valid_length)),1);
for tempi=1:sum(numSeq(1:valid_length))
    temp_score = score_stimuli_to_response(tempi,1:sum(numSeq(1:valid_length)));
    temp_seqProb = responseMatrix_noChoice_valid_n11n(tempi,:);
    
    if isnan(temp_seqProb(1))
        continue
    end
    
    if if_precision_meanProb0_sumProb1 == 0
        x = temp_score';
        y = temp_seqProb';
        % cftool(x,y);
        
    elseif if_precision_meanProb0_sumProb1 == 1
        temp_score_unique = unique(temp_score);
        temp_seqProb_uniqueSum = nan(size(temp_score_unique));
        for tempj=1:length(temp_score_unique)
            temptempBoolIndex = temp_score == temp_score_unique(tempj);
            
            temp_seqProb_uniqueSum(tempj) = sum(temp_seqProb(temptempBoolIndex));
        end
        x = temp_score_unique';
        y = temp_seqProb_uniqueSum';
        
    end
    
    
    mu = max(x);
    
    [xData, yData] = prepareCurveData(x,double(y));
    
    temp_str = sprintf('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
    
    %ft = fittype('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-6)/sigma)^2))',...
    %   'dependent',{'y'},'independent',{'x'},...
    %   'coefficients',{'sigma'});
    ft = fittype(temp_str,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'sigma'});
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = 0.5;
    [fitresult, gof] = fit(xData,yData,ft,opts);
    
    temp_r2 = gof.rsquare;
    temp_sigma = fitresult.sigma;
    
    seqSigma(tempi) = temp_sigma;
    
    if tempi == 25
        a = 1;
    end
end

seqPrecision_behavior = 1./seqSigma;

%%
x = seqPrecision_neuron(~isnan(seqPrecision_neuron));
y = seqPrecision_behavior(1:sum(numSeq(1:valid_length)));
y = y(~isnan(seqPrecision_neuron));
[r_123,p_123] = corr(x,y);



%% Plot
if_plot_seqCorr = 1;
if if_plot_seqCorr == 1
    %% fig, Sequence decoding
    fig = figure('Name','Sequence decoding','NumberTitle','off'); %#ok<*NASGU>
    set(gcf,'Position',[50+0 50+0 240 252*0.925]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    temp_x = seqPrecision_neuron;
    temp_y = seqPrecision_behavior;
    tempBoolIndex = ~isnan(temp_x);

    x = temp_x(~isnan(temp_x));
    y = temp_y(1:sum(numSeq(1:valid_length)));
    y = y(~isnan(temp_x));
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);

    
    
    h = [];
    for tempi=1:valid_length
        temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
        tempIndex = find(tempBoolIndex(temp_range)==true);
        temp_range2 = temp_range(tempIndex);
        
        temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
        temp_h = scatter(temp_x(temp_range2), temp_y(temp_range2), ...
            temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        h = [h temp_h];
    end
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y);
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    h = [h temp_h];
        
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.1 temp_xmax+(temp_xmax-temp_xmin)*0.1]);
    ylim([temp_ymin-(temp_ymax-temp_ymin)*0.1 temp_ymax+(temp_ymax-temp_ymin)*0.1]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 12)
    set(gca,'box','off');% 取消右、上边框
    %set(gca,'XTick',0:0.5:1,'FontSize', 12);%给坐标加标签
    %set(gca,'YTick',0:0.5:1,'FontSize', 12);%给坐标加标签
    
    xlabel(sprintf('Neuron precision'), 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('%s\n%d seqs, r=%.3f, p=%.3f',currentSession,length(y),r,p), 'FontSize', 12, 'FontWeight', 'bold');
    temp_title.Interpreter = 'none';
    
    ylabel(sprintf('Behavior precision'), 'FontSize', 12, 'FontWeight', 'bold');

end


%% End