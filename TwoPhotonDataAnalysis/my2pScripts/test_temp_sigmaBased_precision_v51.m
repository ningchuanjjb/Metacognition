% Chuan's 6th script (20251214)
% This script: One part of 'stepD3_train123_test123.m'.
%% Initialization
close all

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)


if_pProd_seqLevel_n11n = 1;

% spatialInterferenceWeight = 0.5;%0.5

valid_length = 3;

if_precision_meanProb0_sumProb1 = 0;


if_memoryPrecision_accuracy0_sigma1 = if_memoryPrecision_accuracy0_sigma1; %#ok<*ASGSL>

if_behavior_sameWithNeuron = 1;%0

if if_behavior_sameWithNeuron == 1
    if_memoryPrecision_accuracy0_sigma1_forBehavior = if_memoryPrecision_accuracy0_sigma1;
elseif if_behavior_sameWithNeuron == 0
    if_memoryPrecision_accuracy0_sigma1_forBehavior = 0;
end


if_entropy = 1;
if_entropy_behaviorInverted = 0;

% Posterior_2d_model;
Posterior_2d_n11n_mean;

if_fit_gain = 1;

if_plot_beavior = 1;


exampleSeq = 8;%17,8


responseMatrix_noChoice;
responseMatrix_noChoice_valid = responseMatrix_noChoice(1:sum(numSeq(1:4)),1:sum(numSeq(1:4)));
responseMatrix_noChoice_valid_n11n = responseMatrix_noChoice_valid ./ sum(responseMatrix_noChoice_valid,2);

%% To get optimal spatialInterferenceWeight
fun_seqSimilarity_Name_v = autoGetFunName_myScripts('fun_seqSimilarity', [targetPATH '\functions']);
fun_seqSimilarity = str2func(fun_seqSimilarity_Name_v);



temp_spatialInterferenceWeight = 0:0.005:1;
temp_r_mean_iter = nan(length(temp_spatialInterferenceWeight),1);
for tempIter=1:length(temp_spatialInterferenceWeight)
    temp_score_stimuli_to_response = nan(sum(numSeq(1:4)),sum(numSeq(1:4)));
    for tempi=1:sum(numSeq(1:4))
        boolIndex_location_A = boolIndex_location_seq_T(tempi,:);
        
        for tempj=1:sum(numSeq(1:4))
            boolIndex_location_B = boolIndex_location_seq_T(tempj,:);
            
            temp_score_stimuli_to_response(tempi,tempj) = ...
                fun_seqSimilarity(boolIndex_location_A,boolIndex_location_B,temp_spatialInterferenceWeight(tempIter));
            
            if tempi == 24 && tempj == 23
                a = 1; %#ok<*NASGU>
            end
        end
    end
    
    temp_r = nan(size(responseMatrix_noChoice_valid_n11n,1),1);
    for tempi=1:length(temp_r)
        [temp_r(tempi),~] = corr(responseMatrix_noChoice_valid_n11n(tempi,:)',temp_score_stimuli_to_response(tempi,:)');
    end
    temp_r_mean_iter(tempIter) = mean(temp_r);
end
[M,I] = max(temp_r_mean_iter);
temp_r_max_iter = M;
spatialInterferenceWeight = temp_spatialInterferenceWeight(I);
fprintf('spatialInterferenceWeight=%.4f,r=%.4f.\n',spatialInterferenceWeight,temp_r_max_iter);

%% score_stimuli_to_response
score_stimuli_to_response = nan(sum(numSeq(1:4)),sum(numSeq(1:4)));
for tempi=1:sum(numSeq(1:4))
    boolIndex_location_A = boolIndex_location_seq_T(tempi,:);
    for tempj=1:sum(numSeq(1:4))
        %if tempi == 24 && tempj == 23
        if tempi == 1 && tempj == 6
            a = 1; %#ok<*NASGU>
        end
        boolIndex_location_B = boolIndex_location_seq_T(tempj,:);
        score_stimuli_to_response(tempi,tempj) = ...
            fun_seqSimilarity(boolIndex_location_A,boolIndex_location_B,spatialInterferenceWeight);
        
    end
end


% % score_stimuli_to_response = nan(sum(numSeq(1:valid_length)),sum(numSeq(1:valid_length)));
% score_stimuli_to_response = nan(sum(numSeq(1:4)),sum(numSeq(1:4)));
% for tempi=1:sum(numSeq(1:4))
%     boolIndex_location_A = boolIndex_location_seq_T(tempi,:);
%
%     for tempj=1:sum(numSeq(1:4))
%         boolIndex_location_B = boolIndex_location_seq_T(tempj,:);
%
%         score_stimuli_to_response(tempi,tempj) = ...
%             fun_seqSimilarity(boolIndex_location_A,boolIndex_location_B,spatialInterferenceWeight);
%
%         if tempi == 24 && tempj == 23
%             a = 1; %#ok<*NASGU>
%         end
%     end
% end

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
if if_memoryPrecision_accuracy0_sigma1 == 1
    
    if if_entropy == 0
        
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
            
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            
            if if_fit_gain == 0
                temp_str = sprintf('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
                ft = fittype(temp_str,...
                    'dependent',{'y'},'independent',{'x'},...
                    'coefficients',{'sigma'});
                opts.StartPoint = 0.5;
            elseif if_fit_gain == 1
                temp_str = sprintf('gain*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
                ft = fittype(temp_str,...
                    'dependent',{'y'},'independent',{'x'},...
                    'coefficients',{'sigma','gain'});
                opts.StartPoint = [0.5 0.5];
                opts.Lower = [0 0];
            end
            [fitresult, gof] = fit(xData,yData,ft,opts);
            
            temp_r2 = gof.rsquare;
            temp_sigma = fitresult.sigma;
            
            
            
            seqSigma(tempi) = abs(temp_sigma);
            
            %if tempi == 21
            %if tempi == 17
            if tempi == 3
                a = 1;
            end
        end
        
        seqPrecision_neuron = 1./seqSigma;
        % seqPrecision_neuron = 1./(seqSigma.^2);
        % seqPrecision_neuron = seqSigma;
        
    elseif if_entropy == 1
        
        p = ones(1,sum(numSeq(1:valid_length)))./sum(numSeq(1:valid_length));
        p = p + eps;
        p = p./sum(p);
        entropy_max = -sum(p.*log2(p));
        
        seqPrecision_neuron = nan(sum(numSeq(1:valid_length)),1);
        for tempi=1:sum(numSeq(1:valid_length))
            if if_pProd_seqLevel_n11n == 1
                temp_seqProb = seqProb_stimuli_to_response_n11n(tempi,:);
            elseif if_pProd_seqLevel_n11n == 0
                temp_seqProb = seqProb_stimuli_to_response(tempi,:);
            end
            if isnan(temp_seqProb(1))
                continue
            end
            %seqPrecision_neuron(tempi) = temp_seqProb(tempi);
            
            %p = [1 0 0 0 0 0];
            %p = [1 1 1 1 1 1]./6;
            p = temp_seqProb;
            p = p + eps;
            p = p./sum(p);
            entropy = -sum(p.*log2(p));
            
            %seqPrecision_neuron(tempi) = 1./entropy;
            %seqPrecision_neuron(tempi) = entropy_max-entropy;
            seqPrecision_neuron(tempi) = (entropy_max-entropy)./entropy_max;
        end
        
    end
    
    
elseif if_memoryPrecision_accuracy0_sigma1 == 0
    seqPrecision_neuron = nan(sum(numSeq(1:valid_length)),1);
    for tempi=1:sum(numSeq(1:valid_length))
        if if_pProd_seqLevel_n11n == 1
            temp_seqProb = seqProb_stimuli_to_response_n11n(tempi,:);
        elseif if_pProd_seqLevel_n11n == 0
            temp_seqProb = seqProb_stimuli_to_response(tempi,:);
        end
        if isnan(temp_seqProb(1))
            continue
        end
        seqPrecision_neuron(tempi) = temp_seqProb(tempi);
    end
    
end


%% seqPrecision_behavior
temp_valid_length = 3;

score_stimuli_to_response;

responseMatrix_noChoice;
responseMatrix_noChoice_valid = responseMatrix_noChoice(1:sum(numSeq(1:temp_valid_length)),1:sum(numSeq(1:temp_valid_length)));
responseMatrix_noChoice_valid_n11n = responseMatrix_noChoice_valid ./ sum(responseMatrix_noChoice_valid,2);

% if if_memoryPrecision_accuracy0_sigma1 == 1
if if_memoryPrecision_accuracy0_sigma1_forBehavior == 1
    % if false
    if if_entropy == 0
        
        seqSigma = nan(sum(numSeq(1:temp_valid_length)),1);
        for tempi=1:sum(numSeq(1:temp_valid_length))
            temp_score = score_stimuli_to_response(tempi,1:sum(numSeq(1:temp_valid_length)));
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
            
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            
            %     if if_fit_gain == 0
            temp_str = sprintf('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
            ft = fittype(temp_str,...
                'dependent',{'y'},'independent',{'x'},...
                'coefficients',{'sigma'});
            opts.StartPoint = 0.5;
            %     elseif if_fit_gain == 1
            %         temp_str = sprintf('gain*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
            %         ft = fittype(temp_str,...
            %             'dependent',{'y'},'independent',{'x'},...
            %             'coefficients',{'sigma','gain'});
            %         opts.StartPoint = [0.5 0.5];
            %         opts.Lower = [0 0];
            %     end
            [fitresult, gof] = fit(xData,yData,ft,opts);
            
            %     temp_str = sprintf('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
            %     %temp_str = sprintf('gain*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
            %
            %     %ft = fittype('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-6)/sigma)^2))',...
            %     %   'dependent',{'y'},'independent',{'x'},...
            %     %   'coefficients',{'sigma'});
            %     ft = fittype(temp_str,...
            %        'dependent',{'y'},'independent',{'x'},...
            %        'coefficients',{'sigma'});
            %     %ft = fittype(temp_str,...
            %     %    'dependent',{'y'},'independent',{'x'},...
            %     %    'coefficients',{'sigma','gain'});
            %     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            %     opts.Display = 'Off';
            %     opts.StartPoint = 0.5;
            %     %opts.StartPoint = [0.5 0.5];
            %     %opts.Lower = [0.1 0];
            %     [fitresult, gof] = fit(xData,yData,ft,opts);
            
            temp_r2 = gof.rsquare;
            temp_sigma = fitresult.sigma;
            
            seqSigma(tempi) = abs(temp_sigma);
            
            if tempi == 15
                a = 1;
            end
        end
        
        seqPrecision_behavior = 1./seqSigma;
        % seqPrecision_behavior_56 = seqPrecision_behavior;
        % seqPrecision_behavior = seqPrecision_behavior(1:sum(numSeq(1:valid_length)));
        
        a = 1;
        
        % seqPrecision_behavior = 1./(seqSigma.^2);
        % seqPrecision_behavior = seqSigma;
        
    elseif if_entropy == 1
        
        p = ones(1,sum(numSeq(1:valid_length)))./sum(numSeq(1:valid_length));
        p = p + eps;
        p = p./sum(p);
        entropy_max = -sum(p.*log2(p));
        
        seqPrecision_behavior = nan(sum(numSeq(1:valid_length)),1);
        for tempi=1:sum(numSeq(1:valid_length))
            temp_seqProb = responseMatrix_noChoice_valid_n11n(tempi,:);
            
            %seqPrecision_neuron(tempi) = temp_seqProb(tempi);
            
            %p = [1 0 0 0 0 0];
            %p = [1 1 1 1 1 1]./6;
            p = temp_seqProb;
            p = p + eps;
            p = p./sum(p);
            entropy = -sum(p.*log2(p));
            
            %seqPrecision_behavior(tempi) = 1./entropy;
            %seqPrecision_behavior(tempi) = entropy_max-entropy;
            if if_entropy_behaviorInverted == 1
                seqPrecision_behavior(tempi) = (entropy_max-entropy)./entropy_max;
            elseif if_entropy_behaviorInverted == 0
                seqPrecision_behavior(tempi) = entropy./entropy_max;
            end
            
        end
        
        
    end
    
    % end
    
    % elseif if_memoryPrecision_accuracy0_sigma1 == 0
elseif if_memoryPrecision_accuracy0_sigma1_forBehavior == 0
    % if true
    %     seqPrecision_behavior = nan(sum(numSeq(1:valid_length)),1);
    %     for tempi=1:sum(numSeq(1:valid_length))
    %         temp_seqProb = responseMatrix_noChoice_valid_n11n(tempi,:);
    %         %temp_seqProb = responseMatrix_noChoice_valid(tempi,:);
    %         if isnan(temp_seqProb(1))
    %             continue
    %         end
    %         seqPrecision_behavior(tempi) = temp_seqProb(tempi);
    %     end
    
    
    if if_gAcc_behavior0_model1 == 0
        seqPrecision_behavior = gAcc_target_collapsed_inOne(1:sum(numSeq(1:valid_length)))';
        %seqPrecision_behavior = gAcc_noChoice_collapsed_inOne(1:sum(numSeq(1:valid_length)))';
    elseif if_gAcc_behavior0_model1 == 1
        seqPrecision_behavior = gAcc_noChoice_inOne_model(1:sum(numSeq(1:valid_length)))';
    end
    
end


%% seqPrecision_behavior_56
temp_valid_length = 4;

score_stimuli_to_response;

responseMatrix_noChoice;
responseMatrix_noChoice_valid = responseMatrix_noChoice(1:sum(numSeq(1:temp_valid_length)),1:sum(numSeq(1:temp_valid_length)));
responseMatrix_noChoice_valid_n11n = responseMatrix_noChoice_valid ./ sum(responseMatrix_noChoice_valid,2);

% if if_memoryPrecision_accuracy0_sigma1 == 1
if if_memoryPrecision_accuracy0_sigma1_forBehavior == 1
    
    % if false
    if if_entropy == 0
        
        seqSigma = nan(sum(numSeq(1:temp_valid_length)),1);
        for tempi=1:sum(numSeq(1:temp_valid_length))
            temp_score = score_stimuli_to_response(tempi,1:sum(numSeq(1:temp_valid_length)));
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
            
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            
            %     if if_fit_gain == 0
            temp_str = sprintf('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
            ft = fittype(temp_str,...
                'dependent',{'y'},'independent',{'x'},...
                'coefficients',{'sigma'});
            opts.StartPoint = 0.5;
            %     elseif if_fit_gain == 1
            %         temp_str = sprintf('gain*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
            %         ft = fittype(temp_str,...
            %             'dependent',{'y'},'independent',{'x'},...
            %             'coefficients',{'sigma','gain'});
            %         opts.StartPoint = [0.5 0.5];
            %         opts.Lower = [0 0];
            %     end
            [fitresult, gof] = fit(xData,yData,ft,opts);
            
            %     temp_str = sprintf('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
            %     %temp_str = sprintf('gain*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
            %
            %     %ft = fittype('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-6)/sigma)^2))',...
            %     %   'dependent',{'y'},'independent',{'x'},...
            %     %   'coefficients',{'sigma'});
            %     ft = fittype(temp_str,...
            %        'dependent',{'y'},'independent',{'x'},...
            %        'coefficients',{'sigma'});
            %     %ft = fittype(temp_str,...
            %     %    'dependent',{'y'},'independent',{'x'},...
            %     %    'coefficients',{'sigma','gain'});
            %     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            %     opts.Display = 'Off';
            %     opts.StartPoint = 0.5;
            %     %opts.StartPoint = [0.5 0.5];
            %     %opts.Lower = [0.1 0];
            %     [fitresult, gof] = fit(xData,yData,ft,opts);
            
            temp_r2 = gof.rsquare;
            temp_sigma = fitresult.sigma;
            
            seqSigma(tempi) = abs(temp_sigma);
            
            if tempi == 15
                a = 1;
            end
        end
        
        seqPrecision_behavior_56 = 1./seqSigma;
        % seqPrecision_behavior_56 = seqPrecision_behavior;
        % seqPrecision_behavior = seqPrecision_behavior(1:sum(numSeq(1:valid_length)));
        
        
    elseif if_entropy == 1
        
        p = ones(1,sum(numSeq(1:temp_valid_length)))./sum(numSeq(1:temp_valid_length));
        p = p + eps;
        p = p./sum(p);
        entropy_max = -sum(p.*log2(p));
        
        
        seqPrecision_behavior_56 = nan(sum(numSeq(1:temp_valid_length)),1);
        for tempi=1:sum(numSeq(1:temp_valid_length))
            temp_seqProb = responseMatrix_noChoice_valid_n11n(tempi,:);
            
            %p = [1 0 0 0 0 0];
            %p = [1 1 1 1 1 1]./6;
            p = temp_seqProb;
            p = p + eps;
            p = p./sum(p);
            entropy = -sum(p.*log2(p));
            
            %seqPrecision_behavior_56(tempi) = 1./entropy;
            %seqPrecision_behavior_56(tempi) = entropy_max-entropy;
            if if_entropy_behaviorInverted == 1
                seqPrecision_behavior_56(tempi) = (entropy_max-entropy)./entropy_max;
            elseif if_entropy_behaviorInverted == 0
                seqPrecision_behavior_56(tempi) = entropy./entropy_max;
            end
            
        end
        
    end
    
    % end
    
    % elseif if_memoryPrecision_accuracy0_sigma1 == 0
elseif if_memoryPrecision_accuracy0_sigma1_forBehavior == 0
    % if true
    %     seqPrecision_behavior_56 = nan(sum(numSeq(1:temp_valid_length)),1);
    %     for tempi=1:sum(numSeq(1:temp_valid_length))
    %         temp_seqProb = responseMatrix_noChoice_valid_n11n(tempi,:);
    %         %temp_seqProb = responseMatrix_noChoice_valid(tempi,:);
    %         if isnan(temp_seqProb(1))
    %             continue
    %         end
    %         seqPrecision_behavior_56(tempi) = temp_seqProb(tempi);
    %     end
    
    if if_gAcc_behavior0_model1 == 0
        seqPrecision_behavior_56 = gAcc_target_collapsed_inOne';
        %seqPrecision_behavior_56 = gAcc_noChoice_collapsed_inOne';
    elseif if_gAcc_behavior0_model1 == 1
        seqPrecision_behavior_56 = gAcc_noChoice_inOne_model';
    end
    
end


%%
x = seqPrecision_neuron(~isnan(seqPrecision_neuron));
y = seqPrecision_behavior(1:sum(numSeq(1:valid_length)));
y = y(~isnan(seqPrecision_neuron));
[r_123,p_123] = corr(x,y);



%% Plot
if_plot_seqCorr = 1;
if if_plot_seqCorr == 1
    if if_plot_beavior == 1
        %% , Sequence decoding
        fig = figure('Name','Sequence decoding','NumberTitle','off'); %#ok<*NASGU>
        %set(gcf,'Position',[50+0 50+0 240 252*0.925]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 240*1.45 233*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[50+0 50+0 226 186]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        temp_x = seqPrecision_behavior;
        temp_y = offloadingProb_inOne(1:sum(numSeq(1:valid_length)))';
        
        %temp_x = seqPrecision_behavior_56;
        %temp_y = offloadingProb_inOne';
        
        tempBoolIndex = ~isnan(temp_x);
        
        x = temp_x(~isnan(temp_x));
        %y = temp_y(1:sum(numSeq(1:valid_length)));
        y = temp_y;        
        y = y(~isnan(temp_x));
        
        [r,p] = corr(x,y);
        
        mdl = fitglm(x,y);
        
        
        
        h = [];
        for tempi=1:valid_length
            temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
            tempIndex = find(tempBoolIndex(temp_range)==true);
            temp_range2 = temp_range(tempIndex);
            
            %temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
            temp_size = 10;
            temp_h = scatter(temp_x(temp_range2), temp_y(temp_range2), ...
                temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
            h = [h temp_h]; %#ok<*AGROW>
        end
        
        [temp_xmin,temp_xmax] = bounds(x);
        [temp_ymin,temp_ymax] = bounds(y); %#ok<*ASGLU>
        
        x_fit = temp_xmin:0.001:temp_xmax;
        y_fit = predict(mdl,x_fit')';
        
        temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        hold on
        h = [h temp_h];
        
        xlim([temp_xmin-(temp_xmax-temp_xmin)*0.1 temp_xmax+(temp_xmax-temp_xmin)*0.1]);
        %ylim([temp_ymin-(temp_ymax-temp_ymin)*0.1 temp_ymax+(temp_ymax-temp_ymin)*0.1]);
        %ylim([0 temp_ymax+(temp_ymax-temp_ymin)*0.1]);
        ylim([0 1]);
        
        set(gca,'linewidth',1.5)
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        %set(gca,'XTick',0:0.5:1,'FontSize', 12);%给坐标加标签
        %set(gca,'YTick',0:0.5:1,'FontSize', 12);%给坐标加标签
        set(gca,'YTick',0:1);%给坐标加标签
        
        if if_memoryPrecision_accuracy0_sigma1_forBehavior == 0
            xlabel(sprintf('Behavior accuracy'), 'FontSize', 9);
        elseif if_memoryPrecision_accuracy0_sigma1_forBehavior == 1
            xlabel(sprintf('Recall variability'), 'FontSize', 9);
        end
        %temp_title = title(sprintf('%s\n%d seqs, r=%.3f, p=%.3f',currentSession,length(y),r,p), 'FontSize', 9);
        temp_title = title(sprintf('r=%.3f, p=%.3f',r,p), 'FontSize', 9);
        temp_title.Interpreter = 'none';
        
        ylabel(sprintf('Offloading rate'), 'FontSize', 9);
        
    end
    
    %% fig, Sequence decoding
    fig = figure('Name','Sequence decoding','NumberTitle','off'); %#ok<*NASGU>
    %set(gcf,'Position',[50+0 50+0 240 252*0.925]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 240*1.45 233*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 226 186]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 226*0.75 186]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %     if if_memoryPrecision_accuracy0_sigma1 == 0
    %set(gcf,'Position',[50+0 50+0 226*0.75*0.78 186*0.78]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.75 336*1.11*0.5*0.7*0.85*0.95*1.78*0.85*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 226*0.75*0.78*0.94 186*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[300+0 50+0 226*0.75*0.78*0.94*1.25 186*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    %     elseif if_memoryPrecision_accuracy0_sigma1 == 1
    %         set(gcf,'Position',[300+0 50+0 226*0.75 186]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %     end
    
    %t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    %nexttile
    
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
        
        %temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
        temp_size = 10;
        temp_h = scatter(temp_x(temp_range2), temp_y(temp_range2), ...
            temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        h = [h temp_h];
    end
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y);
    if if_memoryPrecision_accuracy0_sigma1 == 0
        temp_xmin = 0;
        temp_xmax = 1;
        temp_ymin = 0;
        temp_ymax = 1;
    elseif if_memoryPrecision_accuracy0_sigma1 == 1
        temp_xmin = 0;
        temp_xmax = 1;
        temp_ymin = 0;
        temp_ymax = 1;
    end
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    h = [h temp_h];
    
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    %     if if_memoryPrecision_accuracy0_sigma1 == 0
    %         %text(temp_xmin+(temp_xmax-temp_xmin)*0.45,0.22,sprintf('r = %.3f',r), 'fontsize',9,'FontWeight','normal');
    %         %text(temp_xmin+(temp_xmax-temp_xmin)*0.45,0.07,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
    text(temp_xmin+(temp_xmax-temp_xmin)*0.55,0.22,sprintf('r = %.3f',r), 'fontsize',9,'FontWeight','normal');
    text(temp_xmin+(temp_xmax-temp_xmin)*0.55,0.07,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
    %     elseif if_memoryPrecision_accuracy0_sigma1 == 1
    %         text(temp_xmin+(temp_xmax-temp_xmin)*0.58,temp_ymin+(temp_ymax-temp_ymin)*0.17,sprintf('r=%.3f',r), 'fontsize',9,'FontWeight','normal');
    %         text(temp_xmin+(temp_xmax-temp_xmin)*0.58,temp_ymin+(temp_ymax-temp_ymin)*0.04,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
    %     end
    
    
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.1 temp_xmax+(temp_xmax-temp_xmin)*0.1]);
    ylim([temp_ymin-(temp_ymax-temp_ymin)*0.1 temp_ymax+(temp_ymax-temp_ymin)*0.1]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    %set(gca,'XTick',0:0.5:1,'FontSize', 12);%给坐标加标签
    %set(gca,'YTick',0:0.5:1,'FontSize', 12);%给坐标加标签
    
    %if if_memoryPrecision_accuracy0_sigma1 == 0
    %xticks([0 1]);
    %yticks([0 1]);
    
    xticks([0:0.2:1]);
    yticks([0:0.2:1]);
    
    xtickangle(0);
    
    %end
    
    %xlabel(sprintf('Neuron precision'), 'FontSize', 11, 'FontWeight', 'bold');
    %xlabel(sprintf('Memory precision'), 'FontSize', 9, 'FontWeight', 'normal');
    xlabel(sprintf('Memory quality'), 'FontSize', 9, 'FontWeight', 'normal');
    %temp_title = title(sprintf('%s\n%d seqs, r=%.3f, p=%.3f',currentSession,length(y),r,p), 'FontSize', 11, 'FontWeight', 'bold');
    %temp_title = title(sprintf('%s\nr=%.3f, p=%.3f',currentSession,r,p), 'FontSize', 11, 'FontWeight', 'bold');
    %temp_title = title(sprintf('%s',currentSession), 'FontSize', 11, 'FontWeight', 'bold');
    
    %     if if_memoryPrecision_accuracy0_sigma1 == 0
    %temp_title = title(sprintf('Sequence-level'), 'FontSize', 9, 'FontWeight', 'normal');
    temp_title = title(sprintf('Across seqs'), 'FontSize', 9, 'FontWeight', 'normal');
    %     elseif if_memoryPrecision_accuracy0_sigma1 == 1
    %         %temp_title = title(sprintf('Precision (1/sigma)'), 'FontSize', 10, 'FontWeight', 'normal');
    %         temp_title = title(sprintf('Neuronal correlation'), 'FontSize', 9, 'FontWeight', 'normal');
    %     end
    temp_title.Interpreter = 'none';
    
    if if_memoryPrecision_accuracy0_sigma1_forBehavior == 0
        ylabel(sprintf('Behavior accuracy'), 'FontSize', 9, 'FontWeight', 'normal');
    elseif if_memoryPrecision_accuracy0_sigma1_forBehavior == 1
        ylabel(sprintf('Recall variability'), 'FontSize', 9, 'FontWeight', 'normal');
        %ylabel(sprintf('Behavior accuracy'), 'FontSize', 9, 'FontWeight', 'normal');
    end
    
    
    %% fig, Sequence decoding
    fig = figure('Name','Sequence decoding','NumberTitle','off'); %#ok<*NASGU>
    %set(gcf,'Position',[50+0 50+0 240 252*0.925]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 240*1.45 233*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 226 186]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[550+0 50+0 226*0.75 186]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    temp_x = seqPrecision_neuron;
    %temp_y = seqPrecision_behavior;
    temp_y = offloadingProb_inOne(1:sum(numSeq(1:valid_length)))';
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
        
        %temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
        temp_size = 10;
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
    %ylim([temp_ymin-(temp_ymax-temp_ymin)*0.1 temp_ymax+(temp_ymax-temp_ymin)*0.1]);
    ylim([0 1]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    %set(gca,'XTick',0:0.5:1,'FontSize', 12);%给坐标加标签
    %set(gca,'YTick',0:0.5:1,'FontSize', 12);%给坐标加标签
    set(gca,'YTick',0:1);%给坐标加标签
    
    xlabel(sprintf('Memory quality'), 'FontSize', 9);
    %temp_title = title(sprintf('%s\n%d seqs, r=%.3f, p=%.3f',currentSession,length(y),r,p), 'FontSize', 11, 'FontWeight', 'bold');
    %temp_title = title(sprintf('%s\nr=%.3f, p=%.3f',currentSession,r,p), 'FontSize', 9);
    temp_title = title(sprintf('r=%.3f, p=%.3f',r,p), 'FontSize', 9);
    temp_title.Interpreter = 'none';
    
    %ylabel(sprintf('Behavior precision'), 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(sprintf('Offloading rate'), 'FontSize', 9);
    
    
    
    
    
    %% fig, Behavior model
    if if_locationDistri_behavior0_model1 == 1
        fig = figure('Name','Sequence decoding','NumberTitle','off'); %#ok<*NASGU>
        set(gcf,'Position',[800+0 50+0 226*0.75*1.2 186*1.15]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        temp_x = gAcc_noChoice_inOne_model(1:sum(numSeq(1:valid_length)))';
        temp_y = gAcc_target_collapsed_inOne(1:sum(numSeq(1:valid_length)))';
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
        
        
        
        tempTxt = 'p>0.05';
        if p < 0.05
            tempTxt = 'p<0.05';
        end
        if p < 0.01
            tempTxt = 'p<0.01';
        end
        if p < 0.001
            tempTxt = 'p<0.001';
        end
        if if_memoryPrecision_accuracy0_sigma1 == 0
            text(temp_xmin+(temp_xmax-temp_xmin)*0.58,0.18,sprintf('r=%.3f',r), 'fontsize',9,'FontWeight','normal');
            text(temp_xmin+(temp_xmax-temp_xmin)*0.58,0.07,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
        elseif if_memoryPrecision_accuracy0_sigma1 == 1
            text(temp_xmin+(temp_xmax-temp_xmin)*0.58,temp_ymin+(temp_ymax-temp_ymin)*0.17,sprintf('r=%.3f',r), 'fontsize',9,'FontWeight','normal');
            text(temp_xmin+(temp_xmax-temp_xmin)*0.58,temp_ymin+(temp_ymax-temp_ymin)*0.04,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
        end
        
        
        %xlim([temp_xmin-(temp_xmax-temp_xmin)*0.1 temp_xmax+(temp_xmax-temp_xmin)*0.1]);
        %ylim([temp_ymin-(temp_ymax-temp_ymin)*0.1 temp_ymax+(temp_ymax-temp_ymin)*0.1]);
        
        xlim([0 1]);
        ylim([0 1]);
        
        set(gca,'linewidth',1.5)
        set(gca, 'FontSize', 10)
        set(gca,'box','off');% 取消右、上边框
        %set(gca,'XTick',0:0.5:1,'FontSize', 12);%给坐标加标签
        %set(gca,'YTick',0:0.5:1,'FontSize', 12);%给坐标加标签
        
        set(gca,'XTick',0:1,'FontSize', 10);%给坐标加标签
        set(gca,'YTick',0:1,'FontSize', 10);%给坐标加标签
        
        %xlabel(sprintf('Behavior precision'), 'FontSize', 11, 'FontWeight', 'bold');
        xlabel(sprintf('Behavior accuracy'), 'FontSize', 10, 'FontWeight', 'normal');
        %temp_title = title(sprintf('%s\nr=%.3f, p=%.3f',currentSession,r,p), 'FontSize', 11, 'FontWeight', 'bold');
        %temp_title.Interpreter = 'none';
        
        %ylabel(sprintf('Behavior precision\n(model)'), 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(sprintf('Behavior accuracy (model)'), 'FontSize', 10, 'FontWeight', 'normal');
        
    end
    
end


%%
behaivorMatrix_model;
responseMatrix_noChoice;
responseMatrix_noChoice_n11n = responseMatrix_noChoice ./ sum(responseMatrix_noChoice,2);

C_min = 0;%0
C_max = 1;


% x1 = responseMatrix_noChoice_n11n;
% x2 = behaivorMatrix_model;
temp_coeff = 1;%0.8
x1 = responseMatrix_noChoice_n11n.^temp_coeff;
x2 = behaivorMatrix_model.^temp_coeff;


%% Plot A stimuli-to-response matrix of behavior data and behavior model
if false
    fig = figure('Name','responseMatrix','NumberTitle','off');
    % set(gcf,'Position',[10 50 533*1.1*1.018 224*1.1-3]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 50 597*0.982 243*2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 50 597*0.982*0.97 243*2*0.85*1.03]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    %t = tiledlayout(1,10,'TileSpacing','Compact','Padding','Loose');
    % t = tiledlayout(1,10,'TileSpacing','Loose','Padding','Loose');
    t = tiledlayout(2,10,'TileSpacing','Loose','Padding','Compact');
    %t = tiledlayout(2,10,'TileSpacing','Tight','Padding','Compact');
    
    my_gray = gray;
    my_gray = my_gray(end:-1:1,:);
    colormap(my_gray);
    
    %% behavior data
    %nexttile
    nexttile([1 5])
    
    C = x1;
    
    imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
    hold on
    
    
    % Plot length bound in stimuli
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([0.5 size(C,2)+0.5], [temptempi temptempi]-0.5, 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    % Plot length bound in response
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    plot(sum(numSeq).*[1 1]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
    hold on
    
    
    set(gca,'linewidth',1.5)
    %set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8) %12
    
    set(gca,'TickLength',[0 0]);
    set(gca,'XTick',1:size(C,1));
    set(gca,'YTick',1:size(C,2));
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    
    tempLabelStr_res = string(1:5);
    xtl=string(tempLabelStr_res);
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    %xtext_xp=xt;
    xtext_xp = nan(1,4);
    for tempi=1:length(xtext_xp)
        xtext_xp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    xtext_xp = [xtext_xp xtext_xp(end)+20/2];
    xtext_xp = xtext_xp+0.5;
    
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(tempLabelStr_res))+0.75-3.75;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_res))-1.6;%-0.75
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',10);%8
    set(gca,'xticklabel','');
    
    
    
    tempLabelStr_sti = string(1:4);
    ytl=string(tempLabelStr_sti);
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_sti))-2.4;%-0.75
    % 设置ttext的y坐标位置
    % ytext_yp=yt;
    ytext_yp = nan(1,4);
    for tempi=1:length(ytext_yp)
        ytext_yp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    ytext_yp = ytext_yp+0.5;
    
    text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',0,'fontsize',10);%8
    set(gca,'yticklabel','');
    
    
    
    
    x_lim = xlim;
    y_lim = ylim;
    
    xlabel('Response seqs','Position',[mean(x_lim) max(y_lim)+2.5],'FontSize',12,'FontWeight','bold');
    ylabel('Stimuli seqs','Position',[min(x_lim)-3 mean(y_lim)], 'FontSize', 12, 'FontWeight', 'bold');
    % xlabel('Response seqs', 'FontSize', 12, 'FontWeight', 'bold')
    % ylabel('Stimuli seqs', 'FontSize', 12, 'FontWeight', 'bold');
    
    temp_title = title(sprintf('Behavior data'),'FontSize',11);
    temp_title.Interpreter = 'none';
    
    %     %c = colorbar('FontSize',8);
    %     %c.Layout.Tile = 'east';%east
    %     c = colorbar('eastoutside','FontSize',8);
    %
    %     % c.Ticks = 0:0.2:1;
    %     c.Ticks = (0:0.2:1).^temp_coeff;
    %     c.TickLabels = {'0','0.2','0.4','0.6','0.8','1'};
    
    
    %% behavior model
    %nexttile
    nexttile([1 5])
    
    C = x2;
    %C_max = 1;
    
    %imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
    hold on
    
    
    % Plot length bound in stimuli
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([0.5 size(C,2)+0.5], [temptempi temptempi]-0.5, 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    % Plot length bound in response
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    plot(sum(numSeq).*[1 1]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
    hold on
    
    
    set(gca,'linewidth',1.5)
    %set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8) %12
    
    set(gca,'TickLength',[0 0]);
    set(gca,'XTick',1:size(C,1));
    set(gca,'YTick',1:size(C,2));
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    
    tempLabelStr_res = string(1:5);
    xtl=string(tempLabelStr_res);
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    %xtext_xp=xt;
    xtext_xp = nan(1,4);
    for tempi=1:length(xtext_xp)
        xtext_xp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    xtext_xp = [xtext_xp xtext_xp(end)+20/2];
    xtext_xp = xtext_xp+0.5;
    
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(tempLabelStr_res))+0.75-3.75;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_res))-1.6;%-0.75
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',10);%8
    set(gca,'xticklabel','');
    
    
    
    x_lim = xlim;
    y_lim = ylim;
    
    xlabel('Response seqs','Position',[mean(x_lim) max(y_lim)+2.5],'FontSize',12,'FontWeight','bold');
    % xlabel('Response seqs', 'FontSize', 12, 'FontWeight', 'bold')
    %ylabel('Stimuli seqs', 'FontSize', 12, 'FontWeight', 'bold');
    
    a = 1;
    
    temp_title = title(sprintf('Behavior model'),'FontSize',11);
    temp_title.Interpreter = 'none';
    
    %c = colorbar('FontSize',8);
    %c.Layout.Tile = 'east';%east
    c = colorbar('eastoutside','FontSize',8);
    
    % c.Ticks = 0:0.2:1;
    c.Ticks = (0:0.2:1).^temp_coeff;
    c.TickLabels = {'0','0.2','0.4','0.6','0.8','1'};
    
end


%% Plot B stimuli-to-response matrix of behavior data and behavior model
if false
    fig = figure('Name','responseMatrix','NumberTitle','off'); %#ok<*UNRCH>
    set(gcf,'Position',[10 50 597*0.982*0.95 243*2*0.9*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,10,'TileSpacing','Loose','Padding','Compact');
    
    my_gray = gray;
    my_gray = my_gray(end:-1:1,:);
    colormap(my_gray);
    
    %% behavior data
    %nexttile
    nexttile([1 5])
    
    C = x1;
    
    imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
    hold on
    
    
    % Plot length bound in stimuli
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([0.5 size(C,2)+0.5], [temptempi temptempi]-0.5, 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    % Plot length bound in response
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    plot(sum(numSeq).*[1 1]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
    hold on
    
    
    set(gca,'linewidth',1.5)
    set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8) %12
    
    set(gca,'TickLength',[0 0]);
    set(gca,'XTick',1:size(C,1));
    set(gca,'YTick',1:size(C,2));
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    
    tempLabelStr_res = string(1:5);
    xtl=string(tempLabelStr_res);
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    %xtext_xp=xt;
    xtext_xp = nan(1,4);
    for tempi=1:length(xtext_xp)
        xtext_xp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    xtext_xp = [xtext_xp xtext_xp(end)+20/2];
    xtext_xp = xtext_xp+0.5;
    
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(tempLabelStr_res))+0.75-3.75;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_res))-1.6;%-0.75
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',10);%8
    set(gca,'xticklabel','');
    
    
    
    tempLabelStr_sti = string(1:4);
    ytl=string(tempLabelStr_sti);
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_sti))-2.4;%-0.75
    % 设置ttext的y坐标位置
    % ytext_yp=yt;
    ytext_yp = nan(1,4);
    for tempi=1:length(ytext_yp)
        ytext_yp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    ytext_yp = ytext_yp+0.5;
    
    text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',0,'fontsize',10);%8
    set(gca,'yticklabel','');
    
    
    
    
    x_lim = xlim;
    y_lim = ylim;
    
    xlabel('Response seqs','Position',[mean(x_lim) max(y_lim)+2.5],'FontSize',12,'FontWeight','bold');
    ylabel('Stimuli seqs','Position',[min(x_lim)-3 mean(y_lim)], 'FontSize', 12, 'FontWeight', 'bold');
    % xlabel('Response seqs', 'FontSize', 12, 'FontWeight', 'bold')
    % ylabel('Stimuli seqs', 'FontSize', 12, 'FontWeight', 'bold');
    
    temp_title = title(sprintf('Behavior data'),'FontSize',11);
    temp_title.Interpreter = 'none';
    
    
    %% behavior model
    %nexttile
    nexttile([1 5])
    
    C = x2;
    %C_max = 1;
    
    %imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
    hold on
    
    
    % Plot length bound in stimuli
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([0.5 size(C,2)+0.5], [temptempi temptempi]-0.5, 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    % Plot length bound in response
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    plot(sum(numSeq).*[1 1]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
    hold on
    
    
    set(gca,'linewidth',1.5)
    set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8) %12
    
    set(gca,'TickLength',[0 0]);
    set(gca,'XTick',1:size(C,1));
    set(gca,'YTick',1:size(C,2));
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    
    tempLabelStr_res = string(1:5);
    xtl=string(tempLabelStr_res);
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    %xtext_xp=xt;
    xtext_xp = nan(1,4);
    for tempi=1:length(xtext_xp)
        xtext_xp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    xtext_xp = [xtext_xp xtext_xp(end)+20/2];
    xtext_xp = xtext_xp+0.5;
    
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(tempLabelStr_res))+0.75-3.75;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_res))-1.6;%-0.75
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',10);%8
    set(gca,'xticklabel','');
    
    
    
    x_lim = xlim;
    y_lim = ylim;
    
    xlabel('Response seqs','Position',[mean(x_lim) max(y_lim)+2.5],'FontSize',12,'FontWeight','bold');
    % xlabel('Response seqs', 'FontSize', 12, 'FontWeight', 'bold')
    %ylabel('Stimuli seqs', 'FontSize', 12, 'FontWeight', 'bold');
    
    a = 1;
    
    temp_title = title(sprintf('Behavior model'),'FontSize',11);
    temp_title.Interpreter = 'none';
    
    %c = colorbar('FontSize',8);
    %c.Layout.Tile = 'east';%east
    %c = colorbar('eastoutside','FontSize',8);
    c = colorbar('southoutside','FontSize',8);
    %     c = colorbar('manual','FontSize',8);
    %c.Layout.Tile = 'south';%east
    
    % c.Ticks = 0:0.2:1;
    c.Ticks = (0:0.2:1).^temp_coeff;
    c.TickLabels = {'0','0.2','0.4','0.6','0.8','1'};
    c.Position = [0.1,0.47,0.8,0.025];
end



%% Plot C stimuli-to-response matrix of behavior data and behavior model
if false

    fig = figure('Name','responseMatrix','NumberTitle','off');
    set(gcf,'Position',[10 50 597*0.982 243*2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,10,'TileSpacing','Loose','Padding','Compact');
    
    my_gray = gray;
    my_gray = my_gray(end:-1:1,:);
    colormap(my_gray);
    
    %% behavior data
    %nexttile
    nexttile([1 5])
    
    C = x1;
    
    imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
    hold on
    
    
    % Plot length bound in stimuli
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([0.5 size(C,2)+0.5], [temptempi temptempi]-0.5, 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    % Plot length bound in response
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    plot(sum(numSeq).*[1 1]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
    hold on
    
    % Highlight example seq in stimuli
    plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq-0.5, 'color', [0.75 0.25 0.25]);
    hold on
    plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq+0.5, 'color', [0.75 0.25 0.25]);
    hold on
    
    set(gca,'linewidth',1.5)
    %set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8) %12
    
    set(gca,'TickLength',[0 0]);
    set(gca,'XTick',1:size(C,1));
    set(gca,'YTick',1:size(C,2));
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    
    tempLabelStr_res = string(1:5);
    xtl=string(tempLabelStr_res);
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    %xtext_xp=xt;
    xtext_xp = nan(1,4);
    for tempi=1:length(xtext_xp)
        xtext_xp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    xtext_xp = [xtext_xp xtext_xp(end)+20/2];
    xtext_xp = xtext_xp+0.5;
    
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(tempLabelStr_res))+0.75-3.75;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_res))-1.6;%-0.75
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',10);%8
    set(gca,'xticklabel','');
    
    
    
    tempLabelStr_sti = string(1:4);
    ytl=string(tempLabelStr_sti);
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_sti))-2.4;%-0.75
    % 设置ttext的y坐标位置
    % ytext_yp=yt;
    ytext_yp = nan(1,4);
    for tempi=1:length(ytext_yp)
        ytext_yp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    ytext_yp = ytext_yp+0.5;
    
    text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',0,'fontsize',10);%8
    set(gca,'yticklabel','');
    
    
    
    
    x_lim = xlim;
    y_lim = ylim;
    
    xlabel('Response seqs','Position',[mean(x_lim) max(y_lim)+2.5],'FontSize',12,'FontWeight','bold');
    ylabel('Stimuli seqs','Position',[min(x_lim)-3 mean(y_lim)], 'FontSize', 12, 'FontWeight', 'bold');
    % xlabel('Response seqs', 'FontSize', 12, 'FontWeight', 'bold')
    % ylabel('Stimuli seqs', 'FontSize', 12, 'FontWeight', 'bold');
    
    temp_title = title(sprintf('Behavior data'),'FontSize',11);
    temp_title.Interpreter = 'none';
    
    %     %c = colorbar('FontSize',8);
    %     %c.Layout.Tile = 'east';%east
    %     c = colorbar('eastoutside','FontSize',8);
    %
    %     % c.Ticks = 0:0.2:1;
    %     c.Ticks = (0:0.2:1).^temp_coeff;
    %     c.TickLabels = {'0','0.2','0.4','0.6','0.8','1'};
    
    
    %% behavior model
    %nexttile
    nexttile([1 5])
    
    C = x2;
    %C_max = 1;
    
    %imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
    hold on
    
    
    % Plot length bound in stimuli
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([0.5 size(C,2)+0.5], [temptempi temptempi]-0.5, 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    % Plot length bound in response
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    plot(sum(numSeq).*[1 1]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
    hold on
    
    
    % Highlight example seq in stimuli
    plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq-0.5, 'color', [0.75 0.25 0.25]);
    hold on
    plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq+0.5, 'color', [0.75 0.25 0.25]);
    hold on
    
    
    set(gca,'linewidth',1.5)
    %set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8) %12
    
    set(gca,'TickLength',[0 0]);
    set(gca,'XTick',1:size(C,1));
    set(gca,'YTick',1:size(C,2));
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    
    tempLabelStr_res = string(1:5);
    xtl=string(tempLabelStr_res);
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    %xtext_xp=xt;
    xtext_xp = nan(1,4);
    for tempi=1:length(xtext_xp)
        xtext_xp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    xtext_xp = [xtext_xp xtext_xp(end)+20/2];
    xtext_xp = xtext_xp+0.5;
    
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(tempLabelStr_res))+0.75-3.75;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_res))-1.6;%-0.75
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',10);%8
    set(gca,'xticklabel','');
    
    
    
    x_lim = xlim;
    y_lim = ylim;
    
    xlabel('Response seqs','Position',[mean(x_lim) max(y_lim)+2.5],'FontSize',12,'FontWeight','bold');
    % xlabel('Response seqs', 'FontSize', 12, 'FontWeight', 'bold')
    %ylabel('Stimuli seqs', 'FontSize', 12, 'FontWeight', 'bold');
    
    a = 1;
    
    temp_title = title(sprintf('Behavior model'),'FontSize',11);
    temp_title.Interpreter = 'none';
    
    %c = colorbar('FontSize',8);
    %c.Layout.Tile = 'east';%east
    c = colorbar('eastoutside','FontSize',8);
    
    % c.Ticks = 0:0.2:1;
    c.Ticks = (0:0.2:1).^temp_coeff;
    c.TickLabels = {'0','0.2','0.4','0.6','0.8','1'};
    
end




%% Plot D stimuli-to-response matrix of behavior data and behavior model
% exampleSeq = 17;
%if true
if if_locationDistri_behavior0_model1 == 1
    fig = figure('Name','responseMatrix','NumberTitle','off');
    %set(gcf,'Position',[10 50 597*0.982 243*2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 50 597*0.982*0.83*1.03*1.1*0.97 243*2*0.5*1.1]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 50 597*0.982*0.83*1.03*1.1*0.97 243*2*0.5*1.1*1.13]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 50 534.8*0.95 302.0*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,10,'TileSpacing','Loose','Padding','Loose');
    
    my_gray = gray;
    my_gray = my_gray(end:-1:1,:);
    colormap(my_gray);
    
    %% behavior data
    %nexttile
    nexttile([1 5])
    
    C = x1;
    
    imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
    hold on
    
    
    % Plot length bound in stimuli
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([0.5 size(C,2)+0.5], [temptempi temptempi]-0.5, 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    % Plot length bound in response
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    plot(sum(numSeq).*[1 1]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
    hold on
    
    %     % Highlight example seq in stimuli
    %     plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq-0.5, 'color', [0.75 0.25 0.25]);
    %     hold on
    %     plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq+0.5, 'color', [0.75 0.25 0.25]);
    %     hold on
    
    set(gca,'linewidth',1.5)
    %set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8) %12
    
    set(gca,'TickLength',[0 0]);
    set(gca,'XTick',1:size(C,1));
    set(gca,'YTick',1:size(C,2));
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    
    tempLabelStr_res = string(1:5);
    tempLabelStr_res(1) = "Len1";
    tempLabelStr_res(2) = "Len2";
    tempLabelStr_res(3) = "Len3";
    tempLabelStr_res(4) = "Len4";
    tempLabelStr_res(5) = "Len5";
    
    xtl=string(tempLabelStr_res);
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    %xtext_xp=xt;
    xtext_xp = nan(1,4);
    for tempi=1:length(xtext_xp)
        xtext_xp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    xtext_xp = [xtext_xp xtext_xp(end)+20/2];
    xtext_xp = xtext_xp+0.5;
    
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(tempLabelStr_res))+0.75-3.75;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_res))-1.6;%-0.75
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',6.5);%8,10
    set(gca,'xticklabel','');
    
    
    
    tempLabelStr_sti = string(1:4);
    tempLabelStr_sti(1) = "1";%Len1
    tempLabelStr_sti(2) = "Len2";
    tempLabelStr_sti(3) = "Len3";
    tempLabelStr_sti(4) = "Len4";
    
    ytl=string(tempLabelStr_sti);
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_sti))-3;%-0.75
    % 设置ttext的y坐标位置
    % ytext_yp=yt;
    ytext_yp = nan(1,4);
    for tempi=1:length(ytext_yp)
        ytext_yp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    ytext_yp = ytext_yp+0.5;
    
    text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',90,'fontsize',6.5);%8,10
    set(gca,'yticklabel','');
    
    
    
    
    x_lim = xlim;
    y_lim = ylim;
    
    xlabel('Response seqs','Position',[mean(x_lim) max(y_lim)+3.5],'FontSize',9,'FontWeight','bold');%12
    ylabel('Stimuli seqs','Position',[min(x_lim)-4 mean(y_lim)], 'FontSize', 9, 'FontWeight', 'bold');%12
    % xlabel('Response seqs', 'FontSize', 12, 'FontWeight', 'bold')
    % ylabel('Stimuli seqs', 'FontSize', 12, 'FontWeight', 'bold');
    
    temp_title = title(sprintf('Behavior data'),'FontSize',8);%11
    temp_title.Interpreter = 'none';
    
    %     %c = colorbar('FontSize',8);
    %     %c.Layout.Tile = 'east';%east
    %     c = colorbar('eastoutside','FontSize',8);
    %
    %     % c.Ticks = 0:0.2:1;
    %     c.Ticks = (0:0.2:1).^temp_coeff;
    %     c.TickLabels = {'0','0.2','0.4','0.6','0.8','1'};
    
    c = colorbar('westoutside','FontSize',8);
    c.Position = c.Position+[-0.05 0.00 -0.01 -0.00];
    c.Ticks = [0 1];
    
    %% behavior model
    %nexttile
    nexttile([1 5])
    
    C = x2;
    %C_max = 1;
    
    %imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
    hold on
    
    
    % Plot length bound in stimuli
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([0.5 size(C,2)+0.5], [temptempi temptempi]-0.5, 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    % Plot length bound in response
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    plot(sum(numSeq).*[1 1]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
    hold on
    
    
    %     % Highlight example seq in stimuli
    %     plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq-0.5, 'color', [0.75 0.25 0.25]);
    %     hold on
    %     plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq+0.5, 'color', [0.75 0.25 0.25]);
    %     hold on
    
    
    set(gca,'linewidth',1.5)
    %set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8) %12
    
    set(gca,'TickLength',[0 0]);
    set(gca,'XTick',1:size(C,1));
    set(gca,'YTick',1:size(C,2));
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    
    tempLabelStr_res = string(1:5);
    tempLabelStr_res(1) = "Len1";
    tempLabelStr_res(2) = "Len2";
    tempLabelStr_res(3) = "Len3";
    tempLabelStr_res(4) = "Len4";
    tempLabelStr_res(5) = "Len5";
    
    xtl=string(tempLabelStr_res);
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    %xtext_xp=xt;
    xtext_xp = nan(1,4);
    for tempi=1:length(xtext_xp)
        xtext_xp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    xtext_xp = [xtext_xp xtext_xp(end)+20/2];
    xtext_xp = xtext_xp+0.5;
    
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(tempLabelStr_res))+0.75-3.75;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_res))-1.6;%-0.75
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',6.5);%8
    set(gca,'xticklabel','');
    
    
    
    x_lim = xlim;
    y_lim = ylim;
    
    xlabel('Response seqs','Position',[mean(x_lim) max(y_lim)+3.5],'FontSize',9,'FontWeight','bold');
    % xlabel('Response seqs', 'FontSize', 12, 'FontWeight', 'bold')
    %ylabel('Stimuli seqs', 'FontSize', 12, 'FontWeight', 'bold');
    
    a = 1;
    
    temp_title = title(sprintf('Behavior model'),'FontSize',8);
    temp_title.Interpreter = 'none';
    
    %c = colorbar('FontSize',8);
    %c.Layout.Tile = 'east';%east
    
    %c = colorbar('eastoutside','FontSize',8);
    %c.Position = c.Position+[0.031 0.00 -0.01 -0.00];
    %c.Ticks = [0 1];
    
end


if false
    fig = figure('Name','responseMatrix','NumberTitle','off');
    %set(gcf,'Position',[10 50 597*0.982 243*2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 50 597*0.982 243*2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,10,'TileSpacing','Loose','Padding','Compact');
    
    my_gray = gray;
    my_gray = my_gray(end:-1:1,:);
    colormap(my_gray);
    
    %% behavior data
    %nexttile
    nexttile([1 5])
    
    C = x1;
    
    imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
    hold on
    
    
    % Plot length bound in stimuli
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([0.5 size(C,2)+0.5], [temptempi temptempi]-0.5, 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    % Plot length bound in response
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    plot(sum(numSeq).*[1 1]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
    hold on
    
    % Highlight example seq in stimuli
    plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq-0.5, 'color', [0.75 0.25 0.25]);
    hold on
    plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq+0.5, 'color', [0.75 0.25 0.25]);
    hold on
    
    set(gca,'linewidth',1.5)
    %set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8) %12
    
    set(gca,'TickLength',[0 0]);
    set(gca,'XTick',1:size(C,1));
    set(gca,'YTick',1:size(C,2));
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    
    tempLabelStr_res = string(1:5);
    tempLabelStr_res(1) = "Len1";
    tempLabelStr_res(2) = "Len2";
    tempLabelStr_res(3) = "Len3";
    tempLabelStr_res(4) = "Len4";
    tempLabelStr_res(5) = "Len5";
    
    xtl=string(tempLabelStr_res);
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    %xtext_xp=xt;
    xtext_xp = nan(1,4);
    for tempi=1:length(xtext_xp)
        xtext_xp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    xtext_xp = [xtext_xp xtext_xp(end)+20/2];
    xtext_xp = xtext_xp+0.5;
    
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(tempLabelStr_res))+0.75-3.75;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_res))-1.6;%-0.75
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',10);%8
    set(gca,'xticklabel','');
    
    
    
    tempLabelStr_sti = string(1:4);
    tempLabelStr_sti(1) = "Len1";
    tempLabelStr_sti(2) = "Len2";
    tempLabelStr_sti(3) = "Len3";
    tempLabelStr_sti(4) = "Len4";
    
    ytl=string(tempLabelStr_sti);
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_sti))-3;%-0.75
    % 设置ttext的y坐标位置
    % ytext_yp=yt;
    ytext_yp = nan(1,4);
    for tempi=1:length(ytext_yp)
        ytext_yp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    ytext_yp = ytext_yp+0.5;
    
    text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',90,'fontsize',10);%8
    set(gca,'yticklabel','');
    
    
    
    
    x_lim = xlim;
    y_lim = ylim;
    
    xlabel('Response seqs','Position',[mean(x_lim) max(y_lim)+3.5],'FontSize',12,'FontWeight','bold');
    ylabel('Stimuli seqs','Position',[min(x_lim)-4 mean(y_lim)], 'FontSize', 12, 'FontWeight', 'bold');
    % xlabel('Response seqs', 'FontSize', 12, 'FontWeight', 'bold')
    % ylabel('Stimuli seqs', 'FontSize', 12, 'FontWeight', 'bold');
    
    temp_title = title(sprintf('Behavior data'),'FontSize',11);
    temp_title.Interpreter = 'none';
    
    %     %c = colorbar('FontSize',8);
    %     %c.Layout.Tile = 'east';%east
    %     c = colorbar('eastoutside','FontSize',8);
    %
    %     % c.Ticks = 0:0.2:1;
    %     c.Ticks = (0:0.2:1).^temp_coeff;
    %     c.TickLabels = {'0','0.2','0.4','0.6','0.8','1'};
    
    
    %% behavior model
    %nexttile
    nexttile([1 5])
    
    C = x2;
    %C_max = 1;
    
    %imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
    hold on
    
    
    % Plot length bound in stimuli
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([0.5 size(C,2)+0.5], [temptempi temptempi]-0.5, 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    % Plot length bound in response
    for temptempi=1:size(C,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    plot(sum(numSeq).*[1 1]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
    hold on
    
    
    % Highlight example seq in stimuli
    plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq-0.5, 'color', [0.75 0.25 0.25]);
    hold on
    plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq+0.5, 'color', [0.75 0.25 0.25]);
    hold on
    
    
    set(gca,'linewidth',1.5)
    %set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8) %12
    
    set(gca,'TickLength',[0 0]);
    set(gca,'XTick',1:size(C,1));
    set(gca,'YTick',1:size(C,2));
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    
    tempLabelStr_res = string(1:5);
    tempLabelStr_res(1) = "Len1";
    tempLabelStr_res(2) = "Len2";
    tempLabelStr_res(3) = "Len3";
    tempLabelStr_res(4) = "Len4";
    tempLabelStr_res(5) = "Len5";
    
    xtl=string(tempLabelStr_res);
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    %xtext_xp=xt;
    xtext_xp = nan(1,4);
    for tempi=1:length(xtext_xp)
        xtext_xp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    xtext_xp = [xtext_xp xtext_xp(end)+20/2];
    xtext_xp = xtext_xp+0.5;
    
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(tempLabelStr_res))+0.75-3.75;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_res))-1.6;%-0.75
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',10);%8
    set(gca,'xticklabel','');
    
    
    
    x_lim = xlim;
    y_lim = ylim;
    
    xlabel('Response seqs','Position',[mean(x_lim) max(y_lim)+3.5],'FontSize',12,'FontWeight','bold');
    % xlabel('Response seqs', 'FontSize', 12, 'FontWeight', 'bold')
    %ylabel('Stimuli seqs', 'FontSize', 12, 'FontWeight', 'bold');
    
    a = 1;
    
    temp_title = title(sprintf('Behavior model'),'FontSize',11);
    temp_title.Interpreter = 'none';
    
    %c = colorbar('FontSize',8);
    %c.Layout.Tile = 'east';%east
    c = colorbar('eastoutside','FontSize',8);
    
    % c.Ticks = 0:0.2:1;
    
    
    %c.Ticks = (0:0.2:1).^temp_coeff;
    %c.TickLabels = {'0','0.2','0.4','0.6','0.8','1'};
    c.Ticks = [0 1];
    
end



if true
    %% Plot behavior matrix of stimuli-to-response
    fig = figure('Name','responseMatrix','NumberTitle','off');
    %set(gcf,'Position',[10 50 597*0.982 243*2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 50 597*0.982 243*2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 50 597*0.982*0.68*1.1*0.95*0.96*0.96 243*2*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 50 597*0.982*0.68*1.1*0.95*0.96*0.96 243*2*0.8*0.85*0.94*0.97*0.985*1.02]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 50 200*1.15 150*2*0.95*1.03*1.03]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 50 200*1.15*1.56*1.04*0.99 150*2*0.95*1.03*1.03*2.01*1.02*1.07]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 50 200*1.15*1.56*1.04*0.99*1.13*0.995 150*2*0.95*1.03*1.03*2.01*1.02*1.07]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,6,'TileSpacing','Loose','Padding','Loose');
    
    %     my_gray = gray;
    %     my_gray = my_gray(end:-1:1,:);
    %     colormap(my_gray);
    
    %     colormap('gray');
    
    colormap(coolwarm());
    %       colormap('jet');
    
    %         colormap parula
    %         colormap(cividis());
    %         colormap(viridis());
    
    %         temp1 = pink;
    %         temp1 = temp1(end:-1:1,:);
    %         colormap(temp1);
    
    %         colormap('pink');
    
    %         colormap hot
    
    %         temp1 = hot;
    %         temp1 = temp1(end:-1:1,:);
    %         colormap(temp1);
    
    %         colormap cool
    
    %         colormap bone
    
    %         temp1 = bone;
    %         temp1 = temp1(end:-1:1,:);
    %         colormap(temp1);
    
    
    %nexttile
    %set(gca, 'visible', 'off')
    
    %% behavior data
    %nexttile
    nexttile([1 5])
    
    C = x1;
    
    imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
    hold on
    
    
    %     % Plot length bound in stimuli
    %     for temptempi=1:size(C,1)
    %         if temptempi > 1
    %             if length(target_seqSet_inOne{temptempi}) > ...
    %                     length(target_seqSet_inOne{temptempi-1})
    %                 plot([0.5 size(C,2)+0.5], [temptempi temptempi]-0.5, 'color', [0.25 0.25 0.25]);
    %                 hold on
    %             end
    %         end
    %     end
    
    %     % Plot length bound in response
    %     for temptempi=1:size(C,1)
    %         if temptempi > 1
    %             if length(target_seqSet_inOne{temptempi}) > ...
    %                     length(target_seqSet_inOne{temptempi-1})
    %                 plot([temptempi temptempi]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
    %                 hold on
    %             end
    %         end
    %     end
    %     plot(sum(numSeq).*[1 1]-0.5, [0.5 size(C,1)+0.5], 'color', [0.25 0.25 0.25]);
    %     hold on
    
    
    temptempColor = [1 1 1]*0.75;%[0.75 0.25 0.25]
    %     temptempColor = [1 1 1]*0.25;%[0.75 0.25 0.25]
    
    % Highlight example seq in stimuli
    plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq-0.5,':', 'color', temptempColor,'linewidth',0.75);
    hold on
    plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq+0.5,':', 'color', temptempColor,'linewidth',0.75);
    hold on
    
    
    exampleSeq2 = 37;
    
    % Highlight example seq in stimuli
    plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq2-0.5,':', 'color', temptempColor,'linewidth',0.75);
    hold on
    plot([0.5 size(C,2)+0.5], [1 1]*exampleSeq2+0.5,':', 'color', temptempColor,'linewidth',0.75);
    hold on
    
    set(gca,'linewidth',1.5)
    %set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8) %12
    
    set(gca,'TickLength',[0 0]);
    set(gca,'XTick',1:size(C,2));
    set(gca,'YTick',1:size(C,1));
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    
    tempLabelStr_res = string(1:5);
    tempLabelStr_res(1) = "Length1";
    tempLabelStr_res(2) = "Length2";
    tempLabelStr_res(3) = "Length3";
    tempLabelStr_res(4) = "Length4";
    tempLabelStr_res(5) = "Length5";
    
    xtl=string(tempLabelStr_res);
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    %xtext_xp=xt;
    xtext_xp = nan(1,4);
    for tempi=1:length(xtext_xp)
        xtext_xp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    xtext_xp = [xtext_xp xtext_xp(end)+20/2];
    xtext_xp = xtext_xp+0.5;
    
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(tempLabelStr_res))+0.75-3.75+6;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_res))-1.6;%-0.75
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',6.5);%8
    set(gca,'xticklabel','');
    
    
%         xtl=string(seqSet_inOne_inOne(1:size(C,1)));%2
%         xt=get(gca,'XTick');
%         yt=get(gca,'YTick');
%         % 设置xtext的x坐标位置
%         xtext_xp=xt;
%         % 设置xtext的y坐标位置
%         xtext_yp=(yt(end))*ones(1,length(xt))+0.75+0.25;
%         % 设置ttext的x坐标位置
%         ytext_xp=(xt(1))*ones(1,length(yt))-1.6;%-0.75
%         % 设置ttext的y坐标位置
%         ytext_yp=yt;
%         text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',90,'fontsize',6);%8
%         set(gca,'xticklabel','');
%         xtickangle(0);

    
    tempLabelStr_sti = string(1:4);
    tempLabelStr_sti(1) = "Length1";
    tempLabelStr_sti(2) = "Length2";
    tempLabelStr_sti(3) = "Length3";
    tempLabelStr_sti(4) = "Length4";
    
    ytl=string(tempLabelStr_sti);
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(tempLabelStr_sti))-3;%-0.75
    % 设置ttext的y坐标位置
    % ytext_yp=yt;
    ytext_yp = nan(1,4);
    for tempi=1:length(ytext_yp)
        ytext_yp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    end
    ytext_yp = ytext_yp+0.5;
    
    text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',90,'fontsize',6.5);%8
    set(gca,'yticklabel','');
    
% 
%         ytl=string(seqSet_inOne_inOne(1:size(C,1)));%2
%         xt=get(gca,'XTick');
%         yt=get(gca,'YTick');
%         % 设置xtext的x坐标位置
%         ytext_xp=yt;
%         % 设置xtext的y坐标位置
%         ytext_yp=(yt(end))*ones(1,length(xt))+0.75-5.5;
%         % 设置ttext的x坐标位置
%         ytext_xp=(xt(1))*ones(1,length(yt))-1;%-0.75
%         % 设置ttext的y坐标位置
%         ytext_yp=yt;
%         text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',6);%8
%         set(gca,'yticklabel','');
%         ytickangle(0);

    
    
    
    x_lim = xlim;
    y_lim = ylim;
    
    xlabel('Response seqs','Position',[mean(x_lim) max(y_lim)+3.5],'FontSize',9,'FontWeight','normal');
    ylabel('Stimuli seqs','Position',[min(x_lim)-4 mean(y_lim)], 'FontSize', 9, 'FontWeight', 'normal');
    % xlabel('Response seqs', 'FontSize', 12, 'FontWeight', 'bold')
    % ylabel('Stimuli seqs', 'FontSize', 12, 'FontWeight', 'bold');
    
    %temp_title = title(sprintf('Stimuli-to-response matrix'),'FontSize',9);
    %temp_title.Interpreter = 'none';
    
    %c = colorbar('FontSize',8);
    %c.Layout.Tile = 'east';%east
    c = colorbar('eastoutside','FontSize',6.5);
    %c = colorbar('westoutside','FontSize',6.5);
    
    % c.Ticks = 0:0.2:1;
    
    %c.Position = c.Position+[-0.08 0.00 -0.04 -0.00];
    %c.Position = c.Position+[-0.08 0.00 -0.02 -0.2];
    c.Position = c.Position+[+0.05 0.218 -0.018 -0.22];
    
    %c.Ticks = (0:0.2:1).^temp_coeff;
    %c.TickLabels = {'0','0.2','0.4','0.6','0.8','1'};
    c.Ticks = [0 1];
    
end


%% Testing for behavior accuracy ~ memory precision + length + memory precision * length
if true
    seqPrecision_behavior;
    seqPrecision_neuron;
    
    temp_if_zscore = 1;
    
    temp_seqLength = sum(boolIndex_location_seq_T,2);
    temp_seqLength = temp_seqLength(1:41);
    
    %x = [seqPrecision_neuron,temp_seqLength];
    %x = [seqPrecision_neuron,1./temp_seqLength];
    %x = [seqPrecision_neuron,-1.*temp_seqLength];
    %x = [seqPrecision_neuron,6-temp_seqLength];
    %x = [1./seqPrecision_neuron,temp_seqLength];
    %x = [-1.*seqPrecision_neuron,temp_seqLength];
    
    x1 = seqPrecision_neuron;
    x2 = 1./temp_seqLength;
    
    x1_n11n = (x1-mean(x1,'omitnan'))/std(x1,'omitnan'); % z-score
    x2_n11n = (x2-mean(x2,'omitnan'))/std(x2,'omitnan'); % z-score
    
    if temp_if_zscore == 0
        x = [x1,x2];
    elseif temp_if_zscore == 1
        x = [x1_n11n,x2_n11n];
    end
    
    y = seqPrecision_behavior;
    
    if temp_if_zscore == 1
        y = (y-mean(y,'omitnan'))/std(y,'omitnan'); % z-score
    end
    
    temp_mdl= fitglm(x,y,'interactions','Distribution','normal','Intercept',true);
    %temp_mdl= fitglm(x,y,'linear','Distribution','normal','Intercept',true)
    temp_glm_beta = temp_mdl.Coefficients.Estimate;
    temp_glm_r2 = temp_mdl.Rsquared.Adjusted;
    
    
end


%% Plot example seq probability distribution
if true == 1
    
    %exampleSeq;
    responseMatrix_noChoice_n11n;
    
    temp_exampleSeq = 8;%8,37
    temp_seqProbDistri = responseMatrix_noChoice_n11n(temp_exampleSeq,:);
    
    fig = figure('Name','Sequence decoding','NumberTitle','off'); %#ok<*NASGU>
    set(gcf,'Position',[50+0 50+0 240*1.15*1.15 252*0.925*0.4*0.9*0.95*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    x = 1:length(temp_seqProbDistri);
    y =  temp_seqProbDistri;
    

    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y);
    
    temp_ymin = 0;
    temp_ymax = 1;
    
    
    for tempi=1:length(x)
        if tempi == temp_exampleSeq
        plot([1 1].*x(tempi),[temp_ymin y(tempi)],'linewidth',1,'color',[206 97 26]./255);
        else
        plot([1 1].*x(tempi),[temp_ymin y(tempi)],'linewidth',1,'color',[0 0 0]);
        end
        hold on
    end
    
    %plot(x,y,'--','linewidth',1,'color',[0 0 0]);
    %hold on
    
    temp_size = 15;%10
    scatter(x,y,temp_size,[0 0 0],'filled');
    hold on  
    
    scatter(x(temp_exampleSeq),y(temp_exampleSeq),temp_size,[206 97 26]./255,'filled');
    hold on      
        
    
    %xlim([temp_xmin-(temp_xmax-temp_xmin)*0.1 temp_xmax+(temp_xmax-temp_xmin)*0.1]);
    ylim([temp_ymin temp_ymax]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    
    xticks(temp_exampleSeq);
    if temp_exampleSeq == 8
        set(gca,'xticklabel',[ '\color[rgb]',sprintf('{%s}13',num2str([206 97 26]./255)) ]);
    elseif temp_exampleSeq == 37
        set(gca,'xticklabel',[ '\color[rgb]',sprintf('{%s}256',num2str([206 97 26]./255)) ]);
    end
    
    set(gca,'YTick',0:1,'FontSize', 8);%给坐标加标签
    ax = gca;
    ax.XAxis.TickLength = [0 0];

    xlabel(sprintf('Response seqs'), 'FontSize', 8);
    ylabel(sprintf('Probability'), 'FontSize', 8);
    %temp_title = title(sprintf('%s\n%d seqs, r=%.3f, p=%.3f',currentSession,length(y),r,p), 'FontSize', 12, 'FontWeight', 'bold');
    %temp_title.Interpreter = 'none';
    
    
end


%% End