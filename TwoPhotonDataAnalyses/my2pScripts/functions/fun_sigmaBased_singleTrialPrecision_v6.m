function temp_precision = fun_sigmaBased_singleTrialPrecision_v6(options)

if_fit_gain = 0;


% tempSeqIndex = options.tempSeqIndex;
boolIndex_location_seq_T = options.boolIndex_location_seq_T;
% Posterior_2d_n11n_mean = options.Posterior_2d_n11n_mean;
numSeq = options.numSeq;
valid_length = options.valid_length;
% score_stimuli_to_response = options.score_stimuli_to_response;
if_precision_meanProb0_sumProb1 = options.if_precision_meanProb0_sumProb1;

temp_locDistri = options.temp_locDistri;
temp_score = options.temp_score;
if_entropy = options.if_entropy;
tempSeqIndex = options.tempSeqIndex;
tempSeqLength = options.tempSeqLength;

% temp_trialNum = size(temp_locDistri,1);

% tempSeqIndex = 10;
% boolIndex_location_seq_T;
% Posterior_2d_n11n_mean;
% numSeq;
% valid_length;
% score_stimuli_to_response;
% if_precision_meanProb0_sumProb1;



%% temp_seqProb_stimuli_to_response_n11n
% temp_locDistri = Posterior_2d_n11n_mean(tempSeqIndex,:);
temp_locDistri; %#ok<*VUNUS>

temp_seqProb_stimuli_to_response = nan(1,sum(numSeq(1:valid_length)));
for tempj=1:sum(numSeq(1:valid_length))
    temp_boolIndex_location_seq = boolIndex_location_seq_T(tempj,:);
    
    temp_p = temp_locDistri;
    temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
    temp_p_seq = prod(temp_p,2);
    
    temp_seqProb_stimuli_to_response(tempj) = temp_p_seq;
end
temp_seqProb_stimuli_to_response_n11n = temp_seqProb_stimuli_to_response./sum(temp_seqProb_stimuli_to_response);




%% single trial precision
if if_entropy == 0
    % temp_score = score_stimuli_to_response(tempSeqIndex,1:sum(numSeq(1:valid_length)));
    temp_score;
    temp_seqProb = temp_seqProb_stimuli_to_response_n11n;
    
    if if_precision_meanProb0_sumProb1 == 0
        x = temp_score';
        y = temp_seqProb';
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
    
    temp_r2 = gof.rsquare; %#ok<*NASGU>
    temp_sigma = fitresult.sigma;
    
    temp_precision = abs(1/temp_sigma);
    % temp_precision = 1/temp_sigma;
    
    
elseif if_entropy == 1
    a = 1;
        
    if_n11n = 1;%1
    
    if if_n11n == 0
        temp_seqProb = temp_seqProb_stimuli_to_response;
    elseif if_n11n == 1
        temp_seqProb = temp_seqProb_stimuli_to_response_n11n;
    end
    
    
    
%     % Assume uniform distribution across all length seqs
    p = sum(temp_seqProb).*ones(1,sum(numSeq(1:valid_length)))./sum(numSeq(1:valid_length));
    
    
%     % Assume uniform distribution across current length seqs
%     p = zeros(1,sum(numSeq(1:valid_length)));    
%     tempSeqLength;
%     temp_range = sum(numSeq(1:(tempSeqLength-1)))+1:sum(numSeq(1:tempSeqLength));
%     p(temp_range) = 1/length(temp_range);

    
    
%     boolIndex_location_seq_T;
%     tempSeqIndex;    
%     temp_boolIndex_location_seq = boolIndex_location_seq_T(tempSeqIndex,:);    
%     temp1 = temp_boolIndex_location_seq & boolIndex_location_seq_T(1:sum(numSeq(1:valid_length)),:);
%     temptempBoolIndex = sum(temp1,2) > 0;
%     temp_range = sum(numSeq(1:(tempSeqLength-1)))+1:sum(numSeq(1:tempSeqLength));
%     temptempBoolIndex(temp_range) = true;
    
    
    
    
%     %Assume uniform distribution across all locations
%     temp_locDistri_random = tempSeqLength .* ones(1,6)./6;
%     p = nan(1,sum(numSeq(1:valid_length)));
%     for tempj=1:sum(numSeq(1:valid_length))
%         temp_boolIndex_location_seq = boolIndex_location_seq_T(tempj,:);
%         
%         temp_p = temp_locDistri_random;
%         temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
%         temp_p_seq = prod(temp_p,2);
%         
%         p(tempj) = temp_p_seq;
%     end
%     p = p./sum(p);
    
    
    p = p + eps;
    if if_n11n == 1
        p = p./sum(p);
    end
    entropy_max = -sum(p.*log2(p));
    
    
    
    
    p = temp_seqProb;
    
%     %Assume no length error for decoder
%     temp_range = sum(numSeq(1:(tempSeqLength-1)))+1:sum(numSeq(1:tempSeqLength));
%     temptempBoolIndex = false(1,sum(numSeq(1:valid_length))); 
%     temptempBoolIndex(temp_range) = true;    
%     p = p.*temptempBoolIndex;
    
    
    p = p + eps;
    if if_n11n == 1
        p = p./sum(p);
    end
    entropy = -sum(p.*log2(p));
    
    %temp_precision = 1./entropy;
    %temp_precision = entropy_max-entropy;
    temp_precision = (entropy_max-entropy)./entropy_max;
    
    a = 1;
end





%% End