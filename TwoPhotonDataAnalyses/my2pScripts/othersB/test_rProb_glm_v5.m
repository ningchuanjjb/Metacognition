%% Initialization
close all

offloadingProb_inOne;
numSeq;
target_seqSet;


F_dff_delay1Bin;
F_dff_delay2Bin;
F_dff_baselineBin;
trial_para;



% threshold_r2 = 0;
% threshold_p = 0.05;

temp_offloadingProb_inOne = offloadingProb_inOne;

temp_F_dff_decisionBin = F_dff_delay1Bin;
% temp_F_dff_decisionBin = F_dff_delay2Bin;
% temp_F_dff_decisionBin = F_dff_baselineBin;


seqIndex = zeros(1,trial_para.trial_count);
for tempi=1:trial_para.trial_count
    currentSequence = trial_para.currentSequence{tempi};
    temp_seq_length = length(currentSequence);
    for tempj=1:numSeq(temp_seq_length)
        if sum(ismember(currentSequence,target_seqSet{temp_seq_length}{tempj})) == temp_seq_length
            break
        end
    end
    seqIndex(tempi) = sum(numSeq(1:temp_seq_length-1)) + tempj;
end
offloadingProb_inSeq = offloadingProb_inOne;
roiNum = size(temp_F_dff_decisionBin,1);


temp_F_dff_decisionBin_seqMerged = zeros(roiNum,sum(numSeq));
for tempi=1:sum(numSeq)
    temp_dff = temp_F_dff_decisionBin(:,seqIndex==tempi);
    temp_F_dff_decisionBin_seqMerged(:,tempi) = mean(temp_dff,2);
end


temp_F_dff_decisionBin = temp_F_dff_decisionBin_seqMerged;
temp_offloadingProb = offloadingProb_inSeq;

%% figglm
x = temp_F_dff_decisionBin';
y = temp_offloadingProb';

mdl = cell(roiNum,1);
beta = zeros(roiNum,2);
r2 = zeros(roiNum,1);
p = zeros(roiNum,1);
t0 = tic;
warning('off');
% for tempi=1:roiNum
parfor tempi=1:roiNum    
    temp_mdl = fitglm(x(:,tempi),y);        
    beta(tempi,:) = temp_mdl.Coefficients.Estimate; %#ok<*NASGU>
    r2(tempi) = temp_mdl.Rsquared.Adjusted;
    p(tempi) = temp_mdl.Coefficients.pValue(2);
    mdl{tempi} = temp_mdl;
end
warning('on');
fprintf('t=%.1f secs.\n',toc(t0));

% [r2_sorted,I_sorted] = sort(r2,'descend');

a = 1;

r2;
p;

% temp1 = r2>0;
% temp1B = sum(temp1);
% temp2 = p<0.05;
% % temp2 = p<0.001;
% temp2B = sum(temp2)
% sum(temp1&temp2)
% 
% temptemp_r2 = r2(temp2);
% min(temptemp_r2)

a = 1;



%% End