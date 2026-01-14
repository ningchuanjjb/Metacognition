function [r2Valid_boolIndex_inOne]=...
    fun_r2Valid_v4(if_glm_12B0_24B1_6B2,r2_threshold,if_singleFOV,B_fileName) %#ok<*INUSL>
%% Initialization

% if_singleFOV = 0;

%%
targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)

% %% get B_fileName
% if if_glm_12B0_24B1_6B2 == 1
%     if if_singleFOV == 1
%         B_fileName = 'glm_beta_singleFOV_allTrial_allLength_full_weighted_24B.mat';
%     elseif if_singleFOV == 0
%         B_fileName = 'glm_beta_multiFOV_allTrial_allLength_full_weighted_24B.mat';
%     end
% elseif if_glm_12B0_24B1_6B2 == 0
%     if if_singleFOV == 1
%         B_fileName = 'glm_beta_singleFOV_allTrial_allLength_full_weighted_12plus12B.mat';
%     elseif if_singleFOV == 0
%         B_fileName = 'glm_beta_multiFOV_allTrial_allLength_full_weighted_12plus12B.mat';
%     end
% elseif if_glm_12B0_24B1_6B2 == 2
%     if if_singleFOV == 1
%         B_fileName = 'glm_beta_singleFOV_allTrial_allLength_full_weighted_6B.mat';
%     elseif if_singleFOV == 0
%         B_fileName = 'glm_beta_multiFOV_allTrial_allLength_full_weighted_6B.mat';
%     end
% end
% fprintf('B_fileName=%s.\n',B_fileName);


%% load
load(['C:\ASDROOT\STUDY\temp','\data\',B_fileName]); %#ok<*LOAD>

glm_beta_lengthx_baselineBin_multiFOV;
glm_beta_lengthx_delay1Bin_multiFOV;
glm_beta_lengthx_delay2Bin_multiFOV;

a = 1;
[num_roi,num_B] = size(glm_beta_lengthx_baselineBin_multiFOV); %#ok<*ASGLU>

%% Compute r2Valid_boolIndex
r2Valid_boolIndex = true(num_roi,3);
for tempIndex=1:3    
    if exist('glm_r2_lengthx_delay1Bin_multiFOV','var') == 1
        %r2_threshold = 0.00;%0.05-->0.1
        if tempIndex == 1
            glm_r2 = glm_r2_lengthx_baselineBin_multiFOV;
        elseif tempIndex == 2
            glm_r2 = glm_r2_lengthx_delay1Bin_multiFOV;
        elseif tempIndex == 3
            glm_r2 = glm_r2_lengthx_delay2Bin_multiFOV;
        end
        temp_r2Valid_boolIndexMulti = glm_r2>r2_threshold;
        %temp_r2Valid_boolIndex = sum(temp_r2Valid_boolIndexMulti,2)>=4;
        temp_r2Valid_boolIndex = temp_r2Valid_boolIndexMulti(:,1);
        
        r2Valid_boolIndex(:,tempIndex) = temp_r2Valid_boolIndex;
    end
end
r2Valid_boolIndex_inOne = r2Valid_boolIndex(:,2) | r2Valid_boolIndex(:,3);
fprintf('r2 valid roi (threshold %.3f) = %d/%d [%d,%d,%d].\n',r2_threshold,sum(r2Valid_boolIndex_inOne),...
    size(r2Valid_boolIndex_inOne,1),sum(r2Valid_boolIndex(:,1)),sum(r2Valid_boolIndex(:,2)),sum(r2Valid_boolIndex(:,3)));

r2Valid_boolIndex_inOne; %#ok<*VUNUS>
a = 1; %#ok<*NASGU>


%% End