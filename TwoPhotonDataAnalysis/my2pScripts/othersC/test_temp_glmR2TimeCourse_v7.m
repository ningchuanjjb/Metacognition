%% Initialization
memoryPrecision_trialLevel;
meta_trialLevel;

F_dff_baselinePeriod;
F_dff_length1_sample;
F_dff_length2_sample;
F_dff_length3_sample;
F_dff_decisionPeriodA;

if_compute = 0;

if_plot = 1;
if_plot_multiPanel = 1;

if_shuffle = 0;
shuffleNum = 100;%10

shufflePrctileThreshold = 95;%99

topX = 300;


if_p_shuffle0_baseline1 = 1;


%% Compute
if if_compute == 1
    for bigLoopCount=1:2
        if bigLoopCount == 1
            trialLabel = memoryPrecision_trialLevel;
        elseif bigLoopCount == 2
            trialLabel = meta_trialLevel;
        end
        
        F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
        temp_r2_roi = nan(1,roiNum);
        parfor tempi=1:roiNum
            x = F_dff_decisionBin1(tempi,:)';
            y = trialLabel;
            temp_mdl = fitglm(x,y,'linear');
            temp_r2_roi(tempi) = temp_mdl.Rsquared.Adjusted;
        end
        [M,I] = sort(temp_r2_roi,'descend');
        I_topX = I(1:topX);
        
        tempCellBoolIndex = false(1,roiNum);
        tempCellBoolIndex(I_topX) = true;
        
        for loopCount=0:4
            
            if loopCount==0
                F_dff_current = F_dff_baselinePeriod;
            elseif loopCount==1
                F_dff_current = F_dff_length1_sample;
            elseif loopCount==2
                F_dff_current = F_dff_length2_sample;
            elseif loopCount==3
                F_dff_current = F_dff_length3_sample;
            elseif loopCount==4
                F_dff_current = F_dff_decisionPeriodA;
            end
            
            temptemp_range = 1:size(F_dff_current,3);
            temp_r2 = nan(1,size(F_dff_current,3));
            parfor tempi=temptemp_range
                x = F_dff_current(tempCellBoolIndex,:,tempi)';
                y = trialLabel;
                
                tempBoolIndex = (~isnan(x(:,1))) & (~isnan(y));
                
                x2 = x(tempBoolIndex,:);
                y2 = y(tempBoolIndex);
                
                temp_mdl = fitglm(x2,y2,'linear');
                temp_r2(tempi) = temp_mdl.Rsquared.Adjusted;
            end
            
            if bigLoopCount == 1
                if loopCount==0
                    memoryPrecision_r2_baselinePeriod = temp_r2;
                elseif loopCount==1
                    memoryPrecision_r2_length1 = temp_r2;
                elseif loopCount==2
                    memoryPrecision_r2_length2 = temp_r2;
                elseif loopCount==3
                    memoryPrecision_r2_length3 = temp_r2;
                elseif loopCount==4
                    memoryPrecision_r2_decisionPeriodA = temp_r2;
                end
                
            elseif bigLoopCount == 2
                if loopCount==0
                    meta_r2_baselinePeriod = temp_r2;
                elseif loopCount==1
                    meta_r2_length1 = temp_r2;
                elseif loopCount==2
                    meta_r2_length2 = temp_r2;
                elseif loopCount==3
                    meta_r2_length3 = temp_r2;
                elseif loopCount==4
                    meta_r2_decisionPeriodA = temp_r2;
                end
                
            end
            
        end
        
    end
    
    
    memoryPrecision_r2_all = {memoryPrecision_r2_baselinePeriod,memoryPrecision_r2_length1,...
        memoryPrecision_r2_length2,memoryPrecision_r2_length3,memoryPrecision_r2_decisionPeriodA};
    
    meta_r2_all = {meta_r2_baselinePeriod,meta_r2_length1,...
        meta_r2_length2,meta_r2_length3,meta_r2_decisionPeriodA};
    
    
    
    
    if if_shuffle == 1
        shuffleNum;
        
        memoryPrecision_r2_all_shuffled = cell(shuffleNum,5);
        meta_r2_all_shuffled = cell(shuffleNum,5);
        
        for temp_shuffleCount=1:shuffleNum
            for bigLoopCount=1:2
                if bigLoopCount == 1
                    trialLabel = memoryPrecision_trialLevel;
                elseif bigLoopCount == 2
                    trialLabel = meta_trialLevel;
                end
                
                temp_trialNum = length(trialLabel);
                tempIndex_shuffle = randperm(temp_trialNum,temp_trialNum);
                trialLabel_shuffle = trialLabel(tempIndex_shuffle);
                
                F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
                temp_r2_roi = nan(1,roiNum);
                parfor tempi=1:roiNum
                    x = F_dff_decisionBin1(tempi,:)';
                    y = trialLabel_shuffle;
                    temp_mdl = fitglm(x,y,'linear');
                    temp_r2_roi(tempi) = temp_mdl.Rsquared.Adjusted;
                end
                [M,I] = sort(temp_r2_roi,'descend');
                I_topX = I(1:topX);
                
                tempCellBoolIndex = false(1,roiNum);
                tempCellBoolIndex(I_topX) = true;
                
                for loopCount=0:4
                    
                    if loopCount==0
                        F_dff_current = F_dff_baselinePeriod;
                    elseif loopCount==1
                        F_dff_current = F_dff_length1_sample;
                    elseif loopCount==2
                        F_dff_current = F_dff_length2_sample;
                    elseif loopCount==3
                        F_dff_current = F_dff_length3_sample;
                    elseif loopCount==4
                        F_dff_current = F_dff_decisionPeriodA;
                    end
                    
                    temptemp_range = 1:size(F_dff_current,3);
                    temp_r2 = nan(1,size(F_dff_current,3));
                    parfor tempi=temptemp_range
                        x = F_dff_current(tempCellBoolIndex,:,tempi)';
                        y = trialLabel_shuffle;
                        
                        tempBoolIndex = (~isnan(x(:,1))) & (~isnan(y));
                        
                        x2 = x(tempBoolIndex,:);
                        y2 = y(tempBoolIndex);
                        
                        temp_mdl = fitglm(x2,y2,'linear');
                        temp_r2(tempi) = temp_mdl.Rsquared.Adjusted;
                    end
                    
                    if bigLoopCount == 1
                        if loopCount==0
                            memoryPrecision_r2_baselinePeriod = temp_r2;
                        elseif loopCount==1
                            memoryPrecision_r2_length1 = temp_r2;
                        elseif loopCount==2
                            memoryPrecision_r2_length2 = temp_r2;
                        elseif loopCount==3
                            memoryPrecision_r2_length3 = temp_r2;
                        elseif loopCount==4
                            memoryPrecision_r2_decisionPeriodA = temp_r2;
                        end
                        
                    elseif bigLoopCount == 2
                        if loopCount==0
                            meta_r2_baselinePeriod = temp_r2;
                        elseif loopCount==1
                            meta_r2_length1 = temp_r2;
                        elseif loopCount==2
                            meta_r2_length2 = temp_r2;
                        elseif loopCount==3
                            meta_r2_length3 = temp_r2;
                        elseif loopCount==4
                            meta_r2_decisionPeriodA = temp_r2;
                        end
                        
                    end
                    
                end
                
            end
            
            
            memoryPrecision_r2_all_shuffled(temp_shuffleCount,:) = {memoryPrecision_r2_baselinePeriod,memoryPrecision_r2_length1,...
                memoryPrecision_r2_length2,memoryPrecision_r2_length3,memoryPrecision_r2_decisionPeriodA};
            
            meta_r2_all_shuffled(temp_shuffleCount,:) = {meta_r2_baselinePeriod,meta_r2_length1,...
                meta_r2_length2,meta_r2_length3,meta_r2_decisionPeriodA};
            
        end
        
        
        memoryPrecision_r2_all_shuffled;
        meta_r2_all_shuffled;
        
        memoryPrecision_r2_all_shuffledB = memoryPrecision_r2_all_shuffled;
        meta_r2_all_shuffledB = meta_r2_all_shuffled;
        
        
    else
        
        memoryPrecision_r2_all_shuffledB = cell(shuffleNum,5);
        meta_r2_all_shuffledB = cell(shuffleNum,5);
        
        for temp_shuffleCount=1:shuffleNum
            memoryPrecision_r2_all_shuffledB(temp_shuffleCount,:) = {memoryPrecision_r2_baselinePeriod.*0,memoryPrecision_r2_length1.*0,...
                memoryPrecision_r2_length2.*0,memoryPrecision_r2_length3.*0,memoryPrecision_r2_decisionPeriodA.*0};
            
            meta_r2_all_shuffledB(temp_shuffleCount,:) = {meta_r2_baselinePeriod.*0,meta_r2_length1.*0,...
                meta_r2_length2.*0,meta_r2_length3.*0,meta_r2_decisionPeriodA.*0};
        end
        
    end
    
end

memoryPrecision_r2_all_shuffled_threshold = memoryPrecision_r2_all;
for temptempi=1:5
    temp1 = memoryPrecision_r2_all_shuffledB(1:end,temptempi);
    temp2 = [];
    for temptempj=1:shuffleNum
        temp2 = [temp2;temp1{temptempj}];
    end
    temp_threshold = prctile(temp2,shufflePrctileThreshold,1);
    memoryPrecision_r2_all_shuffled_threshold{temptempi} = temp_threshold;
end

meta_r2_all_shuffled_threshold = meta_r2_all;
for temptempi=1:5
    temp1 = meta_r2_all_shuffledB(1:end,temptempi);
    temp2 = [];
    for temptempj=1:shuffleNum
        temp2 = [temp2;temp1{temptempj}];
    end
    temp_threshold = prctile(temp2,shufflePrctileThreshold,1);
    meta_r2_all_shuffled_threshold{temptempi} = temp_threshold;
end



memoryPrecision_r2_all_shuffled_threshold;
meta_r2_all_shuffled_threshold;


%%


memoryPrecision_r2_baselinePeriod = memoryPrecision_r2_all{1};
meta_r2_baselinePeriod = meta_r2_all{1};

memoryPrecision_r2_baselinePeriod_threshold = prctile(memoryPrecision_r2_baselinePeriod',shufflePrctileThreshold,1);
meta_r2_baselinePeriod_threshold = prctile(meta_r2_baselinePeriod',shufflePrctileThreshold,1);


temp1 = 2;

memoryPrecision_r2_baselinePeriod_threshold = memoryPrecision_r2_baselinePeriod_threshold * temp1;
meta_r2_baselinePeriod_threshold = meta_r2_baselinePeriod_threshold * temp1;



%% Plot

meta_color = [252,141,98]/255;
memory_color = [141,160,203]/255;


if if_plot == 1
    close all
    if if_plot_multiPanel == 1
        backgounrdColor = [1 1 1];
        
        fig = figure('Name','asd','NumberTitle','off');
        %set(gcf,'Position',[10 50 1650 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 50 720 200*0.9]);
        t = tiledlayout(1,5,'TileSpacing','Compact','Padding','Compact');
        
        t.Title.String = sprintf('Meta-memory and memory precision dynamics, %s',FOVName_currentFOV2);
        
        t.Title.FontSize = 12;
        t.Title.Interpreter = 'none';
        
        
        tempA_r2_all = meta_r2_all;
        tempB_r2_all = memoryPrecision_r2_all;
        
        if if_p_shuffle0_baseline1 == 0
            tempA_r2_all_shuffled = meta_r2_all_shuffled_threshold;
            tempB_r2_all_shuffled = memoryPrecision_r2_all_shuffled_threshold;
        elseif if_p_shuffle0_baseline1 == 1
            tempA_r2_all_shuffled = meta_r2_baselinePeriod_threshold;
            tempB_r2_all_shuffled = memoryPrecision_r2_baselinePeriod_threshold;
        end
        
        
        [yA_min,yA_max] = bounds([tempA_r2_all{1},tempA_r2_all{2},tempA_r2_all{3},tempA_r2_all{4},tempA_r2_all{5}]);
        [yB_min,yB_max] = bounds([tempB_r2_all{1},tempB_r2_all{2},tempB_r2_all{3},tempB_r2_all{4},tempB_r2_all{5}]);
        
        y_min = min([yA_min,yB_min]);
        y_max = max([yA_max,yB_max]);
        
        y_min = y_min-(y_max-y_min)*0.07;
        y_max = y_max+(y_max-y_min)*0.07;
        
        for loopCount=0:4
            
            tempA_r2_crossTime = tempA_r2_all{loopCount+1};
            tempB_r2_crossTime = tempB_r2_all{loopCount+1};
            
            if if_p_shuffle0_baseline1 == 0
                tempA_r2_crossTime_shuffled = tempA_r2_all_shuffled{loopCount+1};
                tempB_r2_crossTime_shuffled = tempB_r2_all_shuffled{loopCount+1};
            elseif if_p_shuffle0_baseline1 == 1
                tempA_r2_crossTime_shuffled = tempA_r2_all_shuffled;
                tempB_r2_crossTime_shuffled = tempB_r2_all_shuffled;
            end
            if loopCount==0
                temp_str = 'Baseline';
                nexttile
                temp_interval = baselinePeriod_interval;
            elseif loopCount==1
                temp_str = 'Sample, length1';
                nexttile
                temp_interval = length1_sample_interval;
            elseif loopCount==2
                temp_str = 'Sample, length2';
                nexttile
                temp_interval = length2_sample_interval;
            elseif loopCount==3
                temp_str = 'Sample, length3';
                nexttile
                temp_interval = length3_sample_interval;
            elseif loopCount==4
                temp_str = 'Delay1';
                nexttile
                temp_interval = decisionPeriodA_interval;
            end
            
            
            range_crossTime = 1:size(tempA_r2_crossTime,2);
            x = range_crossTime;
            
            
            %y = tempA_r2_crossTime_shuffled;
            %plot(x,y,'color',meta_color,'linewidth',1,'LineStyle','--');
            %hold on
            
            %y = tempB_r2_crossTime_shuffled;
            %plot(x,y,'color',memory_color,'linewidth',1,'LineStyle','--');
            %hold on
            
            
            h_line = [];
            
            y = tempA_r2_crossTime;
            h_line = [h_line plot(x,y,'color',meta_color,'linewidth',1)]; %#ok<*AGROW>
            hold on
            
            y = tempB_r2_crossTime;
            h_line = [h_line plot(x,y,'color',memory_color,'linewidth',1)]; %#ok<*AGROW>
            hold on
            
            if if_p_shuffle0_baseline1 == 0
                scatter(x(tempA_r2_crossTime>tempA_r2_crossTime_shuffled),y_min+(y_max-y_min)*0.05,8,meta_color,'*');
                scatter(x(tempB_r2_crossTime>tempB_r2_crossTime_shuffled),y_min+(y_max-y_min)*0.12,8,memory_color,'*');
            
            elseif if_p_shuffle0_baseline1 == 1
%                 scatter(x(tempA_r2_crossTime>meta_r2_baselinePeriod_threshold),y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                 scatter(x(tempB_r2_crossTime>memoryPrecision_r2_baselinePeriod_threshold),y_min+(y_max-y_min)*0.12,8,memory_color,'*');

                x;
                temp1 = 2;%0.5
                for temptempi=1:length(x)
                    [~,temp_p] = ttest2(tempA_r2_crossTime(x(temptempi:end)),meta_r2_baselinePeriod*temp1,'Tail','right');
                    %[~,temp_p] = ttest2(tempA_r2_crossTime(x(temptempi:end)),meta_r2_baselinePeriod);                    
                    if temp_p < 0.001
                        break
                    end
                end
                temptempi;
                scatter(x(temptempi:end),y_min+(y_max-y_min)*0.05,8,meta_color,'*');

                for temptempi=1:length(x)
                    [~,temp_p] = ttest2(tempB_r2_crossTime(x(temptempi:end)),memoryPrecision_r2_baselinePeriod*temp1,'Tail','right');
                    %[~,temp_p] = ttest2(tempB_r2_crossTime(x(temptempi:end)),memoryPrecision_r2_baselinePeriod);                    
                    if temp_p < 0.001
                        break
                    end
                end
                temptempi;
                scatter(x(temptempi:end),y_min+(y_max-y_min)*0.12,8,memory_color,'*');
                
                
                
                
%                 if loopCount==1
%                     [~,tempA_p_1] = ttest2(tempA_r2_crossTime(1:6),meta_r2_baselinePeriod,'Tail','right');
%                     [~,tempB_p_1] = ttest2(tempB_r2_crossTime(1:6),memoryPrecision_r2_baselinePeriod,'Tail','right');                
%                     if tempA_p_1 < 0.001
%                         scatter(4,y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                     end                    
%                     if tempB_p_1 < 0.001
%                         scatter(4,y_min+(y_max-y_min)*0.12,8,memory_color,'*');
%                     end
%                     
%                     [~,tempA_p_2] = ttest2(tempA_r2_crossTime(7:end),meta_r2_baselinePeriod,'Tail','right');
%                     [~,tempB_p_2] = ttest2(tempB_r2_crossTime(7:end),memoryPrecision_r2_baselinePeriod,'Tail','right');                 
%                     if tempA_p_2 < 0.001
%                         scatter(12,y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                     end                    
%                     if tempB_p_2 < 0.001
%                         scatter(12,y_min+(y_max-y_min)*0.12,8,memory_color,'*');
%                     end
% 
%                 elseif loopCount==2
%                     [~,tempA_p_1] = ttest2(tempA_r2_crossTime(1:6),meta_r2_baselinePeriod,'Tail','right');
%                     [~,tempB_p_1] = ttest2(tempB_r2_crossTime(1:6),memoryPrecision_r2_baselinePeriod,'Tail','right');                
%                     if tempA_p_1 < 0.001
%                         scatter(4,y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                     end                    
%                     if tempB_p_1 < 0.001
%                         scatter(4,y_min+(y_max-y_min)*0.12,8,memory_color,'*');
%                     end
%                     
%                     [~,tempA_p_2] = ttest2(tempA_r2_crossTime(7:18),meta_r2_baselinePeriod,'Tail','right');
%                     [~,tempB_p_2] = ttest2(tempB_r2_crossTime(7:18),memoryPrecision_r2_baselinePeriod,'Tail','right');                   
%                     if tempA_p_2 < 0.001
%                         scatter(12,y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                     end                    
%                     if tempB_p_2 < 0.001
%                         scatter(12,y_min+(y_max-y_min)*0.12,8,memory_color,'*');
%                     end      
%                     
%                     
%                     
%                     [~,tempA_p_3] = ttest2(tempA_r2_crossTime((1:6)+18),meta_r2_baselinePeriod,'Tail','right');
%                     [~,tempB_p_3] = ttest2(tempB_r2_crossTime((1:6)+18),memoryPrecision_r2_baselinePeriod,'Tail','right');                  
%                     if tempA_p_3 < 0.001
%                         scatter(4+18,y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                     end                    
%                     if tempB_p_3 < 0.001
%                         scatter(4+18,y_min+(y_max-y_min)*0.12,8,memory_color,'*');
%                     end
%                     
%                     [~,tempA_p_4] = ttest2(tempA_r2_crossTime((7+18):36),meta_r2_baselinePeriod,'Tail','right');
%                     [~,tempB_p_4] = ttest2(tempB_r2_crossTime((7+18):36),memoryPrecision_r2_baselinePeriod,'Tail','right');                   
%                     if tempA_p_4 < 0.001
%                         scatter(12+18,y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                     end                    
%                     if tempB_p_4 < 0.001
%                         scatter(12+18,y_min+(y_max-y_min)*0.12,8,memory_color,'*');
%                     end                      
%                     
%                     
%                     
%                 elseif loopCount==3
%                     [~,tempA_p_1] = ttest2(tempA_r2_crossTime(1:6),meta_r2_baselinePeriod,'Tail','right');
%                     [~,tempB_p_1] = ttest2(tempB_r2_crossTime(1:6),memoryPrecision_r2_baselinePeriod,'Tail','right');                   
%                     if tempA_p_1 < 0.001
%                         scatter(4,y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                     end                    
%                     if tempB_p_1 < 0.001
%                         scatter(4,y_min+(y_max-y_min)*0.12,8,memory_color,'*');
%                     end
%                     
%                     [~,tempA_p_2] = ttest2(tempA_r2_crossTime(7:18),meta_r2_baselinePeriod,'Tail','right');
%                     [~,tempB_p_2] = ttest2(tempB_r2_crossTime(7:18),memoryPrecision_r2_baselinePeriod,'Tail','right');                 
%                     if tempA_p_2 < 0.001
%                         scatter(12,y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                     end                    
%                     if tempB_p_2 < 0.001
%                         scatter(12,y_min+(y_max-y_min)*0.12,8,memory_color,'*');
%                     end      
%                     
%                     
%                     
%                     [~,tempA_p_3] = ttest2(tempA_r2_crossTime((1:6)+18),meta_r2_baselinePeriod,'Tail','right');
%                     [~,tempB_p_3] = ttest2(tempB_r2_crossTime((1:6)+18),memoryPrecision_r2_baselinePeriod,'Tail','right');                   
%                     if tempA_p_3 < 0.001
%                         scatter(4+18,y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                     end                    
%                     if tempB_p_3 < 0.001
%                         scatter(4+18,y_min+(y_max-y_min)*0.12,8,memory_color,'*');
%                     end
%                     
%                     [~,tempA_p_4] = ttest2(tempA_r2_crossTime((7+18):36),meta_r2_baselinePeriod,'Tail','right');
%                     [~,tempB_p_4] = ttest2(tempB_r2_crossTime((7+18):36),memoryPrecision_r2_baselinePeriod,'Tail','right');                 
%                     if tempA_p_4 < 0.001
%                         scatter(12+18,y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                     end                    
%                     if tempB_p_4 < 0.001
%                         scatter(12+18,y_min+(y_max-y_min)*0.12,8,memory_color,'*');
%                     end          
%                     
%                     
%                     
%                     
%                     
%                     [~,tempA_p_5] = ttest2(tempA_r2_crossTime((1:6)+36),meta_r2_baselinePeriod,'Tail','right');
%                     [~,tempB_p_5] = ttest2(tempB_r2_crossTime((1:6)+36),memoryPrecision_r2_baselinePeriod,'Tail','right');             
%                     if tempA_p_5 < 0.001
%                         scatter(4+36,y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                     end                    
%                     if tempB_p_5 < 0.001
%                         scatter(4+36,y_min+(y_max-y_min)*0.12,8,memory_color,'*');
%                     end
%                     
%                     [~,tempA_p_6] = ttest2(tempA_r2_crossTime((7+36):54),meta_r2_baselinePeriod,'Tail','right');
%                     [~,tempB_p_6] = ttest2(tempB_r2_crossTime((7+36):54),memoryPrecision_r2_baselinePeriod,'Tail','right');                   
%                     if tempA_p_6 < 0.001
%                         scatter(12+36,y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                     end                    
%                     if tempB_p_6 < 0.001
%                         scatter(12+36,y_min+(y_max-y_min)*0.12,8,memory_color,'*');
%                     end                      
%                     
%                     
%                 elseif loopCount == 4
%                     [~,tempA_p] = ttest2(tempA_r2_crossTime,meta_r2_baselinePeriod,'Tail','right');
%                     [~,tempB_p] = ttest2(tempB_r2_crossTime,memoryPrecision_r2_baselinePeriod,'Tail','right');
%                     
%                     if tempA_p < 0.001
%                         scatter(x(round(length(x)/2)),y_min+(y_max-y_min)*0.05,8,meta_color,'*');
%                     end
%                     if tempB_p < 0.001
%                         scatter(x(round(length(x)/2)),y_min+(y_max-y_min)*0.12,8,memory_color,'*');
%                     end
%                 end
            end
            
            
            
            if loopCount==1 || loopCount==2 || loopCount==3 || loopCount== 4
                for tempi=1:(length(temp_interval)-1)
                    plot(temp_interval(tempi)*[1 1],[y_min y_max],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                if loopCount==1 || loopCount==2 || loopCount==3
                    for tempi=1:(length(temp_interval)-1)
                        x = [temp_interval(1+tempi-1) y_min;...
                            temp_interval(1+tempi-1) y_max;...
                            temp_interval(1+tempi-1)+6 y_max;...
                            temp_interval(1+tempi-1)+6 y_min];
                        y = [1 2 3 4];
                        patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.1,'EdgeColor','none');
                        hold on
                    end
                end
            end
            
            
            if loopCount == 0
                le = legend(h_line,...
                    sprintf('Meta'),...
                    sprintf('Memory'),...
                    'Location','northeast','fontsize',9,'NumColumns',1);
                le.ItemTokenSize = ones(1,4)*10;%20
                le.Color = backgounrdColor;
            end
            
            if loopCount==0
                xticks(temp_interval(1:end-2));
                xticklabels({'Fixation',''});
            elseif loopCount==1
                xticks(temp_interval(1:end-1));
                xticklabels({'T1'});
            elseif loopCount==2
                xticks(temp_interval(1:end-1));
                xticklabels({'T1','T2'});
            elseif loopCount==3
                xticks(temp_interval(1:end-1));
                xticklabels({'T1','T2','T3'});
            elseif loopCount==4
                xticks(temp_interval(1:end-1));
                xticklabels({'Delay1','ChoiceCue'});
            end
            
            
            xlim([range_crossTime(1) range_crossTime(end)]);
            ylim([y_min y_max]);
            
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 11)
            set(gca,'color',backgounrdColor);
            set(gca,'box','off');% 取消右、上边框
            if loopCount == 0
                ylabel('r2', 'FontSize', 12, 'FontWeight', 'bold');
            end
            xtickangle(0);
            
            temp_title = title(sprintf('%s',temp_str),'FontSize',11,'FontWeight','normal');%'bold
            temp_title.Interpreter = 'none';
            
            
        end
        
    end
    
    %end
    
end

%% End