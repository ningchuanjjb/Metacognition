function modifyCellProb_dff_v77(path)

%temp_spkEndIndex
%temptempIndex_raw

% edit in 20230616
% (1) add spike-sorting-like way to select roi (calcium pulse decay coeff)
% (2) add threshold_validSpkCount_max_per10000 to permit high acative roi that couldn't survive in (1)
% (3) using both F and dF / F to sort roi now

if(~exist('path','var'))
    path = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230620A_Ding_Site05E\Result20230620T231851\plane0';                              
end

a = 1;

fileName_Fall = 'Fall.mat';
fileName_iscell = 'iscell.npy';

if_speedup = 1;
if_profile = 0;

threshold_compact = 1.15;
threshold_aspect_ratio = 1.15;
threshold_cellProb = 0.10;
threshold_npix_smallCell = 6;%10-->15-->10-->6-->10-->6(20230619)

threshold_F = 75;%300-->200-->150-->100-->125-->75
threshold_F_forBound = 150;

threshold_peak_SNR = 2.6;%3-->2.7-->2.6
% threshold_peak_SNR_forSpk = 2.7;% Very sensitive! 2-->2.4-->2.7
% threshold_peak_SNR_forSpk = 2.2;% Very sensitive! 2-->2.4-->2.7-->2.2(20230615)

% threshold_bigF_init = 8000;
% threshold_bigF = 800;%800-->600-->800-->800
% threshold_validSpkCount_max_per10000 = 0.95;%% Very sensitive! 1.75-->1.5-->1.3-->1.1-->0.7
% threshold_peak_SNR_forSpk = 2.7;%% Very sensitive! 2.2-->2-->1.8-->2.2-->2.4-->2.6-->2.7
% threshold_validSpkCount_max_per10000 = 2.15;%0.95-->1.05-->1.15-->1.25-->1.4-->1.6-->1.5-->1.3-->1.5-->1.8-->2.0-->2.15
threshold_validSpkCount_max_per10000 = 1;%1.5-->1.8-->2.0-->2.15-->2.3-->2.2-->2.0-->1.8-->10-->20-->17-->18.5
threshold_peak_SNR_forSpk = 3.0;%% 2.7-->2.8-->2.9-->3.0-->4.5
% threshold_peak_dff = 1.5;%3-->1.5
% threshold_fall_dff = 1.6;%1.5-->1.3-->1.6
% threshold_rise_dff = 1;%0.5-->0.4-->0.6-->1

% threshold_r2_calciumPulse = 0.5;%0.85-->0.5-->0.7-->0.8-->0.5
% threshold_calciumPulse_decayCoeff = 0.6;%0.7-->0.5-->0.45-->0.7-->0.65-->0.675-->0.6-->0.65-->0.6
% threshold_r2_calciumPulse = 0.85;%0.5-->0.6-->0.7-->0.8-->0.82-->0.84-->0.85-->0.86-->0.85
% threshold_calciumPulse_decayCoeff = 0.75;%0.6-->0.65-->0.67-->0.69-->0.70-->0.71-->0.72-->0.75
threshold_r2_calciumPulse = 0.895;%0.85-->0.84-->0.82-->0.80-->0.78-->0.83-->0.82-->0.81-->0.80-->0.77-->0.79-->0.82-->0.805
threshold_calciumPulse_decayCoeff = 0.885;%0.75-->0.74-->0.72-->0.70-->0.68-->0.73-->0.72-->0.71-->0.70-->0.67-->0.69-->0.72-->0.705

compensateValue_forActiveCell = 0.1;%0.1-->0.2-->0.15-->0.07

threshold_spks_start = 10;% Very sensitive! 11-->10
threshold_spks_start2 = 200;%15-->19-->50-->30-->100-->200, edit in 20230615

threshold_spkLateActivity_stdCoeff = 1.4;% Very sensitive! 1.8-->1.5-->1-->1.1-->1.4
threshold_spks_duration = 10;% 8 --> 20, edit in 20230616, 20-->10
threshold_spks_riseTime = 3;% add in 20230615
threshold_spks_fallMinusRise = 4;%4


threshold_validSpkCount_per10000 = 0.01;%0.2-->0.1-->0.01
threshold_validSpkCount_min = 2;%3-->2

threshold_median = 800;
threshold_npix_bigCell = 200;%200-->150-->200
% threshold_npix_hugeCell = 600;

threshold_rescueNum = 0;%100-->80-->70-->42-->0

fullFileName_Fall = [path,'\',fileName_Fall];
fullFileName_iscell = [path,'\',fileName_iscell];

% load(fullFileName_Fall,'stat','F','Fneu','spks','iscell');
load(fullFileName_Fall,'stat','F','Fneu','spks','iscell','F_dff_raw');

F_dff = F_dff_raw;
clear F_dff_raw

toal_validSpkCount = 0; %#ok<*NASGU>
toal_validSpkCount_boolIndex = false(1, length(stat));
threshold_validSpkCount = floor(threshold_validSpkCount_per10000*floor(size(F,2)/10000));
threshold_validSpkCount = max(threshold_validSpkCount, threshold_validSpkCount_min);

threshold_validSpkCount_max = floor(threshold_validSpkCount_max_per10000*floor(size(F,2)/10000));


threshold_xypix_ratio = 10;%10-->8
threshold_xypix_ratio_xyBoundDis = 30;%20

xypix_ratio = zeros(1, length(stat));
xypix_ratio_badBoolIndex = false(1, length(stat));
xpix_median = zeros(1, length(stat));

% windowSize = 2; 
% filter_b = (1/windowSize)*ones(1,windowSize);
% filter_a = 1;

bigCell_bool = false(1, length(stat));
brightCell_bool = false(1, length(stat));

manySpkCountCell_bool = false(1, length(stat));
manySpkCountCell_bool_all = false(1, length(stat));
cell_npix = zeros(1, length(stat));

highFneu_bool = false(1, length(stat));
negativeFpeak_bool = false(1, length(stat));


if if_profile == 1
    profile on
end
t0 = tic;

% windowSize = 10;
% F_raw = F;
% F = smoothdata(F,2,'gaussian',windowSize);


r2_calciumPulse = zeros(1,length(stat));
calciumPulse_decayCoeff = zeros(1,length(stat));


peak_F = zeros(1,length(stat));
peak_F_dff = zeros(1,length(stat));

T = size(F_dff,2);

%% First xypix_ratio_badIndex_suite2p
for tempi=1:length(stat)
        
    temp_stat = stat{tempi};        
        
    % Extract very long and thin roi and make cellProb 0
    [min_xpix,max_xpix] = bounds(temp_stat.xpix);
    [min_ypix,max_ypix] = bounds(temp_stat.ypix);
    range_ypix = double(max_ypix-min_ypix);
    range_xpix = double(max_xpix-min_xpix);
    
    xpix_median(tempi) = median(temp_stat.xpix); %#ok<*SAGROW>
    
    if range_ypix > range_xpix
        temp_max_range_xypix = range_ypix;
        temp_min_range_xypix = range_xpix; %#ok<*NASGU>
        temp_min_range_xyMed = median(temp_stat.xpix);
        
        temp_min_range_xypix = diff(prctile(double(temp_stat.xpix),[25,75]));%[15,85]
    else
        temp_max_range_xypix = range_xpix;
        temp_min_range_xypix = range_ypix;
        temp_min_range_xyMed = median(temp_stat.ypix);
        
        temp_min_range_xypix = diff(prctile(double(temp_stat.ypix),[25,75]));
    end
    temp_min_range_xypix = max([temp_min_range_xypix, 3]); % temp_min_range_xypix should more than 3
    temp_xypix_ratio = temp_max_range_xypix / temp_min_range_xypix;
    
    temp_min_range_xyBoundDis = min(temp_min_range_xyMed,513-temp_min_range_xyMed);
    
    if temp_xypix_ratio > threshold_xypix_ratio ...
            && temp_min_range_xyBoundDis < threshold_xypix_ratio_xyBoundDis
        if range_ypix > range_xpix % only if range_ypix > range_xpix
            xypix_ratio_badBoolIndex(tempi) = true;
        end
    end
    xypix_ratio(tempi) = temp_xypix_ratio;
    
end


temp_binSize = 15;%6-->15

F_2 = spks(:,1:floor(size(F,2)/temp_binSize)*temp_binSize);% F --> spks
%F_2 = F(:,1:floor(size(F,2)/temp_binSize)*temp_binSize);
F_2_3d = reshape(F_2,[size(F_2,1) temp_binSize size(F_2,2)/temp_binSize]);
F_2_bin = squeeze(mean(F_2_3d,2));
F_2_bin_T = F_2_bin';

corr_threshold3 = 100;%0.25-->0.2-->100

xypix_ratio_badIndex = find(xypix_ratio_badBoolIndex==true);
xypix_ratio_badIndex_suite2p = xypix_ratio_badIndex - 1;
fprintf('xypix_ratio_badIndex_suite2p = %s.\n',num2str(xypix_ratio_badIndex_suite2p));
xypix_ratio_badBoolIndex_1 = xypix_ratio_badBoolIndex;
xypix_ratio_badIndex_1 = xypix_ratio_badIndex;

%% Second xypix_ratio_badIndex_suite2p
for tempIndex=1:length(xypix_ratio_badIndex_1)
    tempi = xypix_ratio_badIndex_1(tempIndex);
    
    temp_stat = stat{tempi};     
    temp_min_range_xyMed = median(temp_stat.xpix);
    
    if temp_min_range_xyMed < threshold_xypix_ratio_xyBoundDis
        temp_boundBool_left0_right1 = false;        
        temp_boundBoolIndex = xpix_median < threshold_xypix_ratio_xyBoundDis;        
    else
        temp_boundBool_left0_right1 = true;        
        temp_boundBoolIndex = (513-xpix_median) < threshold_xypix_ratio_xyBoundDis;
    end   
    
    tempID_corr_F= corr(F_2_bin_T(:,tempi),F_2_bin_T); %#ok<*PFBNS>
    tempID_corrBool_F_raw = tempID_corr_F >= corr_threshold3;    
    tempID_corrBool_F = tempID_corrBool_F_raw & temp_boundBoolIndex;
    
    tempCellID_toBeClear = find(tempID_corrBool_F==true);
    tempCellID_suite2p_toBeClear = tempCellID_toBeClear - 1;
    
    xypix_ratio_badBoolIndex(tempCellID_toBeClear) = true;
    
end

xypix_ratio_badIndex = find(xypix_ratio_badBoolIndex==true);
xypix_ratio_badIndex_suite2p = xypix_ratio_badIndex - 1;
fprintf('second xypix_ratio_badIndex_suite2p (num) = %d.\n',length(xypix_ratio_badIndex_suite2p));
xypix_ratio_badBoolIndex_2 = xypix_ratio_badBoolIndex;
xypix_ratio_badIndex_2 = xypix_ratio_badIndex;


% for tempi=1:40
% for tempi=1:length(stat)
parfor tempi=1:length(stat)
  
%     tempi = 73 + 1; %#ok<*FXSET>
    if tempi == 3+1
        a = 1;
    end
    
    temp_spkEndIndex = [];    
    
    temp_stat = stat{tempi};
    
    cell_npix(tempi) = temp_stat.npix;
    
    % Extract dendrite and make cellProb 1
    if (temp_stat.compact>threshold_compact || temp_stat.aspect_ratio>threshold_aspect_ratio) ...
            && iscell(tempi, 2) < threshold_cellProb %#ok<*PFIIN>
        %iscell(tempi, 2) = 0.1;%1
        iscell(tempi, 2) = 1.1;%1
    end    
    
    % Extract very small roi and make cellProb 0
    if temp_stat.npix < threshold_npix_smallCell
        iscell(tempi, 2) = 0;
    end
        
    tempID_F_dff = F_dff(tempi,:);
    tempID_F_dff_mean = mean(tempID_F_dff);
    %peak_tempID_F_dff = max(tempID_F_dff);
    
    tempID_F = F(tempi,:);
    tempID_Fneu = Fneu(tempi,:);
    tempID_F_mean = mean(tempID_F);
    tempID_Fneu_mean = mean(tempID_Fneu);
    peak_tempID_F = max(tempID_F);
    
    peak_F(tempi) = peak_tempID_F;
    %peak_F_dff(tempi) = peak_tempID_F_dff;
    
    peak_SNR = peak_tempID_F/max(tempID_F_mean,tempID_Fneu_mean);
    std_tempID_F = std(tempID_F);
    std_tempID_Fneu = std(tempID_Fneu);
    std_tempID_F_dff = std(tempID_F_dff);
    
    tempID_F_median = median(tempID_F);
    %tempID_F(tempID_F<tempID_F_median) = tempID_F_median;%?????
    
    windowSize = 5;% 10-->5
    tempID_F_smooth = smoothdata(tempID_F,2,'gaussian',windowSize);
    tempID_F_dff_smooth = smoothdata(tempID_F_dff,2,'gaussian',windowSize);
%     tempID_F_smooth = smoothdata(tempID_F,2,'sgolay',windowSize);
%     tempID_F_dff_smooth = smoothdata(tempID_F_dff,2,'sgolay',windowSize);
    
    
    %peak_tempID_F_dff = max(tempID_F_dff);
    peak_tempID_F_dff = max(tempID_F_dff_smooth);
    peak_F_dff(tempi) = peak_tempID_F_dff;
    
    
    %if peak_tempID_F < threshold_F || peak_SNR < threshold_peak_SNR || tempID_F_mean < (tempID_Fneu_mean-std_tempID_Fneu)
    if peak_tempID_F < threshold_F || peak_SNR < threshold_peak_SNR
        iscell(tempi, 2) = 0;
    end
    if xypix_ratio_badBoolIndex(tempi) == true
        if peak_tempID_F < threshold_F_forBound
            iscell(tempi, 2) = 0;
        end
    end
    
    if tempID_F_mean < (tempID_Fneu_mean-std_tempID_Fneu)
    %if (tempID_F_mean-tempID_Fneu_mean) < 3*std_tempID_Fneu
        highFneu_bool(tempi) = true;
    end
    
    a = 1;
    
    temp_coeff = 4;%3-->4
    temp_frameCount = 50;%50    
    temp_diff = tempID_F_mean-temp_coeff*std_tempID_F;
    for temptempi=temp_coeff:-1:1
        temp_coeff = temptempi;
        if temp_diff > 0
            break
        end
        temp_diff = temp_diff + std_tempID_F;
    end
    
    temp_negativeFpeak_frameBoolIndex = tempID_F_smooth < (tempID_F_mean-temp_coeff*std_tempID_F);
    %sum(temp_negativeFpeak_frameBoolIndex)
    %temp_negativeFpeak_frameIndex = find(temp_negativeFpeak_frameBoolIndex==true);    
    %temp1 = tempID_F(temp_negativeFpeak_frameBoolIndex)';
    %temp2 = (temp_negativeFpeak_frameIndex-1)';    
    
    temptemp_F = tempID_F_smooth;
    temptemp_F(~temp_negativeFpeak_frameBoolIndex) = tempID_F_median;
    
    [~,temp_negativeFvalleys_frameIndex] = findpeaks(-temptemp_F); % Time-consuming!!!    
    temp_negativeFvalleys = temptemp_F(temp_negativeFvalleys_frameIndex);
    
    %if length(temp_negativeFvalleys) > temp_frameCount
    %    negativeFpeak_bool(tempi) = true;
    %end
    
    if if_speedup == 1
        if iscell(tempi, 2) < eps
            continue
        end
    end
    tempID_spks = spks(tempi,:);
    temp_validSpkFlag = false;
    temp_validSpkCount = 0;
    
    tempCalciumPulse = zeros((length(tempID_spks)-threshold_spks_start2), threshold_spks_start2);
    tempCalciumPulse_templength = zeros((length(tempID_spks)-threshold_spks_start2), 1);
    
    last_spkEndIndex = 0;
    tempj = 0;
    tempLoopCount = 0;
    %for tempj=1:(length(tempID_spks)-threshold_spks_start2)
    while true
        
        if isempty(last_spkEndIndex) == true
            last_spkEndIndex = 0;
        end
        
        tempj = tempj + 1 + floor(last_spkEndIndex/2);
        
        if tempj > (length(tempID_spks)-threshold_spks_start2)
            break
        end
        
        tempj_minus1 = tempj - 1;           
        
        %if tempID_spks(tempj9 > tempID_F_mean
        %if tempID_spks(tempj) > tempID_Fneu_mean/2
        %if tempID_spks(tempj) > tempID_Fneu_mean/3
        %if tempID_spks(tempj) > tempID_Fneu_mean/4
        if tempID_spks(tempj) > 1
            tempLoopCount = tempLoopCount + 1;
            
%             % To get temp_temp0Flag
%             % Whether spk late activity more than baseline
%             temp0 = tempID_F(tempj+(0:threshold_spks_start));
%             if floor(length(temp0)/2)*2 ~= length(temp0)
%                 temp0(end) = [];
%             end
%             
%             temp0D = filter(filter_b,filter_a,temp0);
%             temp0D(1) = temp0(1);
% 
%             temp_temp0Flag = false;
%             if min(temp0D)>(tempID_F_mean+threshold_spkLateActivity_stdCoeff*std_tempID_F)
%                 temp_temp0Flag = true;
%             end
            
            % To get temp_tem10Flag
            % Whether spk rise faster than fall
            %temp1 = tempID_F(tempj+(0:threshold_spks_start2));
            %temp1 = tempID_F_smooth(tempj+(0:threshold_spks_start2));
            temp1 = tempID_F_dff_smooth(tempj+(0:threshold_spks_start2));
            if floor(length(temp1)/2)*2 ~= length(temp1)
                temp1(end) = [];
            end
            
            temp1D = temp1;
            %temp1D_envelope = abs(hilbert(temp1D));
            %[temp1D_envelope,~] = envelope(temp1D,10,'analytic');
            %[temp1D_envelope,~] = envelope(temp1D,1,'peak');%10
%             [temp_peaks,temp_locs] = findpeaks(temp1D); % Time-consuming!!!
%             temp_peaks = [temp1D(1) temp_peaks temp1D(end)]; %#ok<*AGROW>
%             temp_locs = [1 temp_locs size(temp1D,2)];
%             temp1D_envelope = interp1(temp_locs,temp_peaks,temp_locs(1):temp_locs(end));

            
            temp1D_valid = temp1D;
            %temp1D_valid = temp1D_envelope;
            
            
            %temp1D_raw = tempID_F_smooth(tempj+(0:threshold_spks_start2));
            temp1D_raw = tempID_F_dff(tempj+(0:threshold_spks_start2));
            if floor(length(temp1D_raw)/2)*2 ~= length(temp1D_raw)
                temp1D_raw(end) = [];
            end
            
            
            %temp1 = tempID_F(tempj+(0:threshold_spks_start2));
            %temp1 = tempID_F_dff_smooth(tempj+(0:threshold_spks_start2));
            %if floor(length(temp1)/2)*2 ~= length(temp1)
            %    temp1(end) = [];
            %end
            
            %[temp_max,temp_peakIndex] = max(temp1D);
            %[temp_max,temp_peakIndex] = max(temp1D_valid(1:round(size(temp1D_valid,2)/2)));
            [temp_max,temp_peakIndex] = max(temp1D_valid);
            
            % [valleys2,locs2] = findpeaks(temp_max-temp1D); % Time-consuming!!!
            a = 1;
            

            
            temp_diff = diff(temp1D_valid);
            tempBoolIndex_nagative = temp_diff < 0;
            tempBoolIndex_positive = temp_diff >= 0;
            temp_diff2 = diff(tempBoolIndex_positive);
            temp_valleyBoolIndex = temp_diff2 == 1;
            locs = find(temp_valleyBoolIndex==true)+1;
            valleys = temp_max - temp1D_valid(locs);
            
            
            temp_validIndex_afterPeak = find(temp_peakIndex<locs,1);
            %temp_validIndex_afterPeak = 1;
            
            valleys_prctile = median(valleys(temp_validIndex_afterPeak:end));
            tempBoolIndex = valleys>valleys_prctile;
            tempBoolIndex(1:(temp_validIndex_afterPeak-1)) = false;
            
            tempValleyIndex = find(tempBoolIndex==true,1);
            trueValleyIndex = tempValleyIndex;
            trueValleyIndex2 = tempValleyIndex;
            
            for tempk=temp_validIndex_afterPeak:(length(locs)-2)
                locs(tempk);
                %if temp1D(locs(tempk)) < 0 ...
                %       || temp1D(temp_peakIndex)/temp1D(locs(tempk)) > (threshold_peak_SNR_forSpk*2)
                if temp1D(locs(tempk)) < 0 ...
                  || (temp1D(temp_peakIndex)-temp1D(locs(tempk))) > threshold_peak_SNR_forSpk*std_tempID_F_dff
                %if temp1D(locs(tempk)) < 0 ...
                %  || (temp1D(temp_peakIndex)-temp1D(locs(tempk))) > threshold_peak_SNR_forSpk*1.5*std_tempID_F_dff                
                %if (temp1D_valid(temp_peakIndex)-temp1D_valid(locs(tempk))) > threshold_peak_SNR_forSpk*std_tempID_F_dff
                    %& (temp1D(temp_peakIndex)-max(temp1D(locs(tempk):min(locs(tempk)+8,size(temp1D,2)) ))) > threshold_peak_SNR_forSpk*std_tempID_F_dff
                                
                    trueValleyIndex2 = tempk;
                    break
                end
            end
            
            trueValleyIndex = floor((tempValleyIndex + trueValleyIndex2)/2);
            %trueValleyIndex = trueValleyIndex2;
            %trueValleyIndex = tempValleyIndex;
            %trueValleyIndex = 5;
            
            for tempk=trueValleyIndex:length(locs)
                if locs(tempk) > temp_peakIndex                    
                    break
                end
            end
            trueValleyIndex = tempk;
            
            
            trueValleyLoc = locs(trueValleyIndex);
            trueValleyF = valleys(trueValleyIndex);
            
            temp_spkEndIndex = trueValleyLoc;
                
            temp_duration_rise = temp_peakIndex-1;
            temp_duration_fall = temp_spkEndIndex-temp_peakIndex;
            
            
            
            temp_temp1Flag = false;
            %if (temp_duration_fall-temp_duration_rise > threshold_spks_fallMinusRise) ...
            %        & (temp_duration_fall > threshold_spks_duration) ...
            %        & (temp_duration_rise >= threshold_spks_riseTime) %#ok<*AND2>
            %if true    
            %if temp_peakIndex < 180
            %if temp_peakIndex < 180 & temp_duration_fall > 0
            %if temp_peakIndex < 180 & isempty(temp_duration_fall)==false
            %if temp_peakIndex < 180 & isempty(temp_duration_fall)==false ...
            %        & (temp_duration_fall > threshold_spks_duration)
            if temp_peakIndex < 100 & isempty(temp_duration_fall)==false ...
                   & (temp_duration_fall > threshold_spks_duration)
                
                
                %if (temp1D(temp_peakIndex)/temp1D(temp_spkEndIndex)) > threshold_peak_SNR_forSpk
                %   temp_temp1Flag = true;
                %end   
                
% threshold_peak_dff;
% threshold_fall_dff;
% threshold_rise_dff;

                %if temp1D(temp_peakIndex) > threshold_peak_dff ...
                %if abs((temp1D(temp_peakIndex)/temp1D(1))) > threshold_peak_SNR_forSpk/2 ...
                %        & abs((temp1D(temp_peakIndex)/temp1D(temp_spkEndIndex))) > threshold_peak_SNR_forSpk ...
                %        & temp1D(temp_spkEndIndex) < threshold_fall_dff ...                        
                %        & temp1D(1) < threshold_rise_dff ...
                %        & temp1D_raw(temp_peakIndex) > threshold_F
                
                %if (temp1D(temp_peakIndex)-temp1D(1)) > threshold_peak_SNR_forSpk*std_tempID_F_dff ...
                %       & (temp1D(temp_peakIndex)-temp1D(temp_spkEndIndex)) > threshold_peak_SNR_forSpk*std_tempID_F_dff
                
                %if (temp1D(temp_peakIndex)-temp1D(temp_spkEndIndex)) > threshold_peak_SNR_forSpk*std_tempID_F_dff
                %if (temp1D_valid(temp_peakIndex)-temp1D_valid(temp_spkEndIndex)) > threshold_peak_SNR_forSpk*std_tempID_F_dff ...
                %       & (temp1D_valid(temp_peakIndex)-tempID_F_dff_mean) > threshold_peak_SNR_forSpk*std_tempID_F_dff
                %if (temp1D_valid(temp_peakIndex)-temp1D_valid(temp_spkEndIndex)) > threshold_peak_SNR_forSpk*std_tempID_F_dff ...
                %       & (temp1D_valid(temp_peakIndex)-tempID_F_dff_mean) > threshold_peak_SNR_forSpk*std_tempID_F_dff
                
                %if (temp1D(temp_peakIndex)-max(temp1D(temp_spkEndIndex:min(temp_spkEndIndex+8,size(temp1D,2))))) > threshold_peak_SNR_forSpk*std_tempID_F_dff ...
                %       & (temp1D(temp_peakIndex)-tempID_F_dff_mean) > threshold_peak_SNR_forSpk*std_tempID_F_dff...
                        %& (max(temp1D_raw)-tempID_F_mean) > threshold_peak_SNR_forSpk*std_tempID_F
                        
                        
                        
                %if (temp1D_valid(temp_peakIndex)-temp1D_valid(temp_spkEndIndex)) > threshold_peak_SNR_forSpk*std_tempID_F_dff ...
                %      & (temp1D_valid(temp_peakIndex)-tempID_F_dff_mean) > threshold_peak_SNR_forSpk*std_tempID_F_dff
                     
                %if (temp1D_valid(temp_peakIndex)-tempID_F_dff_mean) > (threshold_peak_SNR_forSpk-0.5)*std_tempID_F_dff
                if (temp1D_valid(temp_peakIndex)-tempID_F_dff_mean) > (threshold_peak_SNR_forSpk+0)*std_tempID_F_dff
                %if (temp1D_valid(temp_peakIndex)-tempID_F_dff_mean) > (threshold_peak_SNR_forSpk+1.5)*std_tempID_F_dff
                
                   temp_temp1Flag = true;
                    
                    %tempCalciumPulse(tempj,1:length(temp_peakIndex:temp_spkEndIndex)) = ...
                    %    temp1D_raw(temp_peakIndex:temp_spkEndIndex);
                    %tempCalciumPulse(tempj,1:length(temp_peakIndex:temp_spkEndIndex)) = ...
                    %   temp1D(temp_peakIndex:temp_spkEndIndex);             
                    tempCalciumPulse(tempj,1:length(temp_peakIndex:(end-1))) = ...
                        temp1D(temp_peakIndex:(end-1));   
                    
                    %if (temp_duration_fall-temp_duration_rise > threshold_spks_fallMinusRise) ...
                    %        & (temp_duration_fall > threshold_spks_duration) ...
                    %        & (temp_duration_rise >= threshold_spks_riseTime) %#ok<*AND2>
                    %    
                    %    tempCalciumPulse_templength(tempj) = temp_duration_fall;
                    %end
                    
                    tempCalciumPulse_templength(tempj) = temp_duration_fall;
                    %tempCalciumPulse_templength(tempj) = 30;
                end
            end
                        
            tempj_minus1; %#ok<*VUNUS>
            tempi;
            
            temptemptemp1 = 331700-20;
            %if tempj_minus1 == temptemptemp1
            if tempj_minus1 > temptemptemp1 - 20 ...
                    & tempj_minus1 < temptemptemp1 + 20 %#ok<*AND2>
                a = 1;
            end            
            
            % To determine whether there is a valid spk or not
            %if temp_temp0Flag == true && temp_temp1Flag == true
            if temp_temp1Flag == true
                temp_validSpkCount = temp_validSpkCount + 1;
                if temp_validSpkCount >= threshold_validSpkCount
                    temp_validSpkFlag = true;
                    %if if_speedup == 1
                    %    break
                    %end                    
                end
                last_spkEndIndex = temp_spkEndIndex;
            else
                last_spkEndIndex = 0;
            end
            
        else
            last_spkEndIndex = 0;
        end
    end
    
    
    a = 1;
    tempLoopCount;
    temp_sum = sum(tempCalciumPulse,2);
    temptempBoolIndex = temp_sum~=0;        
    temptempIndex_raw = find(temptempBoolIndex==true);
    temptempIndex = find(temptempBoolIndex==true);
    
    tempCalciumPulse_templength_valid = tempCalciumPulse_templength(temptempIndex,:);    
    
    tempCalciumPulse_valid = tempCalciumPulse(temptempIndex,:);
    tempCalciumPulse_length = zeros(length(temptempIndex),1);
    for tempkk=1:size(tempCalciumPulse_valid,1)
        tempCalciumPulse_length(tempkk) = find(tempCalciumPulse_valid(tempkk,:)==0,1);
    end
    %temp_validLength = min(tempCalciumPulse_length);
    temp_validLength = round(median(tempCalciumPulse_length));
    if isnan(temp_validLength)
        temp_validLength = 0;
    end
    
    if size(tempCalciumPulse_valid,1) > 1
        tempCalciumPulse_valid2 = tempCalciumPulse_valid(:,1:temp_validLength);        
        tempCalciumPulse_valid3 = tempCalciumPulse_valid2;
    else
        tempCalciumPulse_valid3 = tempCalciumPulse_valid(:,1:temp_validLength);
    end
    
    a = 1;
    
    temp_max = max(tempCalciumPulse_valid3,[],2);
    
    %temp_prctile = 0;
    temp_prctile = min(temp_max)-0.01;
    %temp_prctile = tempID_F_dff_mean+threshold_peak_SNR_forSpk*std_tempID_F_dff;        
    %temp_prctile = prctile(temp_max,80);%70
    
    temp_max_prctile = temp_max(temp_max>temp_prctile);
    temp_max_prctileBoolIndex = temp_max>temp_prctile;
    
    
    temp_validCalciumPluse_frameIndex = temptempIndex_raw(temp_max_prctileBoolIndex);
    %temp_validCalciumPluse_frameIndex = temptempIndex_raw;
    
    temp_validCalciumPluse_frameIndex_valid = ...
        temp_validCalciumPluse_frameIndex(tempCalciumPulse_templength_valid>threshold_spks_duration);
    %temp_validCalciumPluse_frameIndex_valid = temp_validCalciumPluse_frameIndex;
    
    tempCalciumPulse_mean = mean(tempCalciumPulse_valid3(temp_max_prctileBoolIndex,:));
    
    %tempCalciumPulse_mean = mean(tempCalciumPulse_valid3);
        
    
    a = 1;
    if sum(isnan(tempCalciumPulse_mean)) >= 1 ...
            || sum(tempCalciumPulse_mean>100000) >= 1 ...
            || length(tempCalciumPulse_mean) < 3 
            %|| length(temp_validCalciumPluse_frameIndex_valid) < 3
        
    else
        
        tempCalciumPulse_mean_smooth = smoothdata(tempCalciumPulse_mean,2,'gaussian',windowSize);
        
        [temp_min2, temp_max2] = bounds(tempCalciumPulse_mean_smooth); %#ok<*ASGLU>
        [valleys2,locs2] = findpeaks(temp_max2-tempCalciumPulse_mean_smooth); % Time-consuming!!!
        
        tempCalciumPulse_mean_raw = tempCalciumPulse_mean;
        if isempty(locs2) == false
            temp_validLength2 = locs2(1);
        else
            temp_validLength2 = 0;
        end
        
        temp_validLength_final = max(temp_validLength,temp_validLength2);
        tempCalciumPulse_mean = tempCalciumPulse_mean(1:temp_validLength_final);        
        
        if length(tempCalciumPulse_mean) >= 3        
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
            %opts.Lower = [0 0 min(y)-5];
            %opts.Upper = [50 1 min(y)+5];
            % opts.StartPoint = [1.55650486517981 -0.0909553835273515 0];
            opts.StartPoint = [100 0 0];
            [fitresult, gof] = fit(xData,yData,ft,opts);
            
            temp_r2_calciumPulse = gof.rsquare;
            %temp_calciumPulse_decayCoeff = exp(fitresult.b);
            temp_calciumPulse_decayCoeff = fitresult.b;
            
            r2_calciumPulse(tempi) = temp_r2_calciumPulse;
            calciumPulse_decayCoeff(tempi) = temp_calciumPulse_decayCoeff;
            
            threshold_r2_calciumPulse;
            threshold_calciumPulse_decayCoeff;
            
            %if temp_r2_calciumPulse < threshold_r2_calciumPulse
            %if (temp_r2_calciumPulse+temp_calciumPulse_decayCoeff) < (threshold_r2_calciumPulse+threshold_calciumPulse_decayCoeff)\
            if (temp_r2_calciumPulse+temp_calciumPulse_decayCoeff) < (threshold_r2_calciumPulse+threshold_calciumPulse_decayCoeff) ...
                    || temp_r2_calciumPulse < (threshold_r2_calciumPulse-compensateValue_forActiveCell)
            
                temp_validSpkFlag = false;
            end
            if temp_calciumPulse_decayCoeff < threshold_calciumPulse_decayCoeff
                temp_validSpkFlag = false;
            end
            
            %if temp_r2_calciumPulse >= 0.9 ...
            %        && temp_calciumPulse_decayCoeff >= 0.8
            %if temp_r2_calciumPulse >= threshold_r2_calciumPulse ...
            %       && temp_calciumPulse_decayCoeff >= threshold_calciumPulse_decayCoeff            
            %if temp_r2_calciumPulse >= threshold_r2_calciumPulse ...
            %      && temp_calciumPulse_decayCoeff >= threshold_calciumPulse_decayCoeff ...         
            %      && highFneu_bool(tempi) == false
            %if temp_r2_calciumPulse >= threshold_r2_calciumPulse ...
            %      && temp_calciumPulse_decayCoeff >= threshold_calciumPulse_decayCoeff ...         
            %      && negativeFpeak_bool(tempi) == false 
            %if temp_r2_calciumPulse >= threshold_r2_calciumPulse ...
            %     && temp_calciumPulse_decayCoeff >= threshold_calciumPulse_decayCoeff
            %if (temp_r2_calciumPulse+temp_calciumPulse_decayCoeff) >= (threshold_r2_calciumPulse+threshold_calciumPulse_decayCoeff) ...
            %    && temp_calciumPulse_decayCoeff >= threshold_calciumPulse_decayCoeff
            if temp_r2_calciumPulse >= threshold_r2_calciumPulse ...
                && temp_calciumPulse_decayCoeff >= threshold_calciumPulse_decayCoeff
            
                if length(temp_negativeFvalleys) < length(temp_validCalciumPluse_frameIndex_valid)
                    xypix_ratio_badBoolIndex(tempi) = false;
                else
                    negativeFpeak_bool(tempi) = true;
                end
            end
        end
        
    end
    
    a = 1;
    temptempIndex_raw;
    temp_validCalciumPluse_frameIndex_valid;
    tempCalciumPulse_templength_valid;
    tempi;
        
    %if temp_validSpkFlag == false && cell_npix(tempi) >= 10
    %if temp_validSpkFlag == false && cell_npix(tempi) >= 25
    if temp_validSpkFlag == false && cell_npix(tempi) >= threshold_npix_smallCell
        tempBool = false;
        if length(temp_validCalciumPluse_frameIndex_valid) > threshold_validSpkCount_max
            tempBool = true;
        end
        if sum(temp_validCalciumPluse_frameIndex_valid<T/2) > floor(threshold_validSpkCount_max/2)
            tempBool = true;            
        end        
        if sum(temp_validCalciumPluse_frameIndex_valid>T/2) > floor(threshold_validSpkCount_max/2)
            tempBool = true;            
        end
        
        if tempBool == true
            manySpkCountCell_bool_all(tempi) = true;
            
            %if r2_calciumPulse(tempi) >= (threshold_r2_calciumPulse-compensateValue_forActiveCell) ...
            %        && calciumPulse_decayCoeff(tempi) >= (threshold_calciumPulse_decayCoeff-compensateValue_forActiveCell)
            %if (r2_calciumPulse(tempi)+calciumPulse_decayCoeff(tempi)) >= ...
            %        (threshold_r2_calciumPulse+threshold_calciumPulse_decayCoeff-compensateValue_forActiveCell*2) ...
            %       && calciumPulse_decayCoeff(tempi) >= (threshold_calciumPulse_decayCoeff-compensateValue_forActiveCell)
            %if (r2_calciumPulse(tempi)+calciumPulse_decayCoeff(tempi)) >= ...
            %        (threshold_r2_calciumPulse+threshold_calciumPulse_decayCoeff-compensateValue_forActiveCell*2) ...
            %        && r2_calciumPulse(tempi) >= (threshold_r2_calciumPulse-compensateValue_forActiveCell*2) ...
            %        && calciumPulse_decayCoeff(tempi) >= (threshold_calciumPulse_decayCoeff-compensateValue_forActiveCell)
            if (r2_calciumPulse(tempi)+calciumPulse_decayCoeff(tempi)) >= ...
                   (threshold_r2_calciumPulse+threshold_calciumPulse_decayCoeff-compensateValue_forActiveCell*2) ...
                   && r2_calciumPulse(tempi) >= (threshold_r2_calciumPulse-compensateValue_forActiveCell*2) ...
                   && calciumPulse_decayCoeff(tempi) >= (threshold_calciumPulse_decayCoeff-compensateValue_forActiveCell*1)
                
                temp_validSpkFlag = true;
                manySpkCountCell_bool(tempi) = true;
            end
        end
                
    end
    
    %if peak_tempID_F > threshold_bigF
    %    temp_validSpkFlag = true;
    %end
    
    if tempID_F_median > threshold_median
        temp_validSpkFlag = true;
        %bigCell_bool(tempi) = true; %#ok<*PFOUS>
        brightCell_bool(tempi) = true; %#ok<*PFOUS>        
    end
    
    if temp_validSpkFlag == false
        iscell(tempi, 2) = 0;
    elseif temp_validSpkFlag == true
        toal_validSpkCount_boolIndex(tempi) = true;
        if iscell(tempi,2) <= 0.01
            iscell(tempi,2) = 0.011;
        end
        
        if temp_stat.npix > threshold_npix_bigCell
            bigCell_bool(tempi) = true;
        end
    end
    
end

a = 1;
% sum(manySpkCountCell_bool);
manySpkCountCell_bool_raw = manySpkCountCell_bool;
% manySpkCountCell_bool = manySpkCountCell_bool_raw & (cell_npix >= 10);
manySpkCountCell_bool = manySpkCountCell_bool_raw;
manySpkCountCell_index = find(manySpkCountCell_bool==true);
manySpkCountCell_index_suite2p = manySpkCountCell_index-1;

fprintf('Convert %d/%d rois (notcell) with frequent spike.\n',length(manySpkCountCell_index),sum(manySpkCountCell_bool_all));

temp_iscellBoolIndex = iscell(:,2)>0.01;
temp_sum = sum(temp_iscellBoolIndex);

temp_rescueThreshold = ceil(size(iscell,1)/3);
temp_rescueThresholdBoolIndex = false(length(stat),1);
temp_rescueThresholdBoolIndex(1:temp_rescueThreshold) = true;

% peak_F_median_iscell = median(peak_F(temp_iscellBoolIndex));
% peak_F_prctile_iscell = prctile(peak_F(temp_iscellBoolIndex),85);
% tempBoolIndex = (~temp_iscellBoolIndex) & temp_rescueThresholdBoolIndex & (peak_F'>peak_F_prctile_iscell);

% peak_F_dff_prctile_iscell = prctile(peak_F_dff(temp_iscellBoolIndex),60);%85-->75-->83-->90-->85-->75-->85-->75-->70-->60
peak_F_dff_prctile_iscell = prctile(peak_F_dff(temp_iscellBoolIndex),65);%85-->75-->70-->60-->50-->75-->65
tempBoolIndex = (~temp_iscellBoolIndex) & temp_rescueThresholdBoolIndex & (peak_F_dff'>peak_F_dff_prctile_iscell);

% threshold_rescueNum = 100;
% temp_threshold_rescueNum = threshold_rescueNum-length(manySpkCountCell_index);
% temp_threshold_rescueNum = max(threshold_rescueNum-length(manySpkCountCell_index),round(threshold_rescueNum/2));
% temp_threshold_rescueNum = max(threshold_rescueNum-length(manySpkCountCell_index),round(threshold_rescueNum/3));
temp_threshold_rescueNum = max(threshold_rescueNum-length(manySpkCountCell_index),0);
if sum(tempBoolIndex) > temp_threshold_rescueNum
    tempFindIndex = find(tempBoolIndex==true);
    tempFindIndex = tempFindIndex(1:temp_threshold_rescueNum);
    if isempty(tempFindIndex) == true
       tempFindIndex = 0; 
    end
    tempBoolIndex( (tempFindIndex(end)+1):end ) = false;
end

fprintf('Rescue %d rois (notcell) with peak more than x percentile.\n',sum(tempBoolIndex));

iscell(tempBoolIndex,2) = 0.011;

%% Get highCorrCluster with high median roi
corr_threshold = 0.4;%0.4-->0.5-->0.4
corr_threshold2 = 0.7;%0.7

threshold_highCorrCluster_size = 10;%50-->20-->10
threshold2_highCorrCluster_size = 200;%40-->100

iscell_raw = iscell;

highCorrClusterNum = zeros(1,length(stat));
highCorrClusterIndex_cell = cell(1,length(stat));


%%
for tempi=1:length(stat)
    %if bigCell_bool(tempi) == false && brightCell_bool(tempi) == false
    if brightCell_bool(tempi) == false
        continue
    end
        
    tempID_corr_F= corr(F_2_bin_T(:,tempi),F_2_bin_T); %#ok<*PFBNS>
    tempID_corrBool_F = tempID_corr_F >= corr_threshold;
    
    highCorrClusterNum(tempi) = sum(tempID_corrBool_F);
    
    bigCell_bool;
    
    %if bigCell_bool(tempi) == true ...
    %        && sum(tempID_corrBool_F) > threshold_highCorrCluster_size
    if sum(tempID_corrBool_F) > threshold_highCorrCluster_size

        
        if sum(tempID_corrBool_F) < threshold2_highCorrCluster_size
            tempID_corrBool_F = false(1,length(stat));
            fprintf('Now clear the new highCorrCluster (size more than %d) with suite2p id %d, cluster size = %d!\n',...
                threshold_highCorrCluster_size,tempi-1,sum(tempID_corrBool_F));
        elseif sum(tempID_corrBool_F) >= threshold2_highCorrCluster_size
            
            temp_highCorrClusterIndex = find(tempID_corrBool_F==true);
            
            temp_F_highCorrCluster = F_2_bin_T(:,temp_highCorrClusterIndex); %#ok<*FNDSB>
            temp_Fmean_highCorrCluster = mean(temp_F_highCorrCluster,2);
            
            tempID_corr_F3= corr(temp_Fmean_highCorrCluster,F_2_bin_T); %#ok<*PFBNS>
            tempID_corrBool_F3 = abs(tempID_corr_F3) >= corr_threshold2;
            tempCellID_toBeClear = find(tempID_corrBool_F3==true);
            tempCellID_suite2p_toBeClear = tempCellID_toBeClear - 1;
            
            iscell(tempCellID_toBeClear,2) = 0;
            
            fprintf('Now clear the existing (bigcell) highCorrCluster (corr>%.2f) with suite2p id %d, cluster size = %d!\n',...
                corr_threshold2,tempi-1,length(tempCellID_toBeClear));
            
            tempID_corrBool_F = false(1,length(stat));
        end
    end
    
    highCorrClusterIndex = find(tempID_corrBool_F==true)';
    highCorrClusterIndex_suite2p = highCorrClusterIndex - 1;
    for templ=1:sum(tempID_corrBool_F)
        if iscell_raw(highCorrClusterIndex(templ),2) <= 0.01
            iscell(highCorrClusterIndex(templ),2) = 0.011;
            
            highCorrClusterIndex_cell{tempi} = [highCorrClusterIndex_cell{tempi} highCorrClusterIndex(templ)];
        end
    end
    if sum(tempID_corrBool_F) > 1
        fprintf('Get highCorrCluster (corr>%.2f) with suite2p id %d, cluster size = %d!\n',...
            corr_threshold,tempi-1,sum(tempID_corrBool_F));
    end
end    

a = 1;


xypix_ratio_badIndex = find(xypix_ratio_badBoolIndex==true);
iscell(xypix_ratio_badIndex, 2) = 0; %#ok<*FNDSB>
xypix_ratio_badIndex_suite2p = xypix_ratio_badIndex - 1;
fprintf('final xypix_ratio_badIndex_suite2p (num) = %d.\n',length(xypix_ratio_badIndex_suite2p));
xypix_ratio_badBoolIndex_final = xypix_ratio_badBoolIndex;
xypix_ratio_badIndex_final = xypix_ratio_badIndex;


%%

% sum(highFneu_bool);
highFneu_bool_raw = highFneu_bool;
highFneu_bool = highFneu_bool & (iscell(:,2)>0.01)';
highFneu_boolIndex = find(highFneu_bool==true);
highFneu_booIndex_suite2p = highFneu_boolIndex-1;

%% Others
toal_validSpkCount = sum(toal_validSpkCount_boolIndex);
fprintf('Identify %d ROIs with valid spike(s).\n', toal_validSpkCount);
iscell(:,1) = iscell(:,2)>0.01;
cell_count = sum(iscell(:,1));
%cell_count = sum(iscell(:,2)>0.01);
fprintf('Identify %d valid ROIs (including soma and dendrite), time = %.0f secs.\n',cell_count,toc(t0));

writeNPY(iscell, fullFileName_iscell);
if if_profile == 1
    profile viewer
end
a = 1;


%% End
