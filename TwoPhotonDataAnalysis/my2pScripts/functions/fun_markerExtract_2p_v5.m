function output = fun_markerExtract_2p_v5(fullFileName)

marker_2p_fullName = fullFileName;
[~,temp_fname,temp_ext] = fileparts(marker_2p_fullName);
marker_2p_fileName = [temp_fname,temp_ext];

%% Load marker file
fid_2p = fopen(marker_2p_fullName,'r');
current_seek = ftell(fid_2p);
fseek(fid_2p,0,'eof');
file_length_2p_raw = ftell(fid_2p);
fseek(fid_2p,current_seek,'bof');
[marker_data_2p, file_length_2p] = fread(fid_2p,Inf,'*uint8');
fprintf("Load %s ....\n",marker_2p_fileName);
% fprintf("------  Load %s ....\n",marker_2p_fileName);
if file_length_2p == file_length_2p_raw
    clear file_length_2p_raw
else
    fprintf('file_length_2p_raw is different from file_length_2p!\n');
end
fclose(fid_2p);

%% Parse marker file
low_count_continous_all = zeros(1,file_length_2p*8);
low_count_continous_all_sampleIndex = zeros(1,file_length_2p*8);
low_count_continous_2p = zeros(1,file_length_2p*8);
low_count_continous = 0;
low_count = 0;
twoP_marker_count = 0;
twoP_marker_array = false(1,file_length_2p*8);
twoP_marker_index = zeros(1,file_length_2p*8);
sample_count = 0;
last_low_count_continous = 0;
validRecordingFlag = 0; %#ok<*NASGU>
last_validRecordingFlag = 0;

validRecordingBlock_count = 0;
validRecordingBlock_start = 0;
validRecordingBlock_end = 0;

temp_count = 0;

for tempi=1:file_length_2p
    for tempj=1:8
        sample_count = sample_count + 1;
        if bitget(marker_data_2p(tempi),tempj) == 0
            low_count_continous = low_count_continous + 1;
            low_count = low_count + 1;
        else
            if low_count_continous > 0
                %fprintf("There are %d low samples! (sample_count %d to sample_count %d).\n",...
                %  low_count_continous, sample_count-low_count_continous, sample_count-1);
                
                temp_count = temp_count + 1;
                low_count_continous_all(temp_count) = low_count_continous;
                low_count_continous_all_sampleIndex(temp_count) = sample_count-low_count_continous+1;
                
                current_sampleRange = sample_count-low_count_continous+1:sample_count;
                
                if low_count_continous > 2 && low_count_continous <= 10
                %if low_count_continous <= 10
                    validRecordingFlag = 1;
                    twoP_marker_count = twoP_marker_count + 1;                                        
                    twoP_marker_array(current_sampleRange) = true;
                    twoP_marker_index(twoP_marker_count) = current_sampleRange(1);
                    
                    low_count_continous_2p(twoP_marker_count) = low_count_continous;
                else
                    validRecordingFlag = 0;
                end
                                
                if last_validRecordingFlag == 0 && validRecordingFlag == 1
                    % When come into a new valid recording block
                    validRecordingBlock_count = validRecordingBlock_count + 1;
                    validRecordingBlock_start(validRecordingBlock_count) = current_sampleRange(1);%#ok<*AGROW>
                    a = 1;
                elseif last_validRecordingFlag == 1 && validRecordingFlag == 0
                    % When finish a valid recording block
                    validRecordingBlock_end(validRecordingBlock_count) = current_sampleRange(1);                                     
                end
                
                last_validRecordingFlag = validRecordingFlag;                
                low_count_continous = 0;
            end
        end        
    end
end
twoP_marker_index(twoP_marker_index==0) = [];
low_count_continous_all(low_count_continous_all==0) = [];
low_count_continous_all_sampleIndex(low_count_continous_all_sampleIndex==0) = [];
low_count_continous_2p(low_count_continous_2p==0) = [];

if sum(validRecordingBlock_end) == 0
    validRecordingBlock_end = sample_count;
end


a = 1;

% validRecordingBlock_start(validRecordingBlock_start < validRecordingBlock_end(1)) = [];

if length(validRecordingBlock_start) > length(validRecordingBlock_end)
    validRecordingBlock_end = [validRecordingBlock_end sample_count];
end


validRecordingBlock_length = validRecordingBlock_end - validRecordingBlock_start;
[~,I] = max(validRecordingBlock_length);
validRecordingBlock_range = [validRecordingBlock_start(I),validRecordingBlock_end(I)];

% fprintf("validRecordingBlock_range = %d --> %d (samples).\n",...
%     validRecordingBlock_range(1),validRecordingBlock_range(2));

valid_frames_boolIndex = (twoP_marker_index >= validRecordingBlock_range(1)) & ...
                (twoP_marker_index <= validRecordingBlock_range(2));
valid_frames_sampleIndex = twoP_marker_index(valid_frames_boolIndex);

validRecordingBlock_range(2) = valid_frames_sampleIndex(end);

fprintf("validRecordingBlock_range = %d --> %d (samples).\n",...
    validRecordingBlock_range(1),validRecordingBlock_range(2));
fprintf("Get %d valid frames from total %d 2P marker.\n",sum(valid_frames_boolIndex),twoP_marker_count);


valid_frames_sampleInterval = diff(valid_frames_sampleIndex);
[min_sampleInterval,max_sampleInterval] = bounds(valid_frames_sampleInterval);
median_sampleInterval = median(valid_frames_sampleInterval);
mean_sampleInterval = mean(valid_frames_sampleInterval);
min_timeInterval = min_sampleInterval/4000*1000;
max_timeInterval = max_sampleInterval/4000*1000;
median_timeInterval = median_sampleInterval/4000*1000;
mean_timeInterval = mean_sampleInterval/4000*1000;
% fprintf('min_sampleInterval_2p=%d, max_sampleInterval_2p=%d.\n',min_sampleInterval,max_sampleInterval);


output = struct;
output.file_length = file_length_2p;
output.low_count = low_count;
output.marker_data = marker_data_2p;
output.marker_array = twoP_marker_array;
output.marker_count = twoP_marker_count;
output.marker_index = twoP_marker_index;
output.valid_frames_boolIndex = valid_frames_boolIndex;
output.valid_frames_sampleIndex = valid_frames_sampleIndex;
output.valid_frames_sampleInterval = valid_frames_sampleInterval;
output.validRecordingBlock_count = validRecordingBlock_count;
output.validRecordingBlock_range = validRecordingBlock_range;
output.min_timeInterval = min_timeInterval;
output.max_timeInterval = max_timeInterval;
output.median_timeInterval = median_timeInterval;
output.mean_timeInterval = mean_timeInterval;

%% End
