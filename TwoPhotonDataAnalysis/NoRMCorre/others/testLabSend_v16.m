%% Initialization
clear
close all

gpuDevice(1);
gcp;% To open parallel pool

if false
    delete(gcp); %#ok<*UNRCH>
    gcp;
end

if_profileOn = 0;


if if_profileOn == 1
    profile on
end

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
cd(targetPATH)

dataRoot_path = 'C:\ASDROOT\STUDY\twoPhotonData';
currentSession = '113Recording_Bubble_20221211A_fixZtest';
currentSession_path = [dataRoot_path '\' currentSession];

h5_filename = 'testLabSendData.h5';
raw_fileName = 'Image_001_001.raw';

raw_fileFullName = [currentSession_path '\' raw_fileName];
h5_fileFullName = [currentSession_path '\' h5_filename];


workerNum_reader = 2;%3
workerNum_readBuffer = 1;
workerNum_processor = 1;%1
workerNum_writer = 1;%1
numFrames = 6010;%1000


bin_width = 25;%25

%workerNum = 9;%9
spmdCut_reader = workerNum_reader;
spmdCut_readBuffer = spmdCut_reader + workerNum_readBuffer;
spmdCut_processor = spmdCut_readBuffer + workerNum_processor;
spmdCut_writer = spmdCut_processor + workerNum_writer;

workerIndex_reader = 1:spmdCut_reader;
workerIndex_readBuffer = (spmdCut_reader+1):spmdCut_readBuffer;
workerIndex_processor = (spmdCut_readBuffer+1):spmdCut_processor;
workerIndex_writer = (spmdCut_processor+1):spmdCut_writer;


workerNum = workerNum_reader + workerNum_readBuffer + workerNum_writer + workerNum_processor;
fprintf('Worker number = %d.\n',workerNum);


fid = fopen(raw_fileFullName);
imsize = 512*512*2;
current_seek = ftell(fid);
fseek(fid, 0, 1);
file_length = ftell(fid);
fseek(fid, current_seek, -1);
T = file_length/imsize;
if numFrames ~= -1
    T = min(T, numFrames);
end
T = round(T);
sizY = [512,512,T];
fclose(fid);

fprintf('Load %s.\n', raw_fileName);
T = sizY(end);
fprintf('Load %d frames.\n', T);

bin_totalNum = length(1:bin_width:T);
N = bin_totalNum;

i_worker = floor(N/spmdCut_reader) * ones(1, spmdCut_reader);
i_worker(1:N - sum(i_worker)) = i_worker(1:N - sum(i_worker)) + 1;
i_worker_start = zeros(1, spmdCut_reader);
i_worker_end = zeros(1, spmdCut_reader);
for tempj=1:spmdCut_reader
    if tempj == 1
        i_worker_start(1) = 1;
    else
        i_worker_start(tempj) = sum(i_worker(1:tempj-1))+1;
    end
    i_worker_end(tempj) = i_worker_start(tempj)+i_worker(tempj)-1;
end

[pathstr,fname,ext] = fileparts(h5_fileFullName);
h5_fileFullName = fullfile(pathstr,[fname,'_',datestr(now,30),ext]); %#ok<*TNOW1,*DATST> 

spmd(0, workerNum)
   a = 1; 
end
t0 = tic;
spmd(0, workerNum)
    %spmdindex = labindex;    
    %% Reader
    if ismember(spmdIndex,workerIndex_reader)
        fid = fopen(raw_fileFullName);
        imsize = 512*512*2;
        
        tempi_worker = i_worker_start(spmdIndex):i_worker_end(spmdIndex);
        tempi_length = length(tempi_worker);
        
        time_spmd = 0;
        for tempj=1:tempi_length
            t0_spmd = tic;
            
            n = tempi_worker(tempj);
            bin_count = n;
            t = (bin_count-1)*bin_width + 1;
            
            % Read Raw file
            window = min(t+bin_width-1,T)-t+1;
            feek_indicator = (t-1)*imsize;
            fseek(fid,feek_indicator,'bof');
            Y_raw = fread(fid,512*512*window,'uint16=>single',0,'l')';
            
            %time_spmd = time_spmd + toc(t0_spmd);
            
            tag = t;            
            spmdSend(Y_raw,workerIndex_readBuffer,tag);% Wait too loooong for receiving!
            Y_raw = [];
            
            time_spmd = time_spmd + toc(t0_spmd);
        end             
        
        fprintf('Read over!\n');
        fclose(fid); 
        
    %% Read Buffer        
    elseif ismember(spmdIndex,workerIndex_readBuffer)
        if_readBuffOver = 0;
        readOver_markerNum = 0;
        
        totalSend_binNum = 0;
        buffer_binNum = 0;
        buffer = struct;
        buffer.Y = [];
        buffer.t = [];
        buffer.count = [];
        processReady = 0;
        
        %a = 2;
        
        time_spmd = 0;
        while if_readBuffOver == 0
            isDataAvail_reader = false(1, workerNum_reader);
            for tempi=1:workerNum_reader
               isDataAvail_reader(tempi) = spmdProbe(workerIndex_reader(tempi));                
            end
            
            isDataAvail_processor = spmdProbe(workerIndex_processor);            
            
            while sum(isDataAvail_reader) > 0                
            %if sum(isDataAvail_reader) > 0 %&& processReady == 1
                t0_spmd = tic;
                
                temp_index = find(isDataAvail_reader == 1,1,'first');
                %temp_index = workerIndex_reader(find(isDataAvail_reader == 1,1,'last'));
                [data,~,tag] = spmdReceive(workerIndex_reader(temp_index));
                isDataAvail_reader(temp_index) = 0;
                
                %time_spmd = time_spmd + toc(t0_spmd);
                %t0_spmd = tic;
                
                buffer_binNum = buffer_binNum + 1;
                buffer.Y = [buffer.Y data];
                buffer.t = [buffer.t tag];
                buffer.count = [buffer.count length(data)/512/512];
                %fprintf('Receive data from reader!\n');
                time_spmd = time_spmd + toc(t0_spmd);
                
%                 tag = 2;
%                 labSend(buffer,workerIndex_processor,tag);% Wait too loooong for receiving!
%                 totalSend_binNum = totalSend_binNum + buffer_binNum;
%                 buffer_binNum = 0;
%                 buffer.Y = [];
%                 buffer.t = [];
%                 buffer.count = [];
%                 %fprintf('Send to processor!\n');
%                 
%                 if totalSend_binNum == bin_totalNum
%                     if_readBuffOver = 1;
%                     %labSend('readBuffOver',workerIndex_processor);
%                     fprintf('Read Buffer over!\n');
%                 end
                
                %processReady = 0;
                %time_spmd = time_spmd + toc(t0_spmd);                
            end
            
            if isDataAvail_processor == 1                                                
                data = spmdReceive(workerIndex_processor);                                
                if strcmp(data,'processReady') == 1
                    %processReady = 1;
                    t0_spmd = tic;
                    
                    tag = 2;
                    spmdSend(buffer,workerIndex_processor,tag);% Wait too loooong for receiving!
                    
                    
                    time_spmd = time_spmd + toc(t0_spmd);
                    
                    totalSend_binNum = totalSend_binNum + buffer_binNum;
                    buffer_binNum = 0;
                    buffer.Y = [];
                    buffer.t = [];
                    buffer.count = [];
                    %fprintf('Send to processor!\n');
                    
                    if totalSend_binNum == bin_totalNum
                        if_readBuffOver = 1;
                        %labSend('readBuffOver',workerIndex_processor);
                        fprintf('Read Buffer over!\n');
                    end
                    %time_spmd = time_spmd + toc(t0_spmd);
                end
                
            end

        end        
        
    %% Processor    
    elseif ismember(spmdIndex,workerIndex_processor)
        if_processOver = 0;
        if_readBuffOver = 0;
        processor_binNum = 0;
        writeReady = 1;
                
        
        t0_processReady = tic;        
        time_spmd = 0;
        while if_processOver == 0
            isDataAvail_writer = spmdProbe(workerIndex_writer);
            isDataAvail_readBuffer = spmdProbe(workerIndex_readBuffer,2);            
            
            if isDataAvail_writer == 1
                %data = labReceive(workerIndex_writer);
                %if strcmp(data,'writeReady') == 1
                %    writeReady = 1;
                %end     
                writeReady = 1;
            end
            
            if isDataAvail_readBuffer == 1
                %a
                
                %[buffer,~,tag] = labReceive(workerIndex_readBuffer);
                %[buffer,~,tag] = labReceive;
                
                %if tag == 2 && writeReady == 1
                if writeReady == 1
                    t0_spmd = tic;
                    
                    buffer = spmdReceive(workerIndex_readBuffer,2);
                    %time_spmd = time_spmd + toc(t0_spmd);
                    
                    %Y_raw = gpuArray(buffer.Y);
                    %window = length(Y_raw)/512/512;
                    %Y_raw = reshape(Y_raw,[512 512 window]);
                    %Y_raw = pagetranspose(Y_raw);
                    %asd = gather(uint16(Y_raw));
                    
                    Y_raw = gpuArray(buffer.Y);
                    window = length(Y_raw)/512/512;
                    Y_raw = reshape(Y_raw,[512 512 window]);
                    Y_raw = pagetranspose(Y_raw);
                    buffer.Y = gather(uint16(Y_raw));
                    
                    %buffer.count = window;
                    
                    %t = buffer.t;
                    %if length(t) == 1
                    %    buffer.count = window;
                    %else
                    %    buffer.count = zeros(1, length(t));
                    %    for tempi=1:length(t)
                    %        if tempi < length(t)
                    %            buffer.count(tempi) = t(tempi+1) - t(tempi);
                    %        elseif tempi == length(t)
                    %            buffer.count(tempi) = window - sum(buffer.count(1:tempi-1));
                    %        end
                    %    end
                    %end
                    
                    %t0_spmd = tic;
                    tag = 3;                    
                    %time_spmd = time_spmd + toc(t0_spmd);
                    spmdSend(buffer,workerIndex_writer,tag);
                    time_spmd = time_spmd + toc(t0_spmd);
                    processor_binNum = processor_binNum + ceil(window/bin_width);
                    buffer = [];
                    Y_raw = [];
                    if_Y_raw_empty = 1;
                    
                    %fprintf('Process data!\n');
                    
                    if processor_binNum == bin_totalNum
                        if_processOver = 1;
                        fprintf('Process over!\n');
                    end
                    %writeReady = 0;
                    %time_spmd = time_spmd + toc(t0_spmd);
                end                
            elseif isDataAvail_readBuffer == 0
                if if_processOver == 0
                    if toc(t0_processReady) > 10/1000
                      spmdSend('processReady',workerIndex_readBuffer);
                      t0_processReady = tic;
                    end
                end
            end
                        
            
            %elseif isDataAvail == 0
            %    %if if_processOver == 0 && if_readBuffOver == 0                                   
            %    if if_processOver == 0
            %        %if toc(t0_processReady) > 100/1000
            %        %    labSend('processReady',workerIndex_readBuffer);
            %        %    t0_processReady = tic;
            %        %end
            %    end
            %end
            
        end
        
    
    %% Writer
    elseif ismember(spmdIndex,workerIndex_writer)        
        if_exit = 0;
        writer_frameNum = 0;
        
        h5create(h5_fileFullName,'/mov',[512,512,Inf],'Chunksize',[512,512,bin_width],'Datatype','uint16');
        
        t0_writeReady = tic;
        time_spmd = 0;
        while if_exit == 0
            isDataAvail = spmdProbe(workerIndex_processor,3);
            if isDataAvail == 1
                [data,~,tag] = spmdReceive(workerIndex_processor,3);
                
                t0_spmd = tic;
                
                Y = data.Y;
                t = data.t;
                count = data.count;
                data = [];
                
                for tempi=1:length(t)
                    if tempi==1
                        temp_range = 1:count(tempi);
                    else
                        temp_range =  sum(count(1:tempi-1))+(1:count(tempi));
                    end
                    
                    %temp_range = bin_width*(tempi-1) + (1:bin_width);
                    
                    h5write(h5_fileFullName,'/mov',...
                      Y(:,:,temp_range),...
                      [1,1,t(tempi)],...
                      [512,512,count(tempi)]);
                    
                    writer_frameNum = writer_frameNum + count(tempi);
                end
                %fprintf('Write data!\n');
                Y = [];
                t = [];
                count = [];
                
                time_spmd = time_spmd + toc(t0_spmd);
                
                if writer_frameNum == T
                    if_exit = 1;
                    fprintf('Write over, %d frames!\n',writer_frameNum);
                end
            elseif isDataAvail == 0
                if toc(t0_writeReady) > 10/1000 && if_exit == 0
                    %labSend('writeReady',workerIndex_processor);
                    t0_writeReady = tic;
                end
            end
                       
        end  
        
    end
end
toc(t0);
fprintf('All over, exit!\n');
% size_data_all = size_data{spmdCut_writer};
% fprintf('Write %d frames.\n',sum(size_data_all));

time_spmd_all = [];
for tempi=1:workerNum
    time_spmd_all = [time_spmd_all time_spmd{tempi}]; %#ok<*AGROW>
end
            
fprintf("  time_spmd=");
for tempi=1:workerNum
    fprintf('%.1f ',time_spmd_all(tempi));
end
fprintf("\n");

if if_profileOn == 1
    profile viewer
end