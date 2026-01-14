function A = read_raw_file(filename,framenum,window,FOV,bitsize)

%% Starting from frame number: 'framenum', read 'window'
% fid = fopen(filename);
fid = fopen(filename, 'r');
imsize = FOV(1)*FOV(2)*bitsize;                                                   % Bit size of single frame

current_seek = ftell(fid);
fseek(fid, 0, 1);
file_length = ftell(fid);
fseek(fid, current_seek, -1);
frame_length = file_length/imsize;

if bitsize == 2
    bstring = 'uint16=>single';
else
    bstring = 'real*4=>double';
end
% fclose(fid);

%%
% window = min(window,frame_length);
A = zeros(FOV(1),FOV(2),window,'single');

% wokerNum = 1;
% N = window;
% 
% i=1:N;
% i_worker = floor(N/wokerNum) * ones(1, wokerNum);
% i_worker(1:N - sum(i_worker)) = i_worker(1:N - sum(i_worker)) + 1;
% 
% i_worker_start = zeros(1, wokerNum);
% i_worker_end = zeros(1, wokerNum);
% for tempj=1:wokerNum
%     if tempj == 1
%         i_worker_start(1) = 1;
%     else
%         i_worker_start(tempj) = sum(i_worker(1:tempj-1))+1;
%     end
%     i_worker_end(tempj) = i_worker_start(tempj)+i_worker(tempj)-1;
% end
% 
% t0 = tic;
% spmd(0, wokerNum)
%     fids = fopen(filename,'r');
% 
%     %fid;
%     %fprintf('fid %d\n', fids);
%     %fprintf('fid %d\n', fid);
%     
%     %labindex;
%     %numlabs;
%     tempi = i_worker_start(labindex):i_worker_end(labindex);
%     tempi_length = length(tempi);
%     for tempj=1:tempi_length
%         n = tempi(tempj);
%         w = n;
%         t1 = tic;
%                 
%         %fprintf('n=%d\n',n);
%         
%         %fseek(fids,(framenum-1+w)*imsize,'bof');                                 % Position the file-indicator at the beginning of the frame
%         temp = fread(fids,[FOV(1),FOV(2)],bstring,0,'b');                              % Read the frame
%         A(:,:,w+1) = temp';         
%         %fprintf('a=%d\n',1);
%         %size(temp)
%         %fprintf('size = %d\n',size(temp,1));
%         
%             
%         fprintf(['--------------------\n', ...
%             'time = %.1f, Total time = %.1f, %.1f%% \n'],...
%             toc(t1), toc(t0), 100*tempj/tempi_length);
%     end 
%     fclose(fids);
% end
% fprintf('--------------------\n');
% fprintf('All is over! Total time = %.2f secs.\n', toc(t0));
% 
% a = 1;


for w = 0:window-1                                                          % For each frame inside the window
% parfor w = 0:window-1
    fseek(fid,(framenum-1+w)*imsize,'bof');                                 % Position the file-indicator at the beginning of the frame
    temp = fread(fid,[FOV(1),FOV(2)],bstring,0,'b');                              % Read the frame
    A(:,:,w+1) = temp';                                                                    
end
fclose(fid);