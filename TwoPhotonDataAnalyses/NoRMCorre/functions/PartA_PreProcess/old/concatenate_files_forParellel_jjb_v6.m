function [Ycon,ln,error_flag] = concatenate_files_forParellel_jjb_v6(files,filename,output_type)

% concatenate a list of files into a single file
% INPUTS
% files:        list of a files as a struct array. Name of each file is in files(i).name
% filename:     output filename
% type:         output type ('mat' array loaded in RAM (default), 'hdf5' and 'tif' are currently supported)

% OUTPUTS
% Ycon:         concantenated file as an array if output_type == 'mat'
% ln:           length of each file (in frames)

error_flag = 0;

if exist(filename , 'file') == 2
    if strcmp(output_type,'hdf5') == 1
        delete(filename);       
    else
        error_flag = 1;
        Ycon = '';
        ln = 0;
        fprintf('Target file has already existed.\n');
        return
    end
end

numFiles = length(files);
ln = zeros(numFiles,1);

% Y1 = read_file(files(1).name);
% sizY = size(Y1);
% T1 = sizY(end);
% data_type = class(Y1);
% nd = max(ndims(Y1)-1,2);
% sizY = sizY(1:nd);


info = h5info(files(1).name);
dims = info.Datasets.Dataspace.Size;
sframe = 1;
num2read = dims(end)-sframe+1;
if num2read == 1
    error_flag = 1;
    Ycon = '';
    ln = 0;
    fprintf('Too low pictures for Bin!\n');
    return
end
sizY = dims(1:2);
data_type = 'uint16';


Ycon = filename;
h5create(filename,'/mov',[sizY,Inf],'Chunksize',[sizY,1],'Datatype',data_type);
%h5create(filename,'/mov',[sizY,Inf],'Chunksize',[sizY,1],'Datatype',data_type,'Deflate',1,'Shuffle',true);


for i = 1:numFiles
    info = h5info(files(i).name);
    dims = info.Datasets.Dataspace.Size;
    sframe = 1;
    num2read = dims(end)-sframe+1;
    ln(i) = num2read; %#ok<*SAGROW>
    
    if num2read == 1
        error_flag = 1;
        Ycon = '';
        ln = 0;
        fprintf('Too low pictures for Bin!\n');
        return
    end    
    if num2read == 0
        %error_flag = 1;
        %Ycon = '';
        %ln = 0;
        %fprintf('Too low pictures!\n');
        %return
        continue
    end     
end

% Y_temp_1 = h5read(files(1).name,'/mov',[1,1,1],[512,512,ln(1)]);
% numout = 1;
% for i = 1:numFiles       
%     %Y_temp = h5read(files(i).name,'/mov',[1,1,1],[512,512,ln(i)]);
%     
%     if i < numFiles
%         F = parfeval(@h5read,numout,...
%             files(i+1).name,'/mov',[1,1,1],[512,512,ln(i+1)]);
%     end
%     
%     if i == 1
%         Y_temp = Y_temp_1;
%         Y_temp_1 = [];
%     else
%         Y_temp = fetchOutputs(last_F);
%     end
%     
%     h5write(filename,'/mov',Y_temp,[1,1,sum(ln(1:i-1))+1],[512,512,ln(i)]);
%     last_F = F;
% end

%fprintf('Now is concatenating! (Very time-consuming).\n');

for i = 1:numFiles       
    Y_temp = h5read(files(i).name,'/mov',[1,1,1],[512,512,ln(i)]);
    h5write(filename,'/mov',Y_temp,[1,1,sum(ln(1:i-1))+1],[512,512,ln(i)]);
end

% for i = 1:numFiles
% %     if i == 1        
% %         info = h5info(files(i).name);
% %         dims = info.Datasets.Dataspace.Size;
% %         sframe = 1;
% %         num2read = dims(end)-sframe+1;    
% %     end    
% %     Y_temp = read_file(files(i).name,sframe+num2read*(i-1),num2read);
% 
% 
%     info = h5info(files(i).name);
%     dims = info.Datasets.Dataspace.Size;
%     sframe = 1;
%     num2read = dims(end)-sframe+1;
%     if num2read == 1
%         error_flag = 1;
%         Ycon = '';
%         ln = 0;
%         fprintf('Too low pictures for Bin!\n');
%         return
%     end
%     
%     if num2read == 0
%         %error_flag = 1;
%         %Ycon = '';
%         %ln = 0;
%         %fprintf('Too low pictures!\n');
%         %return
%         continue
%     end    
%         
%     Y_temp = read_file(files(i).name);    
%     ln(i) = size(Y_temp,ndims(Y_temp));   
%     
%     size_Y_temp = size(Y_temp);
%     
%     switch lower(output_type)
%         case 'mat'
%             Ycon = cat(nd+1,Ycon,Y_temp); 
%         case {'hdf5','h5'}
%             h5write(filename,'/mov',Y_temp,[ones(1,nd),sum(ln(1:i-1))+1],size_Y_temp);
%         case {'tif','tiff'}
%             saveastiff(cast(Y_temp,data_type),filename,opts_tiff);
%     end         
% end