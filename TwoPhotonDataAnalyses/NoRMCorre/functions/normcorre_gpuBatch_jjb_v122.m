function [M_final,shifts_g,template,options,col_shift,cY,cM] = ...
    normcorre_gpuBatch_jjb_v122(Y,options,template)

% online motion correction through DFT subpixel registration
% Based on the dftregistration.m function from Manuel Guizar and Jim Fienup

% INPUTS
% Y:                Input data, can be already loaded in memory as a 3D
%                   tensor, a memory mapped file, or a pointer to a tiff stack
% options:          options structure for motion correction (optional, rigid registration is performed if not provided)
% template:         provide template (optional)

% OUTPUTS
% M_final:          motion corrected data
% shifts_up:        upsampled shifts
% shifts:           originally calculated shifts
% template:         calculated template
% profile on
%% first determine filetype
if isa(Y,'char')
    [~,~,ext] = fileparts(Y);
    ext = ext(2:end);
    if strcmpi(ext,'tif') || strcmpi(ext,'tiff'); %#ok<*NOSEL>
        tiffInfo = imfinfo(Y);
        filetype = 'tif';
        T = length(tiffInfo);
        sizY = [tiffInfo(1).Height,tiffInfo(1).Width,T];
    elseif strcmpi(ext,'mat')
        filetype = 'mem';
        Y = matfile(Y,'Writable',true);
        details = whos(Y);
        var_sizes = [details.bytes];
        [~,var_ind] = max(var_sizes);
        var_name = details(var_ind).name;
        sizY = size(Y,var_name);
        T = sizY(end);
    elseif strcmpi(ext,'hdf5') || strcmpi(ext,'h5');
        filetype = 'hdf5';
        fileinfo = hdf5info(Y); %#ok<*HDFI>
        data_name = fileinfo.GroupHierarchy.Datasets.Name;
        sizY = fileinfo.GroupHierarchy.Datasets.Dims;
        T = sizY(end);
    elseif strcmpi(ext,'raw')
        filetype = 'raw';
        fid = fopen(Y);
        FOV = [options.d1,options.d2];
        bitsize = options.bitsize;
        imsize = FOV(1)*FOV(2)*bitsize;                                                   % Bit size of single frame
        current_seek = ftell(fid);
        fseek(fid, 0, 1);
        file_length = ftell(fid);
        fseek(fid, current_seek, -1);
        T = file_length/imsize;
        if options.numFrames ~= -1
            T = min(T, options.numFrames);
        end
        %T = int32(T);
        T = round(T);
        sizY = [FOV,T];
        fclose(fid);        
    end    
elseif isobject(Y)
    filetype = 'mem';
    var_name = 'Y';
    sizY = size(Y,var_name);
    T = sizY(end);
else % array loaded in memory
    filetype = 'mat';
    %Y = double(Y);
    sizY = size(Y);
    T = sizY(end);
end

if options.if_metrics == 1
    if options.if_imageAverage == 1
        cY = zeros(1, ceil(T/options.imageAverageBin));
        cM = zeros(1, ceil(T/options.imageAverageBin));
    else
        cY = zeros(1, T);
        cM = zeros(1, T);
    end
end

options.numFrames = T;
nd = length(sizY)-1;                          % determine whether imaging is 2d or 3d
sizY = sizY(1:nd);
%% set default parameters if not present

if ~exist('options','var') || isempty(options)
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2));
    if nd > 2; options.d3 = sizY(3); end
end

memmap = options.memmap;
grid_size = options.grid_size; 
mot_uf = options.mot_uf;
min_patch_size = options.min_patch_size;
overlap_pre = options.overlap_pre;
overlap_post = options.overlap_post;
upd_template = options.upd_template;
bin_width = options.bin_width;
buffer_width = options.buffer_width;
max_dev_g = options.max_dev;
init_batch = options.init_batch;
us_fac = single(options.us_fac);
method = options.method;
filename = options.mem_filename;
iter = options.iter;
add_value = options.add_value;
max_shift = options.max_shift;
if strcmpi(options.boundary,'nan')
    fill_value = NaN;
else
    fill_value = add_value;
end

while mod(T,bin_width) == 1
    if T == 1
        error('Movie appears to have only one frame. Use the function normcorre instead')        
    end
    bin_width = bin_width + 1;
end

%% first check for offset due to bi-directional scanning

if options.correct_bidir && isempty(options.col_shift)
    col_shift = correct_bidirectional_offset(Y,options.nFrames,options.bidir_us);
elseif ~isempty(options.col_shift)
    col_shift = options.col_shift;
else
    col_shift = 0;
end 
options.col_shift = col_shift;
if col_shift 
    fprintf('Offset %1.1d pixels due to bidirectional scanning detected. \n',col_shift); 
    if strcmpi(options.shifts_method,'fft')
        options.shifts_method = 'cubic';
        fprintf('Cubic shifts will be applied. \n'); 
    end
end


%% setup grids for patches

d1 = options.d1;
d2 = options.d2;
d3 = options.d3;

[xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf] = construct_grid(grid_size,mot_uf,d1,d2,d3,min_patch_size);
shifts_g = struct('shifts',cell(T,1),'shifts_up',cell(T,1),'diff',cell(T,1));
% temp_cell = mat2cell_ov(template_in,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap_post,sizY);
temp_cell = mat2cell_ov(magic(d1),xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap_post,sizY);

%% precompute some quantities that are used repetitively for template matching and applying shifts
Nr = cell(size(temp_cell));
Nc = cell(size(temp_cell));
Np = cell(size(temp_cell));
Bs = cell(size(temp_cell));
for i = 1:length(xx_us)
    for j = 1:length(yy_us)
        for k = 1:length(zz_us)
            [nr,nc,np] = size(temp_cell{i,j,k});
            nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
            nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
            np = ifftshift(-fix(np/2):ceil(np/2)-1);
            [Nc{i,j,k},Nr{i,j,k},Np{i,j,k}] = meshgrid(nc,nr,np);
            extended_grid = [max(xx_us(i)-overlap_post(1),1),min(xx_uf(i)+overlap_post(1),d1),max(yy_us(j)-overlap_post(2),1),min(yy_uf(j)+overlap_post(2),d2),max(zz_us(k)-overlap_post(3),1),min(zz_uf(k)+overlap_post(3),d3)];            
            Bs{i,j,k} = permute(construct_weights([xx_us(i),xx_uf(i),yy_us(j),yy_uf(j),zz_us(k),zz_uf(k)],extended_grid),[2,1,3]); 
            Bs{i,j,k} = single(Bs{i,j,k});
        end
    end
end
if nd == 2; Np = cellfun(@(x) 0,Nr,'un',0); end

%% read initial batch and compute template
init_batch = min(T,init_batch);

% init_Y_index = sort(randperm(T-1, init_batch));

% if options.if_imageAverage == 1
%     %temp_init_batch = init_batch * options.imageAverageBin;
%     temp_init_Y_index = sort([init_Y_index, init_Y_index+1]);
% else
%     temp_init_Y_index = init_Y_index;
% end

if options.if_imageAverage == 1
    temp_T = ceil(T/options.imageAverageBin);
    
    % random
    %init_Y_index = sort(randperm(temp_T-1, init_batch));
    
    % interleaved
    %init_Y_index = 1:floor((temp_T-1)/init_batch):(temp_T-1);
    %init_Y_index = init_Y_index(1:init_batch);    
    %init_Y_index = sort([init_Y_index*2,init_Y_index*2+1]);
    %temp_init_batch = length(init_Y_index);
    
    % first
    temp_init_batch = min(T,init_batch*2);
    init_Y_index = 1:temp_init_batch;  
else
    % random
    %init_Y_index = sort(randperm(T, init_batch));
    
    % interleaved
    %init_Y_index = 1:floor(T/init_batch):T;
    %init_Y_index = init_Y_index(1:init_batch);
    %temp_init_batch = init_batch;
    
    % first
    %init_Y_index = 1:init_batch;
    %temp_init_batch = init_batch;
    
    % middle
    init_Y_index = 1:init_batch;
    temp_init_batch = init_batch;
    
    T;
    temp_offset = floor((T-init_batch)/2);
    init_Y_index = init_Y_index + temp_offset;
    
end

Y_temp = zeros(options.d1,options.d2,init_batch,'single');
interval = ceil(T/2-init_batch/2+1):floor(T/2+init_batch/2);

t0_template = tic;
switch filetype
    case 'raw'
        fid = fopen(Y);
        imsize = 512*512*2;
        current_seek = ftell(fid);
        fseek(fid, 0, 1);
        file_length = ftell(fid);
        fseek(fid, current_seek, -1);
        frame_length = file_length/imsize;
        bstring = 'uint16=>single';
        for tempi=1:temp_init_batch
            t = init_Y_index(tempi);
            feek_indicator = (t-1)*imsize;
            fseek(fid,feek_indicator,'bof');
            %Y_temp(:,:,tempi) = fread(fid,[512,512],'uint16=>single',0,'l')';
            Y_temp(:,:,tempi) = (fread(fid,[512,512],'uint16=>single',0,'l')')/2;% Preprocess for suite2p bin (uint16-->int16, although not now)
        end
        fclose(fid);
end
size_3_Y_temp = size(Y_temp,3);
if size_3_Y_temp < 500 % To avoid memory overflow
    backgroundColor_template = median(gpuArray(Y_temp),'all');
else
    %backgroundColor_template = median(gpuArray(Y_temp(:,:,end-500+1:end)),'all');

    temp_Y_temp_interleaved = 1:ceil(init_batch/500):init_batch;    
    backgroundColor_template = median(gpuArray(Y_temp(:,:,temp_Y_temp_interleaved)),'all');    
end
backgroundColor_template = gather(backgroundColor_template);


% % Get a bool mask to denoise images.
% Y_temp = reshape(Y_temp,[512,512*size_3_Y_temp]);
% %BW_raw = gather(Y_temp > 200);
% BW_raw = gather(Y_temp > 100);
% BW = bwpack(BW_raw);
% 
% nhood = [1 1];
% BW1 = imerode(BW,nhood,'ispacked',512);
% BW1_2 = imdilate(BW1,nhood,'ispacked');
% nhood = [1; 1];
% BW2 = imerode(BW,nhood,'ispacked',512);
% BW2_2 = imdilate(BW2,nhood,'ispacked');
% nhood = [0 1; 1 0];
% BW3 = imerode(BW,nhood,'ispacked',512);
% BW3_2 = imdilate(BW3,nhood,'ispacked');
% nhood = [1 0; 0 1];
% BW4 = imerode(BW,nhood,'ispacked',512);
% BW4_2 = imdilate(BW4,nhood,'ispacked');
% 
% closeBW = bitor(bitor(bitor(BW1_2,BW2_2),BW3_2),BW4_2);
% se = strel('rectangle',[3 3]);
% closeBW_dilate = bwunpack(imdilate(closeBW,se,'ispacked'),512);
% 
% Y_temp = Y_temp.*closeBW_dilate + backgroundColor_template.*(~closeBW_dilate);
% Y_temp = reshape(Y_temp,[512,512,size_3_Y_temp]);
% 
% clear BW* closeBW*

if options.if_imageAverage == 1
    Y_temp_mean = squeeze(mean(reshape(Y_temp,[512,512,options.imageAverageBin,temp_init_batch/options.imageAverageBin]), 3));
    Y_temp = Y_temp_mean;
    Y_temp_mean = [];
end
% clear init_Y_index interval
clear interval

data_type = 'uint16';

% fprintf('Registering the random %i frames just to obtain a good template....',init_batch);
% fprintf('Registering the interleaved %i frames just to obtain a good template...',init_batch);
fprintf('Registering the first %i frames just to obtain a good template...',init_batch);

phase_flag = options.phase_flag;
bin_width_forTemplate = 22;%8-->15-->22
bin_totalNum_forTemplate = length(1:bin_width_forTemplate:size_3_Y_temp);

Y_temp_corrected = Y_temp;


% profile viewer
% t0_template = tic;

count_template_threshold = 3;
corr_template_threshold = 0.32;%0.35

count_template = 0;
flag_get_template = false;

test_template_in = cell(1,count_template_threshold);
test_corr_top_measure = zeros(1, count_template_threshold);

while flag_get_template == false
    count_template = count_template + 1;
    
    if count_template >= 2
        init_Y_index = init_Y_index - temp_init_batch;
        
        fid = fopen(Y);
        imsize = 512*512*2;
        current_seek = ftell(fid);
        fseek(fid, 0, 1);
        file_length = ftell(fid);
        fseek(fid, current_seek, -1);
        frame_length = file_length/imsize;
        bstring = 'uint16=>single';
        for tempi=1:temp_init_batch
            t = init_Y_index(tempi);
            feek_indicator = (t-1)*imsize;
            fseek(fid,feek_indicator,'bof');
            %Y_temp(:,:,tempi) = fread(fid,[512,512],'uint16=>single',0,'l')';
            Y_temp(:,:,tempi) = (fread(fid,[512,512],'uint16=>single',0,'l')')/2;% Preprocess for suite2p bin (uint16-->int16, although not now)
        end
        fclose(fid);
        
        size_3_Y_temp = size(Y_temp,3);
    end
    
    if size_3_Y_temp < 500 % To avoid memory overflow
        template_in = median(gpuArray(Y_temp),3);
    else
        %template_in = median(gpuArray(Y_temp(:,:,1:500)),3);
        temp_Y_temp_interleaved = 1:ceil(init_batch/500):init_batch;
        %template_in = median(gpuArray(Y_temp(:,:,temp_Y_temp_interleaved)),3);
        template_in = mean(gpuArray(Y_temp(:,:,temp_Y_temp_interleaved)),3);
    end    
    row_shift_template = zeros(1,size_3_Y_temp,'single','gpuArray');
    col_shift_template = zeros(1,size_3_Y_temp,'single','gpuArray');
    
    for temp_iterIndex=1:options.init_iter
        fprintf('Iter %d,',temp_iterIndex);
        fftTemp_mat = fftn(template_in);
        
        for bin_count = 1:bin_totalNum_forTemplate
            t = (bin_count-1)*bin_width_forTemplate + 1;
            size_3_Y_temp_batch = min(bin_width_forTemplate,size_3_Y_temp-t+1);
            
            Y_temp_batch = gpuArray(Y_temp(:,:,t:(t+size_3_Y_temp_batch-1)));
            lY = size_3_Y_temp_batch;
            
            [minY_batch,maxY_batch] = bounds(Y_temp_batch,[1 2]);
            median_minY_batch = median(minY_batch,'all');
            median_maxY_batch = median(maxY_batch,'all');
            fftY_batch = fft2(Y_temp_batch);
            
            %% dftregistration_min_max
            temp_min_shift = [-20 -20];%20
            temp_max_shift = [20 20];%20
            
            Nr2 = ifftshift(gpuArray(-512:511));
            Nc2 = ifftshift(gpuArray(-512:511));
            tempNr = gpuArray(Nr{1})/512;
            tempNc = gpuArray(Nc{1})/512;
            
            tempNr_batch = repmat(tempNr,[1 1 lY]);
            tempNc_batch = repmat(tempNc,[1 1 lY]);
            
            % Compute cross correlation in frequency domain
            if phase_flag == 1
                %conj1 = fftTemp_mat.*conj(fftY_batch);
                conj2 = fftY_batch.*conj(fftTemp_mat);
                %conj1_pad = sign(FTpad_gpuBatch_jjb_v5(conj1,lY));
                %CC = ifft2(conj1_pad);
                %CCreal = real(CC);
                CCreal = real(ifft2(sign(FTpad_gpuBatch_jjb_v5(fftTemp_mat.*conj(fftY_batch),lY))));%real
            elseif phase_flag == 0
                conj2 = fftY_batch.*conj(fftTemp_mat);
                CCreal = real(ifft2(FTpad_gpuBatch_jjb_v5(fftTemp_mat.*conj(fftY_batch),lY)));%real
            end
            
            [~, III] = max(CCreal, [], [1 2],'linear');
            III = reshape(III,1,lY);
            [row_shift_multi,col_shift_multi,~] = ind2sub([1024 1024 lY],III); %#ok<*ASGLU>
            
            temp_boolIndex = Nr2(row_shift_multi)/2 > temp_max_shift(1) | Nc2(col_shift_multi)/2 > temp_max_shift(2) |...
                Nr2(row_shift_multi)/2 < temp_min_shift(1) | Nc2(col_shift_multi)/2 < temp_min_shift(2);
            
            if sum(temp_boolIndex) > 0
                %temp_boolArray = false(1024,1024,lY,'gpuArray');
                
                temp1 = Nr2/2>temp_max_shift(1);
                temp2 = Nc2/2>temp_max_shift(2);
                temp3 = Nr2/2<temp_min_shift(1);
                temp4 = Nc2/2<temp_min_shift(2);
                
                temp13 = temp1 | temp3;
                temp24 = temp2 | temp4;
                temp13_multi = repmat(temp13',1,1024);
                temp24_multi = repmat(temp24,1024,1);
                temp1234_multi = temp13_multi | temp24_multi;
                
                CCreal2 = CCreal .* ~temp1234_multi;
                
                [~, III] = max(CCreal2, [], [1 2],'linear');
                III = reshape(III,1,lY);
                [row_shift_multi,col_shift_multi,~] = ind2sub([1024 1024 lY],III);
            end
            
            
            row_shift_multi = Nr2(row_shift_multi)/2;
            col_shift_multi = Nc2(col_shift_multi)/2;
            
            % Upsampling, refine estimate with matrix multiply DFT
            row_shift_multi = round(row_shift_multi*us_fac)/us_fac;
            col_shift_multi = round(col_shift_multi*us_fac)/us_fac;
            
            roff_forRefine = 37.5-row_shift_multi*us_fac;
            coff_forRefine = 37.5-col_shift_multi*us_fac;
            
            temp_array = 0:74;
            roff_forRefine_off = repmat(temp_array',1,lY) - roff_forRefine;
            roff_forRefine_off = reshape(roff_forRefine_off,1,75,lY);
            coff_forRefine_off = repmat(temp_array',1,lY) - coff_forRefine;
            coff_forRefine_off = reshape(coff_forRefine_off,1,75,lY);
            
            temp_ifftshift_off = gpuArray(single(ifftshift(0:512-1) - 256));
            temp_ifftshift_off = repmat(temp_ifftshift_off,1,1,lY);
            temp_constant = gpuArray(-1i*2*pi/(512*us_fac));
            temp_ifftshift_off_Xconst = temp_ifftshift_off*temp_constant;
            
            % Upsampled DFT by matrix multiplies
            kernc = exp(pagemtimes(temp_ifftshift_off_Xconst,'transpose',coff_forRefine_off,'none'));
            kernr = exp(pagemtimes(roff_forRefine_off,'transpose',temp_ifftshift_off_Xconst,'none'));
            %CC=kernr*conj2*kernc;
            %CC = pagemtimes(pagemtimes(kernr,conj2),kernc);
            
            %CCreal = real(CC);%real
            CCreal = real(pagemtimes(pagemtimes(kernr,conj2),kernc));%real
            [~, III] = max(CCreal, [], [1 2],'linear');
            III = reshape(III,1,lY);
            [rloc_multi,cloc_multi,~] = ind2sub([75 75 lY],III); %#ok<*ASGLU>
            
            rloc_multi = rloc_multi - 38.5;
            cloc_multi = cloc_multi - 38.5;
            row_shift_multi = row_shift_multi + rloc_multi/us_fac;
            col_shift_multi = col_shift_multi + cloc_multi/us_fac;
            
            if temp_iterIndex == 1
                temp_index = t:(t+size_3_Y_temp_batch-1);
                temp_range = 1:size_3_Y_temp;
                temp_ismember = ismember(temp_range,temp_index);
                for tempi=1:lY
                    row_shift_template = row_shift_template.*(~temp_ismember) + row_shift_multi(tempi).*temp_ismember;
                    col_shift_template = col_shift_template.*(~temp_ismember) + col_shift_multi(tempi).*temp_ismember;
                end
            end
            
            col_InfShift = 512*lY*2;            
            
            displacementField_1d_col = col_shift_multi * -1;
            displacementField_1d_col_int = ceil(abs(displacementField_1d_col)).*sign(displacementField_1d_col);
            col_leftBound_multi = max(1,1 - displacementField_1d_col_int);
            col_rightBound_multi = min(512,512-displacementField_1d_col_int);
            
            temp_grid_col = 1:512;
            temp_grid_col = repmat(temp_grid_col,[512 1]);
            temp_grid_col2 = reshape(temp_grid_col',[1 1 512 512]);
            col_leftBound_multi2 = col_leftBound_multi';
            col_rightBound_multi2 = col_rightBound_multi';
            temp_bool = temp_grid_col2>=col_leftBound_multi2 & temp_grid_col2<=col_rightBound_multi2;
            
            displacementField_3d_col = repmat(col_shift_multi' * -1,[1,1,512,512]);
            displacementField_3d_col_valid = displacementField_3d_col.*temp_bool + col_InfShift.*(~temp_bool);
            
            displacementField_3d_row = repmat(row_shift_multi' * -1,[1,1,512,512]);
            
            displacementField = [displacementField_3d_col_valid,displacementField_3d_row];
            displacementField = permute(displacementField, [4 3 1 2]);
            displacementField = reshape(displacementField,[512,512*lY,2]);
            
            % Shift reconstruct
            Mf = reshape(Y_temp_batch,[512,512*lY]);
            Mf = imwarp(Mf,displacementField,'linear');
            Mf = reshape(Mf,512,512,lY);
            
            overMax_index = Mf > maxY_batch;
            underMin_index = Mf < minY_batch;
            Mf = Mf.*(~overMax_index) + maxY_batch.*overMax_index;
            Mf = Mf.*(~underMin_index) + minY_batch.*underMin_index;
            
            template_in = template_in*(t-lY)/t + mean(Mf,3)*lY/t;
            
            a = 1;
            Y_temp_corrected(:,:,t:(t+size_3_Y_temp_batch-1)) = gather(Mf);
            a = 1;
        end
        
    end
    row_shift_template = gather(row_shift_template);
    col_shift_template = gather(col_shift_template);
            
    temp_corr = motion_metrics_gpuBatch_jjb_v1(Y_temp_corrected,10);
    [B,I] = sort(temp_corr,'descend');
    tempIndex_top = I(1:min(1000,init_batch));%500-->1000-->2000
    Y_temp_corrected_top = Y_temp_corrected(:,:,tempIndex_top);
    template_in = mean(Y_temp_corrected_top,3);
        
    temp_corr_top = temp_corr(tempIndex_top);
    temp_corr_top_mean = mean(temp_corr_top);
    temp_corr_top_median = median(temp_corr_top);
    temp_corr_top_measure = max(temp_corr_top_mean,temp_corr_top_median);
    
    if temp_corr_top_measure > corr_template_threshold
        flag_get_template = true;
        fprintf('temp_corr_top_measure=%.3f, compute template over.\n',temp_corr_top_measure);
    else
        if count_template < count_template_threshold
            fprintf('temp_corr_top_measure=%.3f, compute template again.\n',temp_corr_top_measure);
        else
            flag_get_template = true;
            fprintf('temp_corr_top_measure=%.3f, count_template=%d, compute template must stop.\n',temp_corr_top_measure,count_template);
        end
    end
    
    test_template_in{count_template} = template_in;
    test_corr_top_measure(count_template) = temp_corr_top_measure;
    
    a = 1;
end

a = 1;
[temp_corr_top_measure,I] = max(test_corr_top_measure);
template_in = test_template_in{I};

a = 1;

clear Y_temp
clear Y_temp_corrected_top
clear Y_temp_corrected


clear CCreal CCreal2 cloc_multi coff_forRefine coff_forRefine_off col_shift_multi conj2
clear fftY_batch III kernc kernr maxY_batch
clear median_maxY_batch median_minY_batch Mf  
clear minY_batch Nc2 Nr2 rloc_multi roff_forRefine roff_forRefine_off row_shift_multi 
clear temp1 temp1234_multi temp13 temp13_multi temp2 temp24 temp24_multi temp3 temp4
clear temp_boolArray temp_boolIndex temp_col_shift_multi temp_constant temp_ifftshift_off temp_ifftshift_off_Xconst
clear temp_row_shift_multi tempNc tempNc_batch tempNr tempNr_batch Y_temp_batch

clear col_shift_template overMax_index row_shift_template underMin_index
clear col_leftBound_multi col_leftBound_multi2 col_rightBound_multi col_rightBound_multi2
clear displacementField displacementField_1d_col displacementField_1d_col_int displacementField_3d_col
clear displacementField_3d_col_valid displacementField_3d_row temp_bool
a = 1;
    
fprintf('...done, time=%.2f seconds.\n',toc(t0_template)); 


if options.if_innerSessionNonRigid == 1
    N = 16;%4(185 secs)-->16(290 secs)-->12(250 secs)
    k1 = 3000;%3000-->1500-->3000
    k2 = 300;%300
    
    temp_options = struct;
    temp_options.Nr = Nr;
    temp_options.Nc = Nc;
    
    targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
    innerSessionNonRigid_Name_v = autoGetFunName('innerSessionNonRigid', [targetPATH '\functions']);
    fun_innerSessionNonRigid = str2func(innerSessionNonRigid_Name_v);

    t0_innerSessionNonRigid = tic;
    [D_T_to_nk2, template_nk2, template_nk2_corrected] = ...
        fun_innerSessionNonRigid(template_in,Y,T,options,N,k1,k2,temp_options);
    fprintf('innerSessionNonRigid done, time=%.2f seconds.\n',toc(t0_innerSessionNonRigid));
    a = 1;
    
    [parts_pathstr,~,~] = fileparts(options.h5_filename); %#ok<*ASGLU>
    template_fileFullName = [parts_pathstr,'\template_2.raw'];
    template_nk2_fileFullName = [parts_pathstr,'\template_nk2.raw'];
    template_nk2_corrected_fileFullName = [parts_pathstr,'\template_nk2_corrected.raw'];
    
    a = 1;
    
    template_2 = zeros(512,512,2);
    template_2(:,:,1) = template_in;
    template_2(:,:,2) = template_in;
    template_2_T = pagetranspose(template_2);
    
    fid = fopen(template_fileFullName,'w');
    fwrite(fid,template_2_T,'uint16',0,'l');
    fclose(fid);
    
    template_nk2_T = pagetranspose(template_nk2);
    fid = fopen(template_nk2_fileFullName,'w');
    fwrite(fid,template_nk2_T,'uint16',0,'l');
    fclose(fid);
    
    template_nk2_corrected_T = pagetranspose(template_nk2_corrected);
    fid = fopen(template_nk2_corrected_fileFullName,'w');
    fwrite(fid,template_nk2_corrected_T,'uint16',0,'l');
    fclose(fid);

    a = 1;
    
    T_n = round(linspace(1,T,N+1));
    T_n_offset = T_n+mean(T_n(1:2));
    
    weight_D_T_to_nk2 = zeros(T,N);
    for tempi=1:T
        if tempi < mean(T_n(1:2))
            weight_D_T_to_nk2(tempi,1) = 1;
        elseif tempi > mean(T_n(end-1:end))
            weight_D_T_to_nk2(tempi,N) = 1;
        else
            tempIndex = find(tempi>T_n_offset,1,'last');
            
            temp_dis1 = abs(tempi-T_n_offset(tempIndex));
            temp_dis2 = abs(tempi-T_n_offset(tempIndex+1));
            
            weight_D_T_to_nk2(tempi,tempIndex) = temp_dis2/(temp_dis1+temp_dis2);
            weight_D_T_to_nk2(tempi,tempIndex+1) = 1-weight_D_T_to_nk2(tempi,tempIndex);
        end
    end
    weight_D_T_to_nk2;
    D_T_to_nk2;
    
    D_T_to_nk2_merged = zeros(512,512,2*N,'single');
    for n=1:N
        temp_range = (2*(n-1)+1):2*n;
        D_T_to_nk2_merged(:,:,temp_range) = D_T_to_nk2{n};
    end    
    
end

if T > init_batch
    % To get shiftBounds from interleaved init_batch_shift frames
    init_batch_shift = init_batch;
    %shift_templateBool = getShiftBounds_v3(Y,init_batch_shift,T,options);
    shift_templateBool = getShiftBounds_v4(Y,init_batch_shift,T,options);
    
elseif T <= init_batch
    % To get shiftBounds from template
    outlierIndex_row_shift_template = isoutlier(row_shift_template,'grubbs');
    outlierIndex_col_shift_template = isoutlier(col_shift_template,'grubbs');
    row_shift_template_valid = row_shift_template(~outlierIndex_row_shift_template);
    col_shift_template_valid = col_shift_template(~outlierIndex_col_shift_template);
    
    [min_row_shift_template, max_row_shift_template] = bounds(row_shift_template_valid);
    [min_col_shift_template, max_col_shift_template] = bounds(col_shift_template_valid);
    shiftBounds_template = [min_row_shift_template,max_row_shift_template,...
        min_col_shift_template,max_col_shift_template];
    
    clear row_shift_template col_shift_template
    clear outlierIndex_row_shift_template outlierIndex_col_shift_template
    clear row_shift_template_valid col_shift_template_valid
    
    % To get shift_templateBool
    temp_positiveIndex = shiftBounds_template>0;
    temp_negativeIndex = shiftBounds_template<0;
    shiftBounds_template(temp_positiveIndex) = ceil(shiftBounds_template(temp_positiveIndex));
    shiftBounds_template(temp_negativeIndex) = floor(shiftBounds_template(temp_negativeIndex));
    
    row_shift_templateBool = true(512,512,'gpuArray');
    col_shift_templateBool = true(512,512,'gpuArray');
    row_shift_templateBool(513+shiftBounds_template(1):512,:) = false;
    row_shift_templateBool(1:shiftBounds_template(2),:) = false;
    col_shift_templateBool(:,513+shiftBounds_template(3):512) = false;
    col_shift_templateBool(:,1:shiftBounds_template(4)) = false;
    shift_templateBool = row_shift_templateBool & col_shift_templateBool;
    
    clear row_shift_templateBool col_shift_templateBool outputExistingTemplate
end



% template_in(template_in<0) = 0;
positive_index = (template_in > 0);
template_in = template_in.*positive_index;
template_in_raw = template_in;
% outputExistingTemplate = zeros(1,4); %#ok<*PREALL>
if options.if_useExistingTemplate == 1
    fprintf('Now use existing template!');
    %template_in, template_in_raw; % This is current template.
    %template, template_n11n; % This is existing template.
    
    [min_template_in_raw,max_template_in_raw] = bounds(gather(template_in_raw),'all');
    [min_template,max_template] = bounds(gather(template),'all');
    median_template = median(template,'all');
     
    template_n11n = rescale(template,min_template_in_raw,max_template_in_raw);
    a = 1;
    
    if options.if_nonRigidExistingTemplate == 1
        options_nonrigid = NoRMCorreSetParms(...
            'd1',512,...
            'd2',512,...
            'grid_size',[32,32],...%48
            'overlap_pre',[32,32],...%48
            'overlap_post',[32,32],...%48
            'us_fac',50,...%50
            'max_dev',[12 12],...%8
            'max_shift',[48 48],...  %[48 48]
            'phase_flag',false,... %true-->false
            'upd_template',false,... %true
            'iter',3,... %2
            'boundary','zero',...
            'output_type','mat',...
            'correct_bidir',true);                     
        
        options_nonrigid.shifts_method = 'cubic';
        options_nonrigid.iter = 1;
        
        %gather_template_in = gather(template_in);
        %gather_template_n11n = gather(template_n11n);  

    optimal_NumTiles_histeq = 8;
        
    % To fit nBins_histeq, which is for histogram equlization
    gather_template_n11n = gather(template_n11n);       
    gather_template_in = gather(template_in);
    initial_template_n11n = rescale(gather_template_n11n,0,1);
    initial_template_in = rescale(gather_template_in,0,1);

    temp_start = 8;
    temp_end = 13;
    temp_step = 0.05;
    tempN = floor((temp_end-temp_start)/temp_step);
    temp_nBins_histeq_power = temp_start:temp_step:temp_end;
    rho = zeros(1,tempN);

    temp_options_nonrigid = options_nonrigid; 
    t0_temp = tic;
    parfor tempi=1:tempN
        temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',floor(2^(temp_nBins_histeq_power(tempi))),'Range','original','Distribution','rayleigh','Alpha',0.4);
        temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',floor(2^(temp_nBins_histeq_power(tempi))),'Range','original','Distribution','rayleigh','Alpha',0.4);       
        temp_output = normcorre_jjb_v2(temp_template_in,temp_options_nonrigid,temp_template_n11n);
        rho(tempi) = corr2(temp_output,temp_template_n11n);
    end
    [~,I]=max(rho);
    optimal_nBins_histeq_power = temp_nBins_histeq_power(I);
    optimal_nBins_histeq = floor(2^optimal_nBins_histeq_power);
    fprintf(' optimal_nBins_histeq=%d, time=%.1f secs.',optimal_nBins_histeq,toc(t0_temp));

    temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',0.4);
    temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',0.4);
    gather_template_n11n = rescale(temp_template_n11n,min_template_in_raw,max_template_in_raw); %#ok<*NASGU>
    gather_template_in = rescale(temp_template_in,min_template_in_raw,max_template_in_raw);


    % To fit alpha_histeq, which is for histogram equlization    
    gather_template_n11n = gather(template_n11n);
    gather_template_in = gather(template_in);
    initial_template_n11n = rescale(gather_template_n11n,0,1);
    initial_template_in = rescale(gather_template_in,0,1);

    temp_start = 0.1;%0.1
    temp_end = 1;%10-->5-->3-->1
    temp_step = 0.01;%0.1-->0.025-->0.01
    tempN = floor((temp_end-temp_start)/temp_step);
    temp_alpha_histeq = temp_start:temp_step:temp_end;
    rho = zeros(1,tempN);

    temp_options_nonrigid = options_nonrigid; 
    t0_temp = tic;
    parfor tempi=1:tempN
        temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',temp_alpha_histeq(tempi));
        temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',temp_alpha_histeq(tempi));       
        temp_output = normcorre_jjb_v2(temp_template_in,temp_options_nonrigid,temp_template_n11n);
        rho(tempi) = corr2(temp_output,temp_template_n11n);
    end
    [~,I]=max(rho);
    optimal_alpha_histeq = temp_alpha_histeq(I);
    fprintf(' optimal_alpha_histeq=%.3f, time=%.1f secs.',optimal_alpha_histeq,toc(t0_temp));

    temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);
    temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);
    gather_template_n11n = rescale(temp_template_n11n,min_template_in_raw,max_template_in_raw);
    gather_template_in = rescale(temp_template_in,min_template_in_raw,max_template_in_raw);        
        
    
    % To fit NumTiles_histeq, which is for histogram equlization
    gather_template_n11n = gather(template_n11n);       
    gather_template_in = gather(template_in);
    initial_template_n11n = rescale(gather_template_n11n,0,1);
    initial_template_in = rescale(gather_template_in,0,1);

    temp_start = 8;%2
    temp_end = 64;%128
    temp_step = 1;
    tempN = floor((temp_end-temp_start)/temp_step);
    temp_NumTiles_histeq = temp_start:temp_step:temp_end;
    rho = zeros(1,tempN);

    temp_options_nonrigid = options_nonrigid; 
    t0_temp = tic;
    parfor tempi=1:tempN
        temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*temp_NumTiles_histeq(tempi),'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);
        temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*temp_NumTiles_histeq(tempi),'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);       
        temp_output = normcorre_jjb_v2(temp_template_in,temp_options_nonrigid,temp_template_n11n);
        rho(tempi) = corr2(temp_output,temp_template_n11n);
    end
    [~,I]=max(rho);
    optimal_NumTiles_histeq = temp_NumTiles_histeq(I);
    fprintf(' optimal_NumTiles_histeq=%d, time=%.1f secs.',optimal_NumTiles_histeq,toc(t0_temp));

    temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);
    temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);
    gather_template_n11n = rescale(temp_template_n11n,min_template_in_raw,max_template_in_raw); %#ok<*NASGU>
    gather_template_in = rescale(temp_template_in,min_template_in_raw,max_template_in_raw);
    
        
        % To fit optimal grid size
        temp_start = 24;%32-->24-->16
        temp_end = 64;%100
        temp_step = 1;%1
        tempN = floor((temp_end-temp_start)/temp_step);
        temp_grid_size = temp_start:temp_step:temp_end;
        rho = zeros(1,tempN);
        
        t0_temp = tic;
        parfor tempi=1:tempN
           temp_options_nonrigid = options_nonrigid;
           temp_options_nonrigid.grid_size = [temp_grid_size(tempi),temp_grid_size(tempi),1];
           temp_options_nonrigid.overlap_pre = [temp_grid_size(tempi),temp_grid_size(tempi),1];
           temp_options_nonrigid.overlap_post = [temp_grid_size(tempi),temp_grid_size(tempi),1];
           temp_output = normcorre_jjb_v2(gather_template_in,temp_options_nonrigid,gather_template_n11n);                          
           rho(tempi) = corr2(temp_output,gather_template_n11n);
        end
        [~,I]=max(rho);
        optimal_grid_size = temp_grid_size(I);
        optimal_overlap = optimal_grid_size;
        options_nonrigid.grid_size = [optimal_grid_size,optimal_grid_size,1];
        options_nonrigid.overlap_pre = [optimal_overlap,optimal_overlap,1];
        options_nonrigid.overlap_post = [optimal_overlap,optimal_overlap,1];        
        fprintf(' optimal_grid_size=%d, time=%.1f.',optimal_grid_size,toc(t0_temp));
        
        
        % To fit optimal max_dev
        temp_start = 2;
        temp_end = 15;
        temp_step = 1;
        tempN = floor((temp_end-temp_start)/temp_step);
        temp_max_dev = temp_start:temp_step:temp_end;
        rho = zeros(1,tempN);
        
        t0_temp = tic;
        parfor tempi=1:tempN
           temp_options_nonrigid = options_nonrigid;
           temp_options_nonrigid.max_dev = [temp_max_dev(tempi),temp_max_dev(tempi),1];
           temp_output = normcorre_jjb_v2(gather_template_in,temp_options_nonrigid,gather_template_n11n);                          
           rho(tempi) = corr2(temp_output,gather_template_n11n);
        end
        [~,I]=max(rho);      
        optimal_max_dev = temp_max_dev(I);
        options_nonrigid.max_dev = [optimal_max_dev,optimal_max_dev,1];
        fprintf(' optimal_max_dev=%d, time=%.1f.\n',optimal_max_dev,toc(t0_temp));
        
        
        % To fit optimal contrast_ratio by fminsearch
        range_contrast_ratio = [0.3,3.01];%[0.99,4]-->[0.99,3.01]-->[0.1,3.01]-->[0.3,3.01]
        params_range = range_contrast_ratio;
                
        ModelName = 'nonRigidParamFit_v2';
        Model = str2func(ModelName);
        func = @(params) Model(params,gather_template_in,options_nonrigid,gather_template_n11n,params_range);
        
        % Use fminsearch
        optimset_options = optimset('fminsearch');
        optimset_options.TolFun = 1e-7;%1e-6
        optimset_options.TolX = 1e-3;%1e-1
        %optimset_options.MaxFunEvals = 40;
        optimset_options.MaxIter = 25;%40
        optimset_options.Display = 'off'; % off, iter, final
        optimset_options.FunValCheck = 'on';
        optimset_options.PlotFcns = [];%@optimplotfval
        
        t0 = tic;
                
        initial_contrast_ratio = 1:0.28570:3;%1:0.1428:2
        contrast_ratio_raw = ...
           -log(1./((initial_contrast_ratio-range_contrast_ratio(1))./(range_contrast_ratio(2)-range_contrast_ratio(1)))-1);
        initialParams = contrast_ratio_raw;
        tempN = length(initialParams);
        rho = zeros(1,tempN);
        temp_modelParams = zeros(1,tempN);
        parfor tempi=1:tempN   
        %for tempi=1:tempN   
            temp_optimset_options = optimset_options;
            if tempi==1
                temp_optimset_options.Display = 'iter'; % off, iter, final
            end
            temp_modelParams(tempi) = fminsearch(func, initialParams(tempi), temp_optimset_options);
            [~,~,rho(tempi)] = Model(temp_modelParams(tempi),gather_template_in,options_nonrigid,gather_template_n11n,params_range);                                  
        end
        a = 1;
        [rho2,I]=max(rho);
        modelParams = temp_modelParams(I);
        
        
        contrast_ratio_raw = modelParams(1);
        contrast_ratio = (1/(1+exp(-contrast_ratio_raw)))*(range_contrast_ratio(2)-range_contrast_ratio(1))+range_contrast_ratio(1);% constrain in (1, 6)
        fprintf('fimsearch time=%.1f seconds, rho=%.5f, contrast_ratio=%.3f(from%.3f)\n',...
            toc(t0),rho2,contrast_ratio,initial_contrast_ratio(I));
        
        gather_template_in = gather_template_in.^(1/contrast_ratio);
        gather_template_n11n = gather_template_n11n.^(1/contrast_ratio);
        [min_gather_template_in,max_gather_template_in] = bounds(gather_template_in,'all');

        % To do coarse correction firstly
        options_nonrigid.iter = 3;
        options_nonrigid.shifts_method = 'FFT';
        M1_final = normcorre_jjb_v2(gather_template_in,options_nonrigid,gather_template_n11n);
        
        overMax_index = M1_final > max_gather_template_in;
        underMin_index = M1_final < min_gather_template_in;
        M1_final = M1_final.*(~overMax_index) + max_gather_template_in.*overMax_index;
        M1_final = M1_final.*(~underMin_index) + min_gather_template_in.*underMin_index;
        
        
        if options.if_nonRigidRefine == 1
            % To do refined correction secondly            
            t0_temp = tic;
            %options_nonrigid.shifts_method = 'cubic';
            options_nonrigid.overlap_pre = [16,16,1];%16
            options_nonrigid.overlap_post = [16,16,1];
            
            temp_start = 64;
            temp_end = 32;
            temp_step = -1;
            rho_0 = corr2(M1_final,gather_template_n11n);
            for tempi=temp_start:temp_step:temp_end %64 to 32
                options_nonrigid.grid_size = [tempi,tempi,1];
                temp_output = normcorre_jjb_v2(M1_final,options_nonrigid,gather_template_n11n);
                temp_rho = corr2(temp_output,gather_template_n11n);
                if temp_rho > rho_0
                    M1_final = temp_output;
                    rho_0 = temp_rho;
                    
                    overMax_index = M1_final > max_gather_template_in;
                    underMin_index = M1_final < min_gather_template_in;
                    M1_final = M1_final.*(~overMax_index) + max_gather_template_in.*overMax_index;
                    M1_final = M1_final.*(~underMin_index) + min_gather_template_in.*underMin_index;                    
                end                
            end
            fprintf('Nonrigid refine time=%.2f seconds.\n',toc(t0_temp));
        end

        
        %M1_final = M1_final.^contrast_ratio;
        


        PyramidLevels = 9;
        iterNum = 1000;
        temp_iterRange = zeros(1,PyramidLevels);
        for tempi=1:PyramidLevels
            if tempi == 1
                temp_iterRange(tempi) = iterNum;
            else
                temp_iterRange(tempi) = ceil(temp_iterRange(tempi-1)/1.5);
            end
        end
        %[D,temp_reg] = imregdemons(gather(template_in),gather(M1_final),...%template_in_raw
        [D,temp_reg] = imregdemons(gather_template_in,gather(M1_final),...
            ceil(iterNum:-iterNum/PyramidLevels:1),...
            'AccumulatedFieldSmoothing',1,...
            'PyramidLevels',9,...%9
            'DisplayWaitbar',false);
        D = single(gpuArray(D));        
        
        if options.if_imageAverage == 1
            temp_bin_width = ceil(bin_width/options.imageAverageBin);
        else
            temp_bin_width = bin_width;
        end
        D_batch = repmat(D,[1,temp_bin_width,1]);
        
        fprintf('Template non-rigid done.\n');
        clear M1_final temp_reg
        %clear overMax_index underMin_index
    end
    
    %clear overMax_index underMin_index
    clear max_template_in_raw min_template_in_raw template_n11n
end

clear positive_index


% template_in = template_in.*shift_templateBool + backgroundColor_template.*(~shift_templateBool);

% comment this line in 20230504
% template_in_raw = template_in_raw.*shift_templateBool + backgroundColor_template.*(~shift_templateBool);

clear Y_temp Y_temp2 Greg M_temp
clear temp_cell Bs nr nc

%%
% maxNumCompThreads(1);
% template = mat2cell_ov(template_in,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY);
template = mat2cell_ov(template_in_raw,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY);
% temp_mat = template_in;
% use_windowing = options.use_windowing;
% phase_flag = options.phase_flag;

fftTemp_mat = fftn(template_in_raw);
clear template_in
clear template_in_raw

% if use_windowing
%     fftTemp = cellfun(@fftn,cellfun(@han,template,'un',0),'un',0);
%     fftTemp_mat = fftn(han(temp_mat));
% else
%     fftTemp = cellfun(@fftn,template,'un',0);
%     fftTemp_mat = fftn(temp_mat);
% end
% 
% clear temp_mat


template{1} = gather(template{1});

if upd_template
    if nd == 2; buffer = mat2cell_ov(zeros(d1,d2,bin_width),xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY); end
    if nd == 3; buffer = mat2cell_ov(zeros(d1,d2,d3,bin_width),xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY); end
end
if ~strcmpi(options.output_type,'mat')
    options.mem_batch_size = max(min(round(options.mem_batch_size/bin_width)*bin_width,T),1);
    %if nd == 2; mem_buffer = zeros(d1,d2,options.mem_batch_size,'single'); end
    %if nd == 3; mem_buffer = zeros(d1,d2,d3,options.mem_batch_size,'single'); end
end



a = 1;
% To gather all gpu data, so that gpu memory will release completely.
backgroundColor_template = gather(backgroundColor_template);
fftTemp_mat = gather(fftTemp_mat);
shift_templateBool = gather(shift_templateBool);
if options.if_useExistingTemplate == 1 && options.if_nonRigidExistingTemplate == 1
    D = gather(D);
    D_batch = gather(D_batch);
end
clear overMax_index underMin_index

a = 1;
if options.if_innerSessionNonRigid == 1
    D_T_to_nk2_merged = gpuArray(D_T_to_nk2_merged);
end
backgroundColor_template = gpuArray(backgroundColor_template);
fftTemp_mat = gpuArray(fftTemp_mat);
shift_templateBool = gpuArray(shift_templateBool);
if options.if_useExistingTemplate == 1 && options.if_nonRigidExistingTemplate == 1
    D = gpuArray(D);
    D_batch = gpuArray(D_batch);
    %D_batch_1d = reshape(D_batch,[512*512*temp_bin_width*2,1]);
    %D_batch_3d = reshape(D_batch,[512,512*temp_bin_width,2]);
end
a = 1;

if_spmd = options.if_spmd;
if if_spmd == 0
    switch lower(options.output_type)
        case 'mat'
            M_final = zeros([sizY,T],data_type);
        case 'memmap'
            M_final = matfile(filename,'Writable',true);
            if nd == 2; M_final.Y(d1,d2,T) = zeros(1,data_type); end
            if nd == 3; M_final.Y(d1,d2,d3,T) = zeros(1,data_type); end
            M_final.Yr(d1*d2*d3,T) = zeros(1,data_type);
        case {'hdf5','h5'}
            [h5_pathstr,h5_fname,h5_ext] = fileparts(options.h5_filename_temporary);
            new_filename = fullfile(h5_pathstr,[h5_fname,'_',options.datestr_now30,h5_ext]); %#ok<*TNOW1,*DATST> 
            options.h5_filename = new_filename;

            
            new_filename_bin = fullfile(h5_pathstr,[h5_fname(1:13),'Bin',h5_fname(14:end),'_',options.datestr_now30,h5_ext]);
            options.h5_filename_bin = new_filename_bin;
            
            M_final = options.h5_filename;
            if nd == 2
                %h5create(options.h5_filename,['/',options.h5_groupname],[d1,d2,Inf],'Chunksize',[d1,d2,options.mem_batch_size],'Datatype',data_type);
                h5create(options.h5_filename,['/',options.h5_groupname],[d1,d2,Inf],'Chunksize',[d1,d2,1],'Datatype',data_type);                
                h5create(options.h5_filename_bin,['/',options.h5_groupname],[d1,d2,Inf],'Chunksize',[d1,d2,1],'Datatype',data_type);
            elseif nd == 3
                h5create(options.h5_filename,['/',options.h5_groupname],[d1,d2,d3,Inf],'Chunksize',[d1,d2,d3,options.mem_batch_size],'Datatype',data_type);
            end
        case {'tif','tiff'}
            M_final = ['motion corrected file has been saved as ', options.tiff_filename];
            opts_tiff.append = true;
            opts_tiff.big = true; %#ok<*STRNU> 
            if nd == 3
                error('Saving volumetric tiff stacks is currently not supported. Use a different filetype');
            end
        otherwise
            error('This filetype is currently not supported')
    end
end

cnt_buf = 0;
fprintf('Template initialization complete.  Now registering all the frames with new template. \n')
%%

prevstr = []; %#ok<*NASGU>
bin_count_forSave = 0;
bin_totalNum = length(1:bin_width:T);


if if_spmd == 1
    if options.if_spmd_singleGPU == 0
        iter = 1;
        upd_template = false;
        
        [h5_pathstr,h5_fname,h5_ext] = fileparts(options.h5_filename_temporary_raw);
        new_filename = fullfile(h5_pathstr,[h5_fname,'_',options.datestr_now30,h5_ext]);
        options.h5_filename = new_filename;
        new_filename_bin = fullfile(h5_pathstr,[h5_fname,'Bin','_',options.datestr_now30,h5_ext]);
        options.h5_filename_bin = new_filename_bin;
        M_final = options.h5_filename;
        [h5_pathstr,h5_fname,h5_ext] = fileparts(options.h5_filename);
        
        workerNum = options.workerNum;
        N = bin_totalNum;
        
        
        imageAverageBin = options.imageAverageBin;        
        spmdCut_global = options.spmdCut_global;
        cpuThread_taskCoeff = options.cpuThread_taskCoeff;
        
        temp1=ones(1, workerNum);
        temp1(spmdCut_global+1:end) = temp1(spmdCut_global+1:end)*cpuThread_taskCoeff;
        temp2 = N/sum(temp1);
        temp3 = ceil(temp1*temp2);
        while sum(temp3) > N
            for tempi=1:spmdCut_global
                temp3(tempi) = temp3(tempi) - 1;
                if sum(temp3) == N
                    break
                end
            end
            if sum(temp3) == N
                break
            end
            for tempi=(spmdCut_global+1):workerNum
                temp3(tempi) = temp3(tempi) - 1;
                if sum(temp3) == N
                    break
                end
            end
        end
        
        i_worker = temp3;
        
        if sum(i_worker) < workerNum
            i_worker = zeros(1, workerNum);
            i_worker(1:N) = 1;
        end
        
        options.i_worker = i_worker;
        
        i_worker_start = zeros(1, workerNum);
        i_worker_end = zeros(1, workerNum);
        for tempj=1:workerNum
            if tempj == 1
                i_worker_start(1) = 1;
            else
                i_worker_start(tempj) = sum(i_worker(1:tempj-1))+1;
            end
            i_worker_end(tempj) = i_worker_start(tempj)+i_worker(tempj)-1;
        end
        spmdindex = [];
        
        %fftTemp{1} = gather(fftTemp{1});                
        %fftTemp_mat = gather(fftTemp{1});
        
        for it = 1:iter
            spmd(0, workerNum)
                %labindex;
                %numlabs;
                t0_spmd = tic;
                spmdindex = labindex;
                str_spmdindex = num2str(spmdindex);
                if spmdindex < 10
                    str_spmdindex = ['0' str_spmdindex]; %#ok<*AGROW>
                end
                
                if spmdindex <= options.spmdCut_global
                    gpu_globalEnableFlag = 1;
                else
                    gpu_globalEnableFlag = 0;
                    gpuDevice([]);
                end
                
                
                if spmdindex <= options.spmdCut_4090
                    %gpuDevice(1);
                else
                    %gpuDevice(2);
                end
                
                temp_h5_filename = fullfile(h5_pathstr,['partial_',h5_fname,'_worker',str_spmdindex,h5_ext]);
                temp_h5_filename_bin = fullfile(h5_pathstr,['partial_',h5_fname(1:13),'Bin',h5_fname(14:end),'_worker',str_spmdindex,h5_ext]);
                
                if exist(temp_h5_filename,'file') ~= 0
                    delete(temp_h5_filename);
                end
                if exist(temp_h5_filename_bin,'file') ~= 0
                    delete(temp_h5_filename_bin);
                end                
                
                %h5create(temp_h5_filename,['/',options.h5_groupname],[d1,d2,Inf],'Chunksize',[d1,d2,options.mem_batch_size],'Datatype',data_type);              
                h5create(temp_h5_filename,['/',options.h5_groupname],[d1,d2,Inf],'Chunksize',[d1,d2,1],'Datatype',data_type);                              
                h5create(temp_h5_filename_bin,['/',options.h5_groupname],[d1,d2,Inf],'Chunksize',[d1,d2,1],'Datatype',data_type);


                
                
                fid = fopen(Y);
                imsize = 512*512*2;
                current_seek = ftell(fid);
                fseek(fid, 0, 1);
                file_length = ftell(fid);
                fseek(fid, current_seek, -1);
                frame_length = file_length/imsize;
                %bstring = 'uint16=>single';
                
                last_temp_progress = 0;
                
                tempi_worker = i_worker_start(spmdindex):i_worker_end(spmdindex);
                tempi_length = length(tempi_worker);
                for tempj=1:tempi_length
                    n = tempi_worker(tempj);
                    bin_count = n;
                    t = (bin_count-1)*bin_width + 1;
                    
                    if spmdindex == workerNum
                        temp_progress = floor(tempj/tempi_length*100);
                        %if temp_progress == floor(temp_progress/4)*4 && ...
                        if temp_progress == floor(temp_progress/6)*6 && ...
                                temp_progress ~= last_temp_progress
                            temp_printFlag = 1;
                        else
                            temp_printFlag = 0;
                        end
                        if temp_printFlag == 1
                            str=[int2str(temp_progress),'%% out of ',num2str(T),' frames registered...'];
                            fprintf(['Now Loding...',str,'\n']);
                        end
                        last_temp_progress = temp_progress;
                        
                        
                        %str=[num2str(tempj/tempi_length*100),'%% out of ',num2str(T),' frames registered...'];
                        %fprintf(['Now Loding...',str,'\n']);
                    end
                    
                    if it == iter
                        bin_count_forSave = bin_count;
                    end
                    
                    % Read Raw file
                    %Ytm = single(read_raw_file(Y,t,min(t+bin_width-1,T)-t+1,[512 512],2));
                    window = min(t+bin_width-1,T)-t+1;                    
                    feek_indicator = (t-1)*imsize;
                    fseek(fid,feek_indicator,'bof');
                    Ytm = gpuArray(fread(fid,512*512*window,'uint16=>single',0,'l')');
                    Ytm = Ytm/2;% Preprocess for suite2p bin (uint16-->int16, although not now)
                    Ytm = reshape(Ytm,[512 512 window]);
                    Ytm = pagetranspose(Ytm);                    
                    size_3_Ytm = size(Ytm,ndims(Ytm));
                    
                    
%                     % Get a bool mask to denoise images.
%                     Ytm = reshape(Ytm,[512,512*size_3_Ytm]);
%                     %BW_raw = gather(Ytm > 200);
%                     BW_raw = gather(Ytm > 100);
%                     BW = bwpack(BW_raw);
%                     
%                     nhood = [1 1];
%                     BW1 = imerode(BW,nhood,'ispacked',512);
%                     BW1_2 = imdilate(BW1,nhood,'ispacked');
%                     nhood = [1; 1];
%                     BW2 = imerode(BW,nhood,'ispacked',512);
%                     BW2_2 = imdilate(BW2,nhood,'ispacked');
%                     nhood = [0 1; 1 0];
%                     BW3 = imerode(BW,nhood,'ispacked',512);
%                     BW3_2 = imdilate(BW3,nhood,'ispacked');
%                     nhood = [1 0; 0 1];
%                     BW4 = imerode(BW,nhood,'ispacked',512);
%                     BW4_2 = imdilate(BW4,nhood,'ispacked');
%                     
%                     closeBW = bitor(bitor(bitor(BW1_2,BW2_2),BW3_2),BW4_2);
%                     closeBW = bitor(BW1_2,BW2_2);
%                     se = strel('rectangle',[3 3]);
%                     closeBW_dilate = bwunpack(imdilate(closeBW,se,'ispacked'),512);
%                     
%                     % Ytm(~closeBW_dilate) = backgroundColor_template; % This could lead to gpu memory leakage!
%                     Ytm = Ytm.*closeBW_dilate + backgroundColor_template.*(~closeBW_dilate);
%                     Ytm = reshape(Ytm,[512,512,size_3_Ytm]);

                                      

                    lY_raw = size_3_Ytm;
                    if options.if_imageAverage == 1
                        if ceil(size_3_Ytm/imageAverageBin)*imageAverageBin == size_3_Ytm
                            Ytm_mean = squeeze(mean(reshape(Ytm,[512,512,imageAverageBin,size_3_Ytm/imageAverageBin]), 3));
                        else
                            Ytm_mean = squeeze(mean(reshape(...
                                Ytm(:,:,1:floor(size_3_Ytm/imageAverageBin)*imageAverageBin),[512,512,imageAverageBin,floor(size_3_Ytm/imageAverageBin)]), 3));
                            tempi = (floor(size_3_Ytm/imageAverageBin)*imageAverageBin+1):size_3_Ytm;
                            Ytm_mean(:,:,(floor(size_3_Ytm/imageAverageBin)+1)) = mean(Ytm(:,:,tempi), 3);
                        end
                        
                        Ytm = Ytm_mean;
                        size_3_Ytm = size(Ytm,ndims(Ytm));
                        Ytm_mean = [];
                    end
                    
                    lY = size_3_Ytm;
                    shifts = struct('shifts',cell(lY,1),'shifts_up',cell(lY,1),'diff',cell(lY,1));
                                        
                    [minY_batch,maxY_batch] = bounds(Ytm,[1 2]);
                    median_minY_batch = median(minY_batch,'all');
                    median_maxY_batch = median(maxY_batch,'all');
                    fftY_batch = fft2(Ytm);
                                        
                    if spmdindex <= options.spmdCut_local
                        temp_fftTemp_mat = gpuArray(fftTemp_mat);
                        temp_fftY_batch = gpuArray(fftY_batch);
                        Nr2 = ifftshift(gpuArray(-512:511));
                        Nc2 = ifftshift(gpuArray(-512:511));
                        tempNr = gpuArray(Nr{1})/512;
                        tempNc = gpuArray(Nc{1})/512;     
                        %Nr2 = ifftshift(-fix(512):ceil(512)-1);
                        %Nc2 = ifftshift(-fix(512):ceil(512)-1);
                        %tempNr = Nr{1}/512;
                        %tempNc = Nc{1}/512;                            
                    else
                        temp_fftTemp_mat = fftTemp_mat;
                        temp_fftY_batch = fftY_batch;
                        Nr2 = ifftshift(-fix(512):ceil(512)-1);
                        Nc2 = ifftshift(-fix(512):ceil(512)-1);
                        tempNr = Nr{1}/512;
                        tempNc = Nc{1}/512;                                                 
                    end
                    tempNr_batch = repmat(tempNr,[1 1 lY]);
                    tempNc_batch = repmat(tempNc,[1 1 lY]);
                    
                    
                    if options.if_dftreg_bacth == 0   
                        
                        %if spmdindex <= options.spmdCut_local
                        if gpu_globalEnableFlag == 1   
                                
                            % Compute cross correlation in frequency domain
                            if phase_flag == 1
                                %conj1 = temp_fftTemp_mat.*conj(temp_fftY_batch);
                                conj2 = temp_fftY_batch.*conj(temp_fftTemp_mat);
                                %conj1_pad = sign(FTpad_gpuBatch_jjb_v5(conj1,lY));
                                %CC = ifft2(conj1_pad);
                                %CCreal = real(CC);
                                CCreal = real(ifft2(sign(FTpad_gpuBatch_jjb_v5(temp_fftTemp_mat.*conj(temp_fftY_batch),lY))));
                            elseif phase_flag == 0
                               conj2 = temp_fftY_batch.*conj(temp_fftTemp_mat);
                               CCreal = real(ifft2(FTpad_gpuBatch_jjb_v5(temp_fftTemp_mat.*conj(temp_fftY_batch),lY)));
                            end
                            
                            [~, III] = max(CCreal, [], [1 2],'linear');
                            III = reshape(III,1,lY);
                            [row_shift_multi,col_shift_multi,~] = ind2sub([1024 1024 lY],III); %#ok<*ASGLU>
                            
                            row_shift_multi = Nr2(row_shift_multi)/2;
                            col_shift_multi = Nc2(col_shift_multi)/2;
                            
                             % Upsampling, refine estimate with matrix multiply DFT
                            row_shift_multi = round(row_shift_multi*us_fac)/us_fac;
                            col_shift_multi = round(col_shift_multi*us_fac)/us_fac;
                            
                            roff_forRefine = 37.5-row_shift_multi*us_fac;
                            coff_forRefine = 37.5-col_shift_multi*us_fac;
                            
                            temp_array = 0:74;
                            roff_forRefine_off = repmat(temp_array',1,lY) - roff_forRefine;
                            roff_forRefine_off = reshape(roff_forRefine_off,1,75,lY);
                            coff_forRefine_off = repmat(temp_array',1,lY) - coff_forRefine;
                            coff_forRefine_off = reshape(coff_forRefine_off,1,75,lY);
                            
                            temp_ifftshift_off = gpuArray(single(ifftshift(0:512-1) - 256));
                            temp_ifftshift_off = repmat(temp_ifftshift_off,1,1,lY);
                            temp_constant = gpuArray(-1i*2*pi/(512*us_fac));
                            temp_ifftshift_off_Xconst = temp_ifftshift_off*temp_constant;
                            
                            % Upsampled DFT by matrix multiplies
                            kernc = exp(pagemtimes(temp_ifftshift_off_Xconst,'transpose',coff_forRefine_off,'none'));
                            kernr = exp(pagemtimes(roff_forRefine_off,'transpose',temp_ifftshift_off_Xconst,'none'));
                            %CC=kernr*conj2*kernc;
                            %CC = pagemtimes(pagemtimes(kernr,conj2),kernc);                            
                            %CCreal = real(CC);
                            CCreal = real(pagemtimes(pagemtimes(kernr,conj2),kernc));
                            
                            [~, III] = max(CCreal, [], [1 2],'linear');
                            III = reshape(III,1,lY);
                            [rloc_multi,cloc_multi,~] = ind2sub([75 75 lY],III); %#ok<*ASGLU>
                            rloc_multi = rloc_multi - 38.5;
                            cloc_multi = cloc_multi - 38.5;
                            row_shift_multi = row_shift_multi + rloc_multi/us_fac;
                            col_shift_multi = col_shift_multi + cloc_multi/us_fac;
                            
                            %temp_row_shift_multi = reshape(row_shift_multi,1,1,lY);
                            %temp_row_shift_multi = repmat(temp_row_shift_multi,512,512,1);
                            %temp_col_shift_multi = reshape(col_shift_multi,1,1,lY);
                            %temp_col_shift_multi = repmat(temp_col_shift_multi,512,512,1);
                            
                            % Shift reconstruct
%                             Mf = real(ifft2(...
%                                 temp_fftY_batch.*...
%                                 exp(...
%                                 1i*2*pi*(-temp_row_shift_multi.*tempNr_batch-temp_col_shift_multi.*tempNc_batch)...
%                                 )));
                            
%                 temp_row_shift_multi_2d = reshape(temp_row_shift_multi,[1,512,512*lY]);
%                 temp_col_shift_multi_2d = reshape(temp_col_shift_multi,[1,512,512*lY]);
%                                 
%                 displacementField = [temp_col_shift_multi_2d; temp_row_shift_multi_2d];
%                 displacementField = permute(displacementField, [2 3 1]);
%                 displacementField = displacementField * -1;


                % old methed, which would lead col interference within batch                
                %shift_multi = [col_shift_multi;row_shift_multi];
                %displacementField = repmat(shift_multi * -1,[1,1,512,512]);                                                
                %displacementField = permute(displacementField, [4 3 2 1]);
                %displacementField = reshape(displacementField,[512,512*lY,2]);                              
                
            % new methed
            col_InfShift = 512*lY*2;

            displacementField_1d_col = col_shift_multi * -1;
            displacementField_1d_col_int = ceil(abs(displacementField_1d_col)).*sign(displacementField_1d_col);
            col_leftBound_multi = max(1,1 - displacementField_1d_col_int);
            col_rightBound_multi = min(512,512-displacementField_1d_col_int);          

            temp_grid_col = 1:512;
            temp_grid_col = repmat(temp_grid_col,[512 1]);
            temp_grid_col2 = reshape(temp_grid_col',[1 1 512 512]);
            col_leftBound_multi2 = col_leftBound_multi';
            col_rightBound_multi2 = col_rightBound_multi';
            temp_bool = temp_grid_col2>=col_leftBound_multi2 & temp_grid_col2<=col_rightBound_multi2;

            displacementField_3d_col = repmat(col_shift_multi' * -1,[1,1,512,512]);
            displacementField_3d_col_valid = displacementField_3d_col.*temp_bool + col_InfShift.*(~temp_bool);
            
            displacementField_3d_row = repmat(row_shift_multi' * -1,[1,1,512,512]);

            displacementField = [displacementField_3d_col_valid,displacementField_3d_row];
            displacementField = permute(displacementField, [4 3 1 2]);
            displacementField = reshape(displacementField,[512,512*lY,2]);

                
                
                if options.if_useExistingTemplate == 1 && options.if_nonRigidExistingTemplate == 1
                    if options.if_imageAverage == 1
                        temp_bin_width = ceil(bin_width/options.imageAverageBin);
                    else
                        temp_bin_width = bin_width;
                    end
                    if lY == temp_bin_width
                        temp_D_batch = D_batch;
                    else
                        temp_D_batch = repmat(D,[1,lY,1]);                        
                    end
                    displacementField = displacementField + temp_D_batch;
                end
                
                if options.if_innerSessionNonRigid == 1
                %if false
                    temp_weight_D_T_to_nk2 = weight_D_T_to_nk2(t,:);
                    temp_locIndex_weight_D_T_to_nk2 = find(temp_weight_D_T_to_nk2>0);
                    if length(temp_locIndex_weight_D_T_to_nk2) == 1
                        %temp_D_T_to_nk2 = D_T_to_nk2{temp_locIndex_weight_D_T_to_nk2};
                        
                        n = temp_locIndex_weight_D_T_to_nk2;
                        temp_range = (2*(n-1)+1):2*n;
                        temp_D_T_to_nk2 = D_T_to_nk2_merged(:,:,temp_range);
                    elseif length(temp_locIndex_weight_D_T_to_nk2) == 2
                        %tempD1 = D_T_to_nk2{temp_locIndex_weight_D_T_to_nk2(1)};
                        %tempD2 = D_T_to_nk2{temp_locIndex_weight_D_T_to_nk2(2)};
                        %temp_D_T_to_nk2 = (tempD1 + tempD2)./2;
                        
                        n = temp_locIndex_weight_D_T_to_nk2(1);
                        temp_range1 = (2*(n-1)+1):2*n;
                        n = temp_locIndex_weight_D_T_to_nk2(2);
                        temp_range2 = (2*(n-1)+1):2*n;
                        
                        tempD1 = D_T_to_nk2_merged(:,:,temp_range1);
                        tempD2 = D_T_to_nk2_merged(:,:,temp_range2);
                        tempW1 = temp_weight_D_T_to_nk2(temp_locIndex_weight_D_T_to_nk2(1));
                        tempW2 = temp_weight_D_T_to_nk2(temp_locIndex_weight_D_T_to_nk2(2));
                        
                        temp_D_T_to_nk2 = tempD1.*tempW1 + tempD2.*tempW2;
                    end
                    
                    temp_D_T_to_nk2_batch = repmat(temp_D_T_to_nk2,[1,lY,1]);
                    %temp_D_T_to_nk2_batch = repmat(gpuArray(temp_D_T_to_nk2),[1,lY,1]);
                    displacementField = displacementField + temp_D_T_to_nk2_batch;
                end
                
                
                % Shift reconstruct
                Mf = reshape(Ytm,[512,512*lY]);                
                %Mf = imwarp(Mf,gather(displacementField),'linear');
                Mf = imwarp(Mf,displacementField,'linear');
                Mf = reshape(Mf,512,512,lY); 
                              
                            %Mf = Mf .* shift_templateBool;
                            Mf = Mf.*shift_templateBool + backgroundColor_template.*(~shift_templateBool);

                            
                            
                            
                            shifts_temp_batch_mat = gather([row_shift_multi; col_shift_multi]);
                            %diff_temp_batch_mat;
                            
                            
                            
                            %for ii = 1:lY
                            %    %output_batch{ii} = gather(dftregistration_min_max_gpu_jjb_v10(temp_fftTemp_mat,temp_fftY_batch(:,:,ii),us_fac,temp_min_shift,temp_max_shift,phase_flag,options.if_FTpad,0));   %#ok<*PFBNS>
                            %    [output_batch{ii},Mf(:,:,ii)] = dftregistration_min_max_gpu_jjb_v14(temp_fftTemp_mat,temp_fftY_batch(:,:,ii),us_fac,temp_min_shift,temp_max_shift,phase_flag,options.if_FTpad,1,Nr2,Nc2,tempNr,tempNc,conj2(:,:,ii),CCreal(:,:,ii));   %#ok<*PFBNS>
                            %    %Mf(:,:,ii) = gather(Greg_real);
                            %end
                            
                        %elseif spmdindex > options.spmdCut_local
                        elseif gpu_globalEnableFlag == 0
                            
                            % Compute cross correlation in frequency domain
                            if phase_flag == 1
                                %conj1 = temp_fftTemp_mat.*conj(temp_fftY_batch);
                                conj2 = temp_fftY_batch.*conj(temp_fftTemp_mat);
                                %conj1_pad = sign(FTpad_gpuBatch_jjb_v5(conj1,lY));
                                %CC = ifft2(conj1_pad);
                                %CCreal = real(CC);
                                CCreal = real(ifft2(sign(FTpad_batch_jjb_v4(temp_fftTemp_mat.*conj(temp_fftY_batch),lY))));
                            elseif phase_flag == 0
                                conj2 = temp_fftY_batch.*conj(temp_fftTemp_mat);
                                CCreal = real(ifft2(FTpad_batch_jjb_v4(temp_fftTemp_mat.*conj(temp_fftY_batch),lY)));
                            end
                            
                            [~, III] = max(CCreal, [], [1 2],'linear');
                            III = reshape(III,1,lY);
                            [row_shift_multi,col_shift_multi,~] = ind2sub([1024 1024 lY],III); %#ok<*ASGLU>
                            
                            row_shift_multi = Nr2(row_shift_multi)/2;
                            col_shift_multi = Nc2(col_shift_multi)/2;
                            
                            row_shift_multi = round(row_shift_multi*us_fac)/us_fac;
                            col_shift_multi = round(col_shift_multi*us_fac)/us_fac;
                            
                            roff_forRefine = 37.5-row_shift_multi*us_fac;
                            coff_forRefine = 37.5-col_shift_multi*us_fac;
                            
                            temp_array = 0:74;
                            roff_forRefine_off = repmat(temp_array',1,lY) - roff_forRefine;
                            roff_forRefine_off = reshape(roff_forRefine_off,1,75,lY);
                            coff_forRefine_off = repmat(temp_array',1,lY) - coff_forRefine;
                            coff_forRefine_off = reshape(coff_forRefine_off,1,75,lY);
                            
                            temp_ifftshift_off = single(ifftshift(0:512-1) - 256);
                            temp_ifftshift_off = repmat(temp_ifftshift_off,1,1,lY);
                            temp_constant = -1i*2*pi/(512*us_fac);
                            temp_ifftshift_off_Xconst = temp_ifftshift_off*temp_constant;
                            
                            % Upsampled DFT by matrix multiplies
                            kernc = exp(pagemtimes(temp_ifftshift_off_Xconst,'transpose',coff_forRefine_off,'none'));
                            kernr = exp(pagemtimes(roff_forRefine_off,'transpose',temp_ifftshift_off_Xconst,'none'));
                            %CC=kernr*conj2*kernc;
                            %CC = pagemtimes(pagemtimes(kernr,conj2),kernc);
                            
                            %CCreal = real(CC);
                            CCreal = real(pagemtimes(pagemtimes(kernr,conj2),kernc));
                            [~, III] = max(CCreal, [], [1 2],'linear');
                            III = reshape(III,1,lY);
                            [rloc_multi,cloc_multi,~] = ind2sub([75 75 lY],III); %#ok<*ASGLU>
                            rloc_multi = rloc_multi - 38.5;
                            cloc_multi = cloc_multi - 38.5;
                            row_shift_multi = row_shift_multi + rloc_multi/us_fac;
                            col_shift_multi = col_shift_multi + cloc_multi/us_fac;
                            
                            %temp_row_shift_multi = reshape(row_shift_multi,1,1,lY);
                            %temp_row_shift_multi = repmat(temp_row_shift_multi,512,512,1);
                            %temp_col_shift_multi = reshape(col_shift_multi,1,1,lY);
                            %temp_col_shift_multi = repmat(temp_col_shift_multi,512,512,1);
                            
                            % Shift reconstruct
                            %Mf = real(ifft2(...
                            %    fftY_batch.*...
                            %    exp(...
                            %    1i*2*pi*(-temp_row_shift_multi.*tempNr_batch-temp_col_shift_multi.*tempNc_batch)...
                            %    )));
                            
                            
%                 temp_row_shift_multi_2d = reshape(temp_row_shift_multi,[1,512,512*lY]);
%                 temp_col_shift_multi_2d = reshape(temp_col_shift_multi,[1,512,512*lY]);
%                                 
%                 displacementField = [temp_col_shift_multi_2d; temp_row_shift_multi_2d];
%                 displacementField = permute(displacementField, [2 3 1]);
%                 displacementField = displacementField * -1;
                

                % old methed, which would lead col interference within batch                
                %shift_multi = [col_shift_multi;row_shift_multi];
                %displacementField = repmat(shift_multi * -1,[1,1,512,512]);                                                
                %displacementField = permute(displacementField, [4 3 2 1]);
                %displacementField = reshape(displacementField,[512,512*lY,2]);                              
                
            % new methed
            col_InfShift = 512*lY*2;

            displacementField_1d_col = col_shift_multi * -1;
            displacementField_1d_col_int = ceil(abs(displacementField_1d_col)).*sign(displacementField_1d_col);
            col_leftBound_multi = max(1,1 - displacementField_1d_col_int);
            col_rightBound_multi = min(512,512-displacementField_1d_col_int);          

            temp_grid_col = 1:512;
            temp_grid_col = repmat(temp_grid_col,[512 1]);
            temp_grid_col2 = reshape(temp_grid_col',[1 1 512 512]);
            col_leftBound_multi2 = col_leftBound_multi';
            col_rightBound_multi2 = col_rightBound_multi';
            temp_bool = temp_grid_col2>=col_leftBound_multi2 & temp_grid_col2<=col_rightBound_multi2;

            displacementField_3d_col = repmat(col_shift_multi' * -1,[1,1,512,512]);
            displacementField_3d_col_valid = displacementField_3d_col.*temp_bool + col_InfShift.*(~temp_bool);
            
            displacementField_3d_row = repmat(row_shift_multi' * -1,[1,1,512,512]);

            displacementField = [displacementField_3d_col_valid,displacementField_3d_row];
            displacementField = permute(displacementField, [4 3 1 2]);
            displacementField = reshape(displacementField,[512,512*lY,2]);

                
                if options.if_useExistingTemplate == 1 && options.if_nonRigidExistingTemplate == 1
                    if options.if_imageAverage == 1
                        temp_bin_width = ceil(bin_width/options.imageAverageBin);
                    else
                        temp_bin_width = bin_width;
                    end
                    if lY == temp_bin_width
                        temp_D_batch = D_batch;
                    else
                        temp_D_batch = repmat(D,[1,lY,1]);                        
                    end
                    displacementField = displacementField + temp_D_batch;
                end
                
                if options.if_innerSessionNonRigid == 1
                    temp_weight_D_T_to_nk2 = weight_D_T_to_nk2(t,:);
                    temp_locIndex_weight_D_T_to_nk2 = find(temp_weight_D_T_to_nk2>0);
                    if length(temp_locIndex_weight_D_T_to_nk2) == 1
                        %temp_D_T_to_nk2 = D_T_to_nk2{temp_locIndex_weight_D_T_to_nk2};
                        
                        n = temp_locIndex_weight_D_T_to_nk2;
                        temp_range = (2*(n-1)+1):2*n;
                        temp_D_T_to_nk2 = D_T_to_nk2_merged(:,:,temp_range);
                    elseif length(temp_locIndex_weight_D_T_to_nk2) == 2
                        %tempD1 = D_T_to_nk2{temp_locIndex_weight_D_T_to_nk2(1)};
                        %tempD2 = D_T_to_nk2{temp_locIndex_weight_D_T_to_nk2(2)};
                        %temp_D_T_to_nk2 = (tempD1 + tempD2)./2;
                        
                        n = temp_locIndex_weight_D_T_to_nk2(1);
                        temp_range1 = (2*(n-1)+1):2*n;
                        n = temp_locIndex_weight_D_T_to_nk2(2);
                        temp_range2 = (2*(n-1)+1):2*n;
                        
                        tempD1 = D_T_to_nk2_merged(:,:,temp_range1);
                        tempD2 = D_T_to_nk2_merged(:,:,temp_range2);
                        tempW1 = temp_weight_D_T_to_nk2(temp_locIndex_weight_D_T_to_nk2(1));
                        tempW2 = temp_weight_D_T_to_nk2(temp_locIndex_weight_D_T_to_nk2(2));
                        
                        temp_D_T_to_nk2 = tempD1.*tempW1 + tempD2.*tempW2;
                    end
                    
                    temp_D_T_to_nk2_batch = repmat(temp_D_T_to_nk2,[1,lY,1]);
                    %temp_D_T_to_nk2_batch = repmat(gpuArray(temp_D_T_to_nk2),[1,lY,1]);
                    displacementField = displacementField + temp_D_T_to_nk2_batch;
                end
                
                % Shift reconstruct
                Mf = reshape(Ytm,[512,512*lY]);                
                %Mf = imwarp(Mf,gather(displacementField),'linear');
                Mf = imwarp(Mf,displacementField,'linear');
                Mf = reshape(Mf,512,512,lY);                            
                            
                            
                            shifts_temp_batch_mat = [row_shift_multi; col_shift_multi];
                            %diff_temp_batch_mat;                            
                        
                        
                            %for ii = 1:lY
                            %    %output_batch{ii} = dftregistration_min_max_jjb_v8(temp_fftTemp_mat,temp_fftY_batch(:,:,ii),us_fac,temp_min_shift,temp_max_shift,phase_flag,options.if_FTpad,0);   %#ok<*PFBNS>
                            %    [output_batch{ii},Mf(:,:,ii)] = dftregistration_min_max_jjb_v11(temp_fftTemp_mat,temp_fftY_batch(:,:,ii),us_fac,temp_min_shift,temp_max_shift,phase_flag,options.if_FTpad,1,Nr2,Nc2,tempNr,tempNc,conj2(:,:,ii),CCreal(:,:,ii));   %#ok<*PFBNS>
                            %    %Mf(:,:,ii) = Greg_real;
                            %end
                        end   
                        
                    else
                        output_batch = cell(1, size_3_Ytm);
                        lb = -max_shift(1,2);
                        ub = max_shift(1,2);
                        max_dev = 0*max_dev_g;


                        %fftTemp_mat = fftTemp{1};
                        temp_min_shift = lb-max_dev(1:2);
                        temp_max_shift = ub+max_dev(1:2);
                        
                        mini_batch_num = options.mini_batch_num;
                        for tempi=1:mini_batch_num
                            temp_mini_batch_width = ceil(lY/mini_batch_num);
                            mini_frame_range = (tempi-1)*temp_mini_batch_width + (1:temp_mini_batch_width);
                            if max(mini_frame_range) > lY
                                mini_frame_range = min(mini_frame_range):lY;
                            end
                            temp_mini_batch_width = length(mini_frame_range);
                            
                            %output_batch(mini_frame_range) = ...
                            %    dftregistration_min_max_batch_jjb_v6(...
                            %    fftTemp{1},fftY_batch(:,:,mini_frame_range),us_fac,lb-max_dev(1:2),ub+max_dev(1:2),phase_flag,temp_mini_batch_width...
                            %    );   %#ok<*PFBNS>
                            if spmdindex <= options.spmdCut_local
                                output_batch(mini_frame_range) = dftregistration_min_max_gpuBatch_jjb_v8(...
                                    temp_fftTemp_mat,temp_fftY_batch(:,:,mini_frame_range),us_fac,temp_min_shift,temp_max_shift,phase_flag,temp_mini_batch_width,Nr2,Nc2);   %#ok<*NODEF,*PFBNS>                                
                            else
                                output_batch(mini_frame_range) = dftregistration_min_max_batch_jjb_v7(...
                                    temp_fftTemp_mat,temp_fftY_batch(:,:,mini_frame_range),us_fac,temp_min_shift,temp_max_shift,phase_flag,temp_mini_batch_width,Nr2,Nc2);   %#ok<*PFBNS>                                                                    
                            end
                            
                            a = 1;
                        end
                        a = 1;
                    end
                    a = 1;
                    
                    
                    %for ii = 1:lY
                    %    %shifts_temp_batch_mat(:,ii) = output_batch{ii}(3:4);
                    %    %diff_temp_batch_mat(ii) = output_batch{ii}(2);
                    %   
                    %    shifts(ii).shifts = shifts_temp_batch_mat(:,ii);
                    %    %shifts(ii).diff = diff_temp_batch_mat(ii);
                    %    shifts(ii).diff = 0;
                    %    shifts(ii).shifts_up = shifts(ii).shifts;
                    %end
                                        
                    % Note: process underMin could cause gpu memory leak!
                    % But when use uin16 translation, any value less than 0 would
                    % be 0 instead. So it's ok to leave underMin alone!
                        %underMin_index = logical(-(sign((Mf-minY_batch)+0.5)-1)/2);
                    %overMax_index = logical((sign((Mf-maxY_batch)-0.5)+1)/2);
                    %if sum(overMax_index==1,'all') > 0
                    %  temp_array = reshape(gpuArray(1:lY),[1 1 lY]);
                    %  temp_mesh = ones(512,512,lY,'single','gpuArray') .* temp_array;
                    %  overMax_batchIndex = temp_mesh(overMax_index);
                    %  Mf(overMax_index) = maxY_batch(overMax_batchIndex);
                    %end
                    
                    %%temp_max = maxY_batch(1);
                    %%temp_max = median(maxY_batch);
                    %temp_max = maxY_batch_median;
                    %overMax_index = logical((sign((Mf-temp_max)-0.5)+1)/2);
                    %%Mf(overMax_index) = temp_max;   
                    %Mf = Mf.*(~overMax_index) + temp_max.*overMax_index;

                    
%             if options.if_useExistingTemplate == 1 && options.if_nonRigidExistingTemplate == 1
%                 Mf = reshape(Mf,[512,512*lY]);                
%                 if options.if_imageAverage == 1
%                     temp_bin_width = ceil(bin_width/options.imageAverageBin);
%                 else
%                     temp_bin_width = bin_width;
%                 end
%                 if lY == temp_bin_width
%                     Mf = imwarp(Mf,D_batch,'linear'); 
%                 else
%                     temp_D_batch = repmat(D,[1,lY,1]);
%                     Mf = imwarp(Mf,temp_D_batch,'linear');
%                 end
%                 Mf = reshape(Mf,[512,512,lY]);                              
%             end     
            
                %overMax_index = Mf > median_maxY_batch;
                %underMin_index = Mf < median_minY_batch;
                %Mf = Mf.*(~overMax_index) + median_maxY_batch.*overMax_index;
                %Mf = Mf.*(~underMin_index) + median_minY_batch.*underMin_index;              
                
                overMax_index = Mf > maxY_batch;
                underMin_index = Mf < minY_batch;
                Mf = Mf.*(~overMax_index) + maxY_batch.*overMax_index;
                Mf = Mf.*(~underMin_index) + minY_batch.*underMin_index;
                    
                                  
                    mem_buffer_cpu = gather(pagetranspose(Mf));
                    
                    if options.if_imageAverage == 1
                        t_write = ceil(t/imageAverageBin);
                        temp_bin_width = ceil(bin_width/imageAverageBin);
                    else
                        t_write = t;
                        temp_bin_width = bin_width;
                    end
                    
                    
                    h5write(temp_h5_filename,['/',options.h5_groupname],...
                        uint16(mem_buffer_cpu),...
                        ...%gather(uint16(pagetranspose(Mf))),...
                        [ones(1,nd),t_write-(tempi_worker(1)-1)*temp_bin_width],...
                        [sizY(1:nd),lY]);
                    
                    h5write(temp_h5_filename_bin,['/',options.h5_groupname],...
                        uint16(mean(mem_buffer_cpu,3)),...
                        ...%gather(uint16(mean(pagetranspose(Mf),3))),...
                        [ones(1,nd),bin_count_forSave-(tempi_worker(1)-1)],...
                        [sizY(1:nd),1]);
                    
                    if options.if_metrics == 1
                        cY(t_write:t_write+lY-1) = motion_metrics_gpuBatch_jjb_v1(Ytm,10);
                        cM(t_write:t_write+lY-1) = motion_metrics_gpuBatch_jjb_v1(Mf,10);
                    end
                    
                end
                fclose(fid);
                time_spmd = toc(t0_spmd);
            end
            a = 1;
            temp_fftTemp_mat;
            if it == iter
                template = cellfun(@(x) x - add_value,template,'un',0);
                template = cell2mat_ov(template,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY);
            end
            if memmap;
                M_final.shifts = shifts_g;
                M_final.template = template;
            end
            spmdindex; %#ok<*VUNUS>
            %M_final = M_final{1};
            if options.if_metrics == 1
                if options.if_imageAverage == 1
                    temp_T = ceil(T/options.imageAverageBin);
                    temp_bin_width = ceil(bin_width/options.imageAverageBin);
                else
                    temp_T = T;
                    temp_bin_width = bin_width;
                end
                
                cY_all = [];
                cM_all = [];
                for tempi=1:length(spmdindex)
                    temp_cY = cY{tempi};
                    temp_cM = cM{tempi};
                    if tempi == 1
                        temp_range = 1:i_worker_end(tempi)*temp_bin_width;
                    else
                        temp_range = (i_worker_end(tempi-1))*temp_bin_width+1:i_worker_end(tempi)*temp_bin_width;
                        %a = 1;
                    end
                    temp_range(temp_range>temp_T) = [];
                    cY_all = [cY_all temp_cY(temp_range)];
                    cM_all = [cM_all temp_cM(temp_range)];
                end
                cY = cY_all;
                cM = cM_all;
            else
                cY = 0;
                cM = 0;
            end
            
            %shifts_g_all = [];
            time_spmd_all = [];
            for tempi=1:length(spmdindex)
                %temp_shifts_g = shifts_g{tempi};

                if tempi == 1
                    temp_range = 1:i_worker_end(tempi)*bin_width;
                else
                    temp_range = (i_worker_end(tempi-1))*bin_width+1:i_worker_end(tempi)*bin_width;
                    %a = 1;
                end
                temp_range(temp_range>T) = [];
                %shifts_g_all = [shifts_g_all; temp_shifts_g(temp_range)]; %#ok<*AGROW>
                time_spmd_all = [time_spmd_all time_spmd{tempi}];
            end
            %shifts_g = shifts_g_all;
            a = 1;
        end
        
        home;% To scroll down in command window
        str=['All ', num2str(T), ' frames registered.'];
        fprintf("%s\n",str);
        
        fprintf("time_spmd=");
        for tempi=1:workerNum
            fprintf('%.1f ',time_spmd_all(tempi));
        end
        %fprintf("\n");
        %options = options{1};
        options.time_spmd = time_spmd_all;
        
    elseif options.if_spmd_singleGPU == 1
        % abandoned now
    end
    
    
end

if if_spmd == 0
% elseif if_spmd == 0
    for it = 1:iter
        %for t = 1:bin_width:T
        
        fid = fopen(Y);
        imsize = 512*512*2;
        current_seek = ftell(fid);
        fseek(fid, 0, 1);
        file_length = ftell(fid);
        fseek(fid, current_seek, -1);
        frame_length = file_length/imsize;
        bstring = 'uint16=>single';
                        
        imageAverageBin = options.imageAverageBin;
        
        for bin_count = 1:bin_totalNum
            t = (bin_count-1)*bin_width + 1;
            
            if it == iter
                bin_count_forSave = bin_count;
            end
            % Little slow
            %switch filetype
            %    case 'tif'
            %        Ytm = single(read_file(Y, t, min(t+bin_width-1,T)-t+1, [], tiffInfo));
            %    case 'hdf5'
            %        Ytm = single(h5read(Y,data_name,[ones(1,nd),t],[sizY(1:nd),min(t+bin_width-1,T)-t+1]));
            %    case 'mem'
            %        if nd == 2; Ytm = single(Y.(var_name)(:,:,t:min(t+bin_width-1,T))); end
            %        if nd == 3; Ytm = single(Y.(var_name)(:,:,:,t:min(t+bin_width-1,T))); end
            %    case 'mat'
            %        if nd == 2; Ytm = single(Y(:,:,t:min(t+bin_width-1,T))); end
            %        if nd == 3; Ytm = single(Y(:,:,:,t:min(t+bin_width-1,T))); end
            %    case 'raw'
            %        Ytm = gpuArray(single(read_raw_file(Y,t,min(t+bin_width-1,T)-t+1,FOV,bitsize)));
            %end
            
            % Read Raw file
            %Ytm = single(read_raw_file(Y,t,min(t+bin_width-1,T)-t+1,[512 512],2));
            window = min(t+bin_width-1,T)-t+1;            
            feek_indicator = (t-1)*imsize;
            fseek(fid,feek_indicator,'bof');
            Ytm = gpuArray(fread(fid,512*512*window,'uint16=>single',0,'l')');
            Ytm = Ytm/2;% Preprocess for suite2p bin (uint16-->int16, although not now)
            Ytm = reshape(Ytm,[512 512 window]);
            Ytm = pagetranspose(Ytm);
            size_3_Ytm = size(Ytm,ndims(Ytm));
            
            
            
%             % Get a bool mask to denoise images.
%             Ytm = reshape(Ytm,[512,512*size_3_Ytm]);
%             %BW_raw = gather(Ytm > 200);
%             BW_raw = gather(Ytm > 100);
%             BW = bwpack(BW_raw);
%             
%             nhood = [1 1];
%             BW1 = imerode(BW,nhood,'ispacked',512);
%             BW1_2 = imdilate(BW1,nhood,'ispacked');
%             nhood = [1; 1];
%             BW2 = imerode(BW,nhood,'ispacked',512);
%             BW2_2 = imdilate(BW2,nhood,'ispacked');
%             nhood = [0 1; 1 0];
%             BW3 = imerode(BW,nhood,'ispacked',512);
%             BW3_2 = imdilate(BW3,nhood,'ispacked');
%             nhood = [1 0; 0 1];
%             BW4 = imerode(BW,nhood,'ispacked',512);
%             BW4_2 = imdilate(BW4,nhood,'ispacked');
%             
%             closeBW = bitor(bitor(bitor(BW1_2,BW2_2),BW3_2),BW4_2);
%             %closeBW = bitor(BW1_2,BW2_2);
%             se = strel('rectangle',[3 3]);
%             closeBW_dilate = bwunpack(imdilate(closeBW,se,'ispacked'),512);
%             
%             %Ytm(~closeBW_dilate) = backgroundColor_template;
%             Ytm = Ytm.*closeBW_dilate + backgroundColor_template.*(~closeBW_dilate);
%             Ytm = reshape(Ytm,[512,512,size_3_Ytm]);
            
            
            lY_raw = size_3_Ytm;
            if options.if_imageAverage == 1                
                if ceil(size_3_Ytm/imageAverageBin)*imageAverageBin == size_3_Ytm
                    Ytm_mean = squeeze(mean(reshape(Ytm,[512,512,imageAverageBin,size_3_Ytm/imageAverageBin]), 3));                    
                else
                    Ytm_mean = squeeze(mean(reshape(...
                        Ytm(:,:,1:floor(size_3_Ytm/imageAverageBin)*imageAverageBin),[512,512,imageAverageBin,floor(size_3_Ytm/imageAverageBin)]), 3));
                    tempi = (floor(size_3_Ytm/imageAverageBin)*imageAverageBin+1):size_3_Ytm;
                    Ytm_mean(:,:,(floor(size_3_Ytm/imageAverageBin)+1)) = mean(Ytm(:,:,tempi), 3);                                          
                end

                Ytm = Ytm_mean;                
                size_3_Ytm = size(Ytm,ndims(Ytm));
                Ytm_mean = [];
            end                        
            
            Mf = zeros(size(Ytm), 'single', 'gpuArray');
            lY = size_3_Ytm;
            shifts = struct('shifts',cell(lY,1),'shifts_up',cell(lY,1),'diff',cell(lY,1));                                                
            shifts_temp_batch_mat = zeros(2, size_3_Ytm, 'single', 'gpuArray');    
            
            [minY_batch,maxY_batch] = bounds(Ytm,[1 2]);
            median_minY_batch = median(minY_batch,'all');
            median_maxY_batch = median(maxY_batch,'all');
            fftY_batch = fft2(Ytm);
            lb = -max_shift(1,nd);
            ub = max_shift(1,nd);
            max_dev = 0*max_dev_g;
            
            %% dftregistration_min_max
            if options.if_dftreg_bacth == 0
                %fftTemp_mat = fftTemp{1};
                temp_min_shift = lb-max_dev(1:2);
                temp_max_shift = ub+max_dev(1:2);
                
                Nr2 = ifftshift(gpuArray(-512:511));
                Nc2 = ifftshift(gpuArray(-512:511));
                tempNr = gpuArray(Nr{1})/512;
                tempNc = gpuArray(Nc{1})/512;
                
                tempNr_batch = repmat(tempNr,[1 1 lY]);
                tempNc_batch = repmat(tempNc,[1 1 lY]);
                
                % Compute cross correlation in frequency domain
                if phase_flag == 1
                    %conj1 = fftTemp_mat.*conj(fftY_batch);
                    conj2 = fftY_batch.*conj(fftTemp_mat);
                    %conj1_pad = sign(FTpad_gpuBatch_jjb_v5(conj1,lY));
                    %CC = ifft2(conj1_pad);
                    %CCreal = real(CC);
                    CCreal = real(ifft2(sign(FTpad_gpuBatch_jjb_v5(fftTemp_mat.*conj(fftY_batch),lY))));
                elseif phase_flag == 0
                    conj2 = fftY_batch.*conj(fftTemp_mat);
                    CCreal = real(ifft2(FTpad_gpuBatch_jjb_v5(fftTemp_mat.*conj(fftY_batch),lY)));
                end
                
                [~, III] = max(CCreal, [], [1 2],'linear');
                III = reshape(III,1,lY);
                [row_shift_multi,col_shift_multi,~] = ind2sub([1024 1024 lY],III); %#ok<*ASGLU>
                
                row_shift_multi = Nr2(row_shift_multi)/2;
                col_shift_multi = Nc2(col_shift_multi)/2;
                
                 % Upsampling, refine estimate with matrix multiply DFT
                row_shift_multi = round(row_shift_multi*us_fac)/us_fac;
                col_shift_multi = round(col_shift_multi*us_fac)/us_fac;
                
                roff_forRefine = 37.5-row_shift_multi*us_fac;
                coff_forRefine = 37.5-col_shift_multi*us_fac;
                
                temp_array = 0:74;
                roff_forRefine_off = repmat(temp_array',1,lY) - roff_forRefine;
                roff_forRefine_off = reshape(roff_forRefine_off,1,75,lY);
                coff_forRefine_off = repmat(temp_array',1,lY) - coff_forRefine;
                coff_forRefine_off = reshape(coff_forRefine_off,1,75,lY);
                
                temp_ifftshift_off = gpuArray(single(ifftshift(0:512-1) - 256));
                temp_ifftshift_off = repmat(temp_ifftshift_off,1,1,lY);
                temp_constant = gpuArray(-1i*2*pi/(512*us_fac));
                temp_ifftshift_off_Xconst = temp_ifftshift_off*temp_constant;
                
                % Upsampled DFT by matrix multiplies
                kernc = exp(pagemtimes(temp_ifftshift_off_Xconst,'transpose',coff_forRefine_off,'none'));
                kernr = exp(pagemtimes(roff_forRefine_off,'transpose',temp_ifftshift_off_Xconst,'none'));
                %CC=kernr*conj2*kernc;
                %CC = pagemtimes(pagemtimes(kernr,conj2),kernc);
                
                %CCreal = real(CC);
                CCreal = real(pagemtimes(pagemtimes(kernr,conj2),kernc));
                [~, III] = max(CCreal, [], [1 2],'linear');
                III = reshape(III,1,lY);
                [rloc_multi,cloc_multi,~] = ind2sub([75 75 lY],III); %#ok<*ASGLU>
                rloc_multi = rloc_multi - 38.5;
                cloc_multi = cloc_multi - 38.5;
                row_shift_multi = row_shift_multi + rloc_multi/us_fac;
                col_shift_multi = col_shift_multi + cloc_multi/us_fac;
                
                %temp_row_shift_multi = reshape(row_shift_multi,1,1,lY);
                %temp_row_shift_multi = repmat(temp_row_shift_multi,512,512,1);
                %temp_col_shift_multi = reshape(col_shift_multi,1,1,lY);
                %temp_col_shift_multi = repmat(temp_col_shift_multi,512,512,1);
                
                % Shift reconstruct
                %Mf = gather(real(ifft2(...
                %    fftY_batch.*...
                %    exp(...
                %        1i*2*pi*(-temp_row_shift_multi.*tempNr_batch-temp_col_shift_multi.*tempNc_batch)...
                %    ))));
                
                % Shift reconstruct
%                 Mf = real(ifft2(...
%                    fftY_batch.*...
%                    exp(...
%                    1i*2*pi*(-temp_row_shift_multi.*tempNr_batch-temp_col_shift_multi.*tempNc_batch)...
%                    )));                
                
                % old methed, which would lead col interference within batch                
                %shift_multi = [col_shift_multi;row_shift_multi];
                %displacementField = repmat(shift_multi * -1,[1,1,512,512]);                                                
                %displacementField = permute(displacementField, [4 3 2 1]);
                %displacementField = reshape(displacementField,[512,512*lY,2]);                              
                
            % new methed
            col_InfShift = 512*lY*2;

            displacementField_1d_col = col_shift_multi * -1;
            displacementField_1d_col_int = ceil(abs(displacementField_1d_col)).*sign(displacementField_1d_col);
            col_leftBound_multi = max(1,1 - displacementField_1d_col_int);
            col_rightBound_multi = min(512,512-displacementField_1d_col_int);          

            temp_grid_col = 1:512;
            temp_grid_col = repmat(temp_grid_col,[512 1]);
            temp_grid_col2 = reshape(temp_grid_col',[1 1 512 512]);
            col_leftBound_multi2 = col_leftBound_multi';
            col_rightBound_multi2 = col_rightBound_multi';
            temp_bool = temp_grid_col2>=col_leftBound_multi2 & temp_grid_col2<=col_rightBound_multi2;

            displacementField_3d_col = repmat(col_shift_multi' * -1,[1,1,512,512]);
            displacementField_3d_col_valid = displacementField_3d_col.*temp_bool + col_InfShift.*(~temp_bool);
            
            displacementField_3d_row = repmat(row_shift_multi' * -1,[1,1,512,512]);

            displacementField = [displacementField_3d_col_valid,displacementField_3d_row];
            displacementField = permute(displacementField, [4 3 1 2]);
            displacementField = reshape(displacementField,[512,512*lY,2]);
                
                
                if options.if_useExistingTemplate == 1 && options.if_nonRigidExistingTemplate == 1
                    if options.if_imageAverage == 1
                        temp_bin_width = ceil(bin_width/options.imageAverageBin);
                    else
                        temp_bin_width = bin_width;
                    end
                    if lY == temp_bin_width
                        temp_D_batch = D_batch;
                    else
                        temp_D_batch = repmat(D,[1,lY,1]);                        
                    end
                    displacementField = displacementField + temp_D_batch;
                end                
                
                
                if options.if_innerSessionNonRigid == 1
                    temp_weight_D_T_to_nk2 = weight_D_T_to_nk2(t,:);
                    temp_locIndex_weight_D_T_to_nk2 = find(temp_weight_D_T_to_nk2>0);
                    if length(temp_locIndex_weight_D_T_to_nk2) == 1
                        %temp_D_T_to_nk2 = D_T_to_nk2{temp_locIndex_weight_D_T_to_nk2};
                        
                        n = temp_locIndex_weight_D_T_to_nk2;
                        temp_range = (2*(n-1)+1):2*n;
                        temp_D_T_to_nk2 = D_T_to_nk2_merged(:,:,temp_range);
                    elseif length(temp_locIndex_weight_D_T_to_nk2) == 2
                        %tempD1 = D_T_to_nk2{temp_locIndex_weight_D_T_to_nk2(1)};
                        %tempD2 = D_T_to_nk2{temp_locIndex_weight_D_T_to_nk2(2)};
                        %temp_D_T_to_nk2 = (tempD1 + tempD2)./2;
                        
                        n = temp_locIndex_weight_D_T_to_nk2(1);
                        temp_range1 = (2*(n-1)+1):2*n;
                        n = temp_locIndex_weight_D_T_to_nk2(2);
                        temp_range2 = (2*(n-1)+1):2*n;
                        
                        tempD1 = D_T_to_nk2_merged(:,:,temp_range1);
                        tempD2 = D_T_to_nk2_merged(:,:,temp_range2);
                        tempW1 = temp_weight_D_T_to_nk2(temp_locIndex_weight_D_T_to_nk2(1));
                        tempW2 = temp_weight_D_T_to_nk2(temp_locIndex_weight_D_T_to_nk2(2));
                        
                        temp_D_T_to_nk2 = tempD1.*tempW1 + tempD2.*tempW2;
                    end
                    
                    temp_D_T_to_nk2_batch = repmat(temp_D_T_to_nk2,[1,lY,1]);
                    %temp_D_T_to_nk2_batch = repmat(gpuArray(temp_D_T_to_nk2),[1,lY,1]);
                    displacementField = displacementField + temp_D_T_to_nk2_batch;
                end
                
                
                % Shift reconstruct
                Mf = reshape(Ytm,[512,512*lY]);                
                %Mf = imwarp(Mf,gather(displacementField),'linear');
                Mf = imwarp(Mf,displacementField,'linear');
                Mf = reshape(Mf,512,512,lY);
                
                a = 1;
                
                %for ii=1:lY
                %   Mf(:,:,ii) = remove_boundaries_jjb_v2(Mf(:,:,ii),[row_shift_multi(ii),col_shift_multi(ii)]);
                %end
                
%                 row_shift_multi_round = round(row_shift_multi);
%                 col_shift_multi_round = round(col_shift_multi);
%                 
%                 row_shift_multi_round_batch = true(512,512,lY,'gpuArray');
%                 col_shift_multi_round_batch = true(512,512,lY,'gpuArray');
%                 
%                 temp_index_row_shift = [];
%                 for ii=1:lY
%                     if row_shift_multi_round(ii) == 0
%                         continue
%                     end
%                     if row_shift_multi_round(ii) < 0
%                         temp_array_start = 512+row_shift_multi_round(ii)+1;
%                         temp1 = gpuArray(1:512:(511*512+1)) + (temp_array_start-1);
%                         temp2 = repmat(temp1,-row_shift_multi_round(ii),1);
%                         temp3 = temp2 + ((1:-row_shift_multi_round(ii))-1)';
%                         temp4 = reshape(temp3,1,numel(temp3));
%                         temp5 = temp4 + (ii-1)*512*512;
%                         temp_index_row_shift = [temp_index_row_shift temp5];
%                     elseif row_shift_multi_round(ii) > 0
%                         temp_array_end = row_shift_multi_round(ii);
%                         temp1 = gpuArray(1:512:(511*512+1));
%                         temp2 = repmat(temp1,temp_array_end,1);
%                         temp3 = temp2 + ((1:temp_array_end)-1)';
%                         temp4 = reshape(temp3,1,numel(temp3));
%                         temp5 = temp4 + (ii-1)*512*512;
%                         temp_index_row_shift = [temp_index_row_shift temp5];
%                     end
%                 end
%                 row_shift_multi_round_batch(temp_index_row_shift) = false;
%                 
%                 temp_index_col_shift = [];
%                 for ii=1:lY
%                     if col_shift_multi_round(ii) == 0
%                         continue
%                     end
%                     if col_shift_multi_round(ii) < 0
%                         temp_array_start = 512+col_shift_multi_round(ii)+1;
%                         temp1 = ((temp_array_start-1)*512+1):512*512;
%                         temp2 = temp1 + (ii-1)*512*512;
%                         temp_index_col_shift = [temp_index_col_shift temp2];
%                     elseif col_shift_multi_round(ii) > 0
%                         temp_array_end = col_shift_multi_round(ii);
%                         temp1 = 1:temp_array_end*512;
%                         temp2 = temp1 + (ii-1)*512*512;
%                         temp_index_col_shift = [temp_index_col_shift temp2];
%                     end
%                 end
%                 col_shift_multi_round_batch(temp_index_col_shift) = false;
%                 
%                 
%                 row_shift_multi_round_batch;
%                 col_shift_multi_round_batch;
%                 Mf = Mf .* row_shift_multi_round_batch .* col_shift_multi_round_batch;
                

                %shift_templateBool;
                %Mf = Mf .* shift_templateBool;                             
                Mf = Mf.*shift_templateBool + backgroundColor_template.*(~shift_templateBool);

                shifts_temp_batch_mat = [row_shift_multi; col_shift_multi];
                %diff_temp_batch_mat;
               
            else
                output_batch = cell(1, size_3_Ytm);
                
                mini_batch_num = options.mini_batch_num;
                for mini_batch_count=1:mini_batch_num
                    temp_mini_batch_width = ceil(lY/mini_batch_num);
                    mini_frame_range = (mini_batch_count-1)*temp_mini_batch_width + (1:temp_mini_batch_width);
                    if max(mini_frame_range) > lY
                        mini_frame_range = min(mini_frame_range):lY;
                    end
                    temp_mini_batch_width = length(mini_frame_range);
                    
                    %output_batch(mini_frame_range) = ...
                    %   dftregistration_min_max_batch_jjb_v6(...
                    %   fftTemp{1},fftY_batch(:,:,mini_frame_range),us_fac,lb-max_dev(1:2),ub+max_dev(1:2),phase_flag,temp_mini_batch_width...
                    %   );   %#ok<*PFBNS>
                    
                    output_batch(mini_frame_range) = ...
                        dftregistration_min_max_gpuBatch_jjb_v7(...
                        fftTemp_mat,fftY_batch(:,:,mini_frame_range),us_fac,lb-max_dev(1:2),ub+max_dev(1:2),phase_flag,temp_mini_batch_width...
                        );   %#ok<*PFBNS>
                    
                    a = 1;
                end
                
                %output_batch = dftregistration_min_max_batch_jjb_v6(fftTemp{1},fftY_batch,us_fac,lb-max_dev(1:2),ub+max_dev(1:2),phase_flag,lY);   %#ok<*PFBNS>
                a = 1;
            end
            a = 1;
            %[tform,peakcorr] = imregcorr(gather(Ytm(:,:,ii)),gather(template{1}));
            %B = imwarp(A,tform)
            %parfor ii = 1:lY
            for ii = 1:lY
                %shifts_temp_batch_mat(:,ii) = output_batch{ii}(3:4);
                %diff_temp_batch_mat(ii) = output_batch{ii}(2);
                
                shifts(ii).shifts = shifts_temp_batch_mat(:,ii);
                %shifts(ii).diff = diff_temp_batch_mat(ii);
                shifts(ii).diff = 0;
                shifts(ii).shifts_up = shifts(ii).shifts;
            end
            
            %             for ii = 1:lY
            %                temp_Mf = Mf(:,:,ii);
            %                temp_Mf(temp_Mf<minY_batch(ii))=minY_batch(ii);
            %                temp_Mf(temp_Mf>maxY_batch(ii))=maxY_batch(ii);
            %                Mf(:,:,ii) = temp_Mf;
            %             end
            
            % Note: process underMin could cause gpu memory leak!
            % But when use uin16 translation, any value less than 0 would
            % be 0 instead. So it's ok to leave underMin alone!
            %underMin_index = logical(-(sign((Mf-minY_batch)+0.5)-1)/2);
            %overMax_index = logical((sign((Mf-maxY_batch)-0.5)+1)/2);
            %if sum(overMax_index==1,'all') > 0
            %    temp_array = reshape(gpuArray(1:lY),[1 1 lY]);
            %    temp_mesh = ones(512,512,lY,'single','gpuArray') .* temp_array;
            %    overMax_batchIndex = temp_mesh(overMax_index);
            %    Mf(overMax_index) = maxY_batch(overMax_batchIndex);
            %end
            
            
            %             if options.if_imageAverage == 1
            %                 if ceil(size_3_Ytm/imageAverageBin)*imageAverageBin == size_3_Ytm
            %                     Mf_mean = squeeze(mean(reshape(Mf,[512,512,imageAverageBin,size_3_Ytm/imageAverageBin]), 3));
            %                 else
            %                     Mf_mean = squeeze(mean(reshape(...
            %                         Mf(:,:,1:floor(size_3_Ytm/imageAverageBin)*imageAverageBin),[512,512,imageAverageBin,floor(size_3_Ytm/imageAverageBin)]), 3));
            %                     tempi = (floor(size_3_Ytm/imageAverageBin)*imageAverageBin+1):size_3_Ytm;
            %                     Mf_mean(:,:,(floor(size_3_Ytm/imageAverageBin)+1)) = mean(Mf(:,:,tempi), 3);
            %                 end
            %
            %                 Mf = Mf_mean;
            %                 Mf_mean = [];
            %                 lY = size(Mf,ndims(Mf));
            %             end
            
            
            
%             if options.if_useExistingTemplate == 1 && options.if_nonRigidExistingTemplate == 1
%                 %for tempi=1:lY
%                 %    Mf(:,:,tempi) = imwarp(Mf(:,:,tempi),D,'linear'); 
%                 %end       
%                 
%                 Mf = reshape(Mf,[512,512*lY]);                
%                 if options.if_imageAverage == 1
%                     temp_bin_width = ceil(bin_width/options.imageAverageBin);
%                 else
%                     temp_bin_width = bin_width;
%                 end
%                 if lY == temp_bin_width
%                     Mf = imwarp(Mf,D_batch,'linear'); 
%                 else
%                     temp_D_batch = repmat(D,[1,lY,1]);
%                     Mf = imwarp(Mf,temp_D_batch,'linear');
%                 end
%                 Mf = reshape(Mf,[512,512,lY]);
%             end
                            
            %overMax_index = Mf > median_maxY_batch;
            %underMin_index = Mf < median_minY_batch;
            %Mf = Mf.*(~overMax_index) + median_maxY_batch.*overMax_index;
            %Mf = Mf.*(~underMin_index) + median_minY_batch.*underMin_index;            
            
            overMax_index = Mf > maxY_batch;
            underMin_index = Mf < minY_batch;
            Mf = Mf.*(~overMax_index) + maxY_batch.*overMax_index;
            Mf = Mf.*(~underMin_index) + minY_batch.*underMin_index;
            
            mem_buffer_cpu = gather(pagetranspose(Mf));
            
            if options.if_imageAverage == 1
                %options.imageAverageBin;
                %t_write = t;
                t_write = ceil(t/options.imageAverageBin);
            else
                t_write = t;
            end
            
            h5write(options.h5_filename,['/',options.h5_groupname],...
                uint16(mem_buffer_cpu),...
                [ones(1,nd),t_write],...
                [sizY(1:nd),lY]);
            
            h5write(options.h5_filename_bin,['/',options.h5_groupname],...
                uint16(mean(mem_buffer_cpu,3)),...
                [ones(1,nd),bin_count_forSave],...
                [sizY(1:nd),1]);
            
            if options.if_metrics == 1
                cY(t_write:t_write+lY-1) = motion_metrics_gpuBatch_jjb_v1(Ytm,10);
                %cY(t_write:t_write+lY-1) = motion_metrics_gpuBatch_jjb_v1(Ytm_mean,10);
                cM(t_write:t_write+lY-1) = motion_metrics_gpuBatch_jjb_v1(Mf,10);
            end
            
            
            str=[num2str(t+lY_raw-1), ' out of ', num2str(T), ' frames registered, iteration ', num2str(it), ' out of ', num2str(iter), '..'];
            refreshdisp(str, prevstr, t);
            prevstr=str;
            
            
        end
        fclose(fid);
        
        if it == iter
            template = cellfun(@(x) x - add_value,template,'un',0);
            template = cell2mat_ov(template,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY);
        end
        if memmap;
            M_final.shifts = shifts_g;
            M_final.template = template;
        end
        %maxNumCompThreads('automatic');       
    end    
    
    if options.if_metrics == 0
       cY = 0;
       cM = 0;
    end
end

if options.if_useExistingTemplate == 1 && options.if_nonRigidExistingTemplate == 1
    template = imwarp(template,gather(D),'cubic');%'nearest' would lose sub-pixel displayment!
    overMax_index = template > max_template;
    underMin_index = template < min_template;
    template = template.*(~overMax_index) + max_template.*overMax_index;
    template = template.*(~underMin_index) + min_template.*underMin_index;
end

