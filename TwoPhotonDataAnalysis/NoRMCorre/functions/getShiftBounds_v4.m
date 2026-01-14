function shift_templateBool = getShiftBounds_v4(Y,init_batch_shift,T,options)
% init_batch_shift = 2000;
% T;
% options;
% Y;

t0_fun = tic;

init_batch = min(T,init_batch_shift);

us_fac = options.us_fac;

% interleaved
init_Y_index = 1:floor(T/init_batch):T;
init_Y_index = init_Y_index(1:init_batch);
temp_init_batch = init_batch;

Y_temp = zeros(options.d1,options.d2,init_batch,'single');

fid = fopen(Y);
imsize = 512*512*2;
current_seek = ftell(fid);
fseek(fid, 0, 1);
% file_length = ftell(fid);
fseek(fid, current_seek, -1);
% frame_length = file_length/imsize;
% bstring = 'uint16=>single';
for tempi=1:temp_init_batch
    t = init_Y_index(tempi);
    feek_indicator = (t-1)*imsize;
    fseek(fid,feek_indicator,'bof');
    Y_temp(:,:,tempi) = fread(fid,[512,512],'uint16=>single',0,'l')';
end
fclose(fid);

size_3_Y_temp = size(Y_temp,3);

% data_type = 'uint16';

% fprintf('Registering the interleaved %i frames just to obtain a good template...',init_batch);

if size_3_Y_temp < 500 % To avoid memory overflow
    template_in = median(Y_temp,3);
else
    %template_in = median(Y_temp(:,:,1:500),3);
    
    temp_Y_temp_interleaved = 1:ceil(init_batch/500):init_batch;
    template_in = median(Y_temp(:,:,temp_Y_temp_interleaved),3);    
end

phase_flag = options.phase_flag;

row_shift_template = zeros(1,size_3_Y_temp,'single','gpuArray');
col_shift_template = zeros(1,size_3_Y_temp,'single','gpuArray');

bin_width_forTemplate = 22;%8-->15-->22
bin_totalNum_forTemplate = length(1:bin_width_forTemplate:size_3_Y_temp);

fftTemp_mat = fftn(template_in);
for bin_count = 1:bin_totalNum_forTemplate
    t = (bin_count-1)*bin_width_forTemplate + 1;
    size_3_Y_temp_batch = min(bin_width_forTemplate,size_3_Y_temp-t+1);
    
    Y_temp_batch = gpuArray(Y_temp(:,:,t:(t+size_3_Y_temp_batch-1)));
    lY = size_3_Y_temp_batch;
    
    [minY_batch,maxY_batch] = bounds(Y_temp_batch,[1 2]);
    %median_minY_batch = median(minY_batch,'all');
    %median_maxY_batch = median(maxY_batch,'all');
    fftY_batch = fft2(Y_temp_batch);
    
    %% dftregistration_min_max
    temp_min_shift = [-20 -20];%20
    temp_max_shift = [20 20];%20
    
    Nr2 = ifftshift(gpuArray(-512:511));
    Nc2 = ifftshift(gpuArray(-512:511));
    %         tempNr = gpuArray(Nr{1})/512;
    %         tempNc = gpuArray(Nc{1})/512;
    
    %         tempNr_batch = repmat(tempNr,[1 1 lY]);
    %         tempNc_batch = repmat(tempNc,[1 1 lY]);
    
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
    
    %         temp_row_shift_multi = reshape(row_shift_multi,1,1,lY);
    %         temp_row_shift_multi = repmat(temp_row_shift_multi,512,512,1);
    %         temp_col_shift_multi = reshape(col_shift_multi,1,1,lY);
    %         temp_col_shift_multi = repmat(temp_col_shift_multi,512,512,1);
    
    temp_index = t:(t+size_3_Y_temp_batch-1);
    temp_range = 1:size_3_Y_temp;
    temp_ismember = ismember(temp_range,temp_index);
    for tempi=1:lY
        row_shift_template = row_shift_template.*(~temp_ismember) + row_shift_multi(tempi).*temp_ismember;
        col_shift_template = col_shift_template.*(~temp_ismember) + col_shift_multi(tempi).*temp_ismember;
    end
    
end


row_shift_template = gather(row_shift_template);
col_shift_template = gather(col_shift_template);

clear CCreal CCreal2 cloc_multi coff_forRefine coff_forRefine_off col_shift_multi conj2
clear fftY_batch III kernc kernr maxY_batch
clear median_maxY_batch median_minY_batch Mf
clear minY_batch Nc2 Nr2 rloc_multi roff_forRefine roff_forRefine_off row_shift_multi
clear temp1 temp1234_multi temp13 temp13_multi temp2 temp24 temp24_multi temp3 temp4
clear temp_boolArray temp_boolIndex temp_col_shift_multi temp_constant temp_ifftshift_off temp_ifftshift_off_Xconst
clear temp_row_shift_multi tempNc tempNc_batch tempNr tempNr_batch Y_temp_batch

% fprintf('...done, time=%.2f seconds.\n',toc(t0_fun));

% To get shiftBounds_template
outlierIndex_row_shift_template = isoutlier(row_shift_template,'grubbs');
outlierIndex_col_shift_template = isoutlier(col_shift_template,'grubbs');
row_shift_template_valid = row_shift_template(~outlierIndex_row_shift_template);
col_shift_template_valid = col_shift_template(~outlierIndex_col_shift_template);

[min_row_shift_template, max_row_shift_template] = bounds(row_shift_template_valid);
[min_col_shift_template, max_col_shift_template] = bounds(col_shift_template_valid);
shiftBounds_template = [min_row_shift_template,max_row_shift_template,...
    min_col_shift_template,max_col_shift_template];

shiftBounds_template = shiftBounds_template ./ 2; % add this line in 20230718, to loose the criteria

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

fprintf('Get shift bounds from interleaved %d frames, time=%.2f seconds.\n',init_batch,toc(t0_fun));
