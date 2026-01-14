function [D_T_to_nk2, template_nk2, template_nk2_corrected] = ...
    innerSessionNonRigid_v11(template_T,Y,T,options,N,k1,k2,temp_options)


a = 1;

if_plot = 0;

t0_fun = tic;

Nr = temp_options.Nr;
Nc = temp_options.Nc;


T_n = round(linspace(1,T,N+1));

k1 = min(T_n(2),k1);

us_fac = options.us_fac;

template_nk2 = zeros(options.d1,options.d2,N,'single');
fprintf('Get template_nk2, ');
for n=1:N
    % interleaved
    Y_index_nk1 = round(linspace(T_n(n),T_n(n+1),k1));
    
    Y_temp = zeros(options.d1,options.d2,k1,'single');
    
    fid = fopen(Y);
    imsize = 512*512*2;
    current_seek = ftell(fid);
    fseek(fid, 0, 1);
    fseek(fid, current_seek, -1);
    for tempi=1:k1
        t = Y_index_nk1(tempi);
        feek_indicator = (t-1)*imsize;
        fseek(fid,feek_indicator,'bof');
        Y_temp(:,:,tempi) = fread(fid,[512,512],'uint16=>single',0,'l')';
    end
    fclose(fid);
    
    size_3_Y_temp = size(Y_temp,3);
    
    
    phase_flag = options.phase_flag;
    
    bin_width_forTemplate = 21;%8-->15-->22-->21
    bin_totalNum_forTemplate = length(1:bin_width_forTemplate:size_3_Y_temp);
    a = 1;
    
    %Y_temp_corrected = Y_temp;
    Y_temp_corrected = Y_temp;
    % profile viewer
    % t0_template = tic;
    
    fftTemp_mat = fftn(template_T);
    
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
        
        %temp_row_shift_multi = reshape(row_shift_multi,1,1,lY);
        %temp_row_shift_multi = repmat(temp_row_shift_multi,512,512,1);
        %temp_col_shift_multi = reshape(col_shift_multi,1,1,lY);
        %temp_col_shift_multi = repmat(temp_col_shift_multi,512,512,1);
                
        % % Shift reconstruct
        %Mf = real(ifft2(...
        %   fftY_batch.*...
        %   exp(...
        %   1i*2*pi*(-temp_row_shift_multi.*tempNr_batch-temp_col_shift_multi.*tempNc_batch)...
        %   )));
        
    col_InfShift = 512*lY*2;

%     row_shift_multi = gather(row_shift_multi);
%     col_shift_multi = gather(col_shift_multi);    
    
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
        
        a = 1; %#ok<*NASGU>
        Y_temp_corrected(:,:,t:(t+size_3_Y_temp_batch-1)) = gather(Mf);
        a = 1;
    end
        
    
    clear CCreal CCreal2 cloc_multi coff_forRefine coff_forRefine_off col_shift_multi conj2
    clear fftY_batch III kernc kernr maxY_batch
    clear median_maxY_batch median_minY_batch Mf
    clear minY_batch Nc2 Nr2 rloc_multi roff_forRefine roff_forRefine_off row_shift_multi
    clear temp1 temp1234_multi temp13 temp13_multi temp2 temp24 temp24_multi temp3 temp4
    clear temp_boolArray temp_boolIndex temp_col_shift_multi temp_constant temp_ifftshift_off temp_ifftshift_off_Xconst
    clear temp_row_shift_multi tempNc tempNc_batch tempNr tempNr_batch Y_temp_batch
    
    a = 1;
    
    %temp_corr = gather(motion_metrics_gpuBatch_jjb_v1(gpuArray(Y_temp_corrected),10));
    temp_corr = motion_metrics_gpuBatch_jjb_v1(Y_temp_corrected,10);
    [B,I] = sort(temp_corr,'descend');
    tempIndex_top = I(1:min(k2,k1));
    Y_temp_corrected_top = Y_temp_corrected(:,:,tempIndex_top);
    template_in = mean(Y_temp_corrected_top,3);
    
    template_nk2(:,:,n) = template_in;
    
    fprintf('n=%d, ',n);
    
    a = 1;
    clear Y_temp
    clear Y_temp_corrected_top
    clear Y_temp_corrected
    clear col_shift_template overMax_index row_shift_template underMin_index
    clear col_leftBound_multi col_leftBound_multi2 col_rightBound_multi col_rightBound_multi2
    clear displacementField displacementField_1d_col displacementField_1d_col_int displacementField_3d_col
    clear displacementField_3d_col_valid displacementField_3d_row temp_bool
    a = 1;
end
fprintf('time (template_nk2) = %.1f secs.\n',toc(t0_fun));

a = 1;

template_nk2; %#ok<*VUNUS>
a = 1;

template_nk2;
template_T;
% D_T_to_nk2 = zeros(options.d1,options.d2,N,'single');
D_T_to_nk2 = cell(1,N);
template_nk2_corrected = zeros(options.d1,options.d2,N,'single');


targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
singleImage_nonrigid_jjb_Name_v = autoGetFunName('singleImage_nonrigid_jjb', [targetPATH '\functions']);
fun_singleImage_nonrigid_jjb= str2func(singleImage_nonrigid_jjb_Name_v);
    
% imregdemons
t1_fun = tic;
% for n=1:N
parfor n=1:N
    template = template_nk2(:,:,n);
    [temp_min_template,temp_max_template] = bounds(template,'all');
    
    a = 1;
    
    [D,~] = fun_singleImage_nonrigid_jjb(template,template_T); %#ok<*PFBNS>
    
    template_corrected = imwarp(template,gather(D),'cubic');%'nearest' would lose sub-pixel displayment!
    overMax_index = template_corrected > temp_max_template;
    underMin_index = template_corrected < temp_min_template;
    template_corrected = template_corrected.*(~overMax_index) + temp_max_template.*overMax_index;
    template_corrected = template_corrected.*(~underMin_index) + temp_min_template.*underMin_index;
    
    template_nk2_corrected(:,:,n) = template_corrected;
    D_T_to_nk2{n} = D;
    a = 1;    
end
fprintf('time (D_T_to_nk2) = %.1f secs.\n',toc(t1_fun));


% if_plot = 1;

if if_plot == 1
    close all
    max_template = max(template_nk2,[],'all');
    for n=1:N
        figure(n);
        set(gcf,'Position',[35+400 35+0 850 850]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set (gca,'position',[0.01,0.01,0.98,0.98])
        template = template_nk2(:,:,n);
        template_n11n = template./max_template;
        
        imshow(template_n11n.^0.4);
        
        colormap(gray);
        axis off equal
    end
end

if if_plot == 1
    max_template_corrected = max(template_nk2_corrected,[],'all');
    for n=1:N
        figure(n+N);
        set(gcf,'Position',[35+400 35+0 850 850]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set (gca,'position',[0.01,0.01,0.98,0.98])
        template = template_nk2_corrected(:,:,n);
        template_n11n = template./max_template_corrected;
        
        imshow(template_n11n.^0.4);
        
        colormap(gray);
        axis off equal
    end
end

a = 1;

%% End