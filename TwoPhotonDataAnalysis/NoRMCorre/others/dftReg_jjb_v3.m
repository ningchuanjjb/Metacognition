function [mem_buffer_cpu,bin_count,t,size_3_Ytm,readerIndex] = dftReg_jjb_v3(...
    Ytm,bin_count,t,size_3_Ytm,readerIndex,Nr,Nc,phase_flag,fftTemp_mat,us_fac)

%input
% Ytm;
% size_3_Ytm;
% max_shift;
% max_dev_g;
% Nr;
% Nc;
% phase_flag;
% fftTemp_mat;
% us_fac;

%output
% mem_buffer_cpu;

%% Preprocess data in gpu
Ytm = gpuArray(Ytm);
fftTemp_mat = gpuArray(fftTemp_mat);
lY = size_3_Ytm;
[minY_batch,maxY_batch] = bounds(Ytm,[1 2]);
fftY_batch = fft2(Ytm);


%% dftregistration_min_max
Nr2 = ifftshift(gpuArray(-512:511));
Nc2 = ifftshift(gpuArray(-512:511));
tempNr = gpuArray(Nr{1}/512);
tempNc = gpuArray(Nc{1}/512);

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

temp_row_shift_multi = reshape(row_shift_multi,1,1,lY);
temp_row_shift_multi = repmat(temp_row_shift_multi,512,512,1);
temp_col_shift_multi = reshape(col_shift_multi,1,1,lY);
temp_col_shift_multi = repmat(temp_col_shift_multi,512,512,1);

% Shift reconstruct
Mf = real(ifft2(...
    fftY_batch.*...
    exp(...
    1i*2*pi*(-temp_row_shift_multi.*tempNr_batch-temp_col_shift_multi.*tempNc_batch)...
    )));


temp_max = maxY_batch(1);
overMax_index = logical((sign((Mf-temp_max)-0.5)+1)/2);
Mf(overMax_index) = temp_max;

mem_buffer_cpu = gather(pagetranspose(Mf));