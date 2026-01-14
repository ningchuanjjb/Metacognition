function I = shift_reconstruct_batch_jjb_v1(Y,shifts,diffphase,~,Nr,Nc,~,method,add_value)

% applies 3-d sub-pixel shifts to an input image
% INPUTS:
% Y:            input image (double 3d tensor) in space (real) or frequency (complex) domain
% shifts:       shifts 
% diffphase:    phase difference
% us_fac:       upsampling factor for subpixel shifts
% method:       method for treating boundaries
% add_value:    value to add when zero-ing out boundaries

% OUTPUT:
% I:        output image

% Written by Eftychios A. Pnevmatikakis, Simons Foundation, 2016

if isreal(Y)
    buf2ft = fftn(Y);
else
    buf2ft = Y;
end


[nr,nc,np]=size(Y);

row_shift = shifts(1,:);
col_shift = shifts(2,:);
batch_size = length(col_shift);

temp_r = zeros([nr,nc,np], 'single');
temp_c = zeros([nr,nc,np], 'single');
for tempi=1:batch_size
    temp_r(:,:,tempi) = row_shift(tempi)*Nr/nr;
    temp_c(:,:,tempi) = col_shift(tempi)*Nc/nc;    
end

Greg = buf2ft.*exp(1i*2*pi*(-temp_r-temp_c));

for tempi=1:batch_size
    Greg(:,:,tempi) = Greg(:,:,tempi)*exp(1i*diffphase(tempi));   
end

I = real(ifftn(Greg));
for tempi=1:batch_size
    I(:,:,tempi) = remove_boundaries(I(:,:,tempi),shifts(:,tempi),method,add_value);
end

