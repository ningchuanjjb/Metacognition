function output_batch = dftregistration_min_max_gpuBatch_jjb_v7(buf1ft,buf2ft_multi,usfac,min_shift,max_shift,phase_flag,batch_size,Nr2,Nc2)
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory 
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a 
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007
%
% Rewrote all code not authored by either Manuel Guizar or Jim Fienup
% Manuel Guizar - May 13, 2016
%
% Modified by Eftychios A. Pnevmatikakis to include upper bound on possible
% shifts - November 1, 2016
%
% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).
%
% Inputs
% buf1ft    Fourier transform of reference image, 
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register, 
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)
% max_shift Maximum shift in each dimension (2x1 vector). (default = Inf, no max)
%
% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.
%
%
% Copyright (c) 2016, Manuel Guizar Sicairos, James R. Fienup, University of Rochester
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the University of Rochester nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

if ~exist('phase_flag','var')
    phase_flag = true;
end

if ~exist('usfac','var')
    usfac = 1;
end

if ~exist('max_shift','var')
    min_shift = -Inf(1,2);
end

if ~exist('max_shift','var')
    max_shift = Inf(1,2);
end


[nr,nc]=size(buf1ft);
d1 = nr;
d2 = nc;
pad_size_coeff = 2;
pad_size = single(pad_size_coeff*[d1 d2]);


buf_pad_multi = zeros([pad_size_coeff*[d1 d2] batch_size], 'single', 'gpuArray');
CC_multi = zeros([pad_size_coeff*[d1 d2]  batch_size], 'single', 'gpuArray'); %#ok<*PREALL>
CCabs_multi = zeros([pad_size_coeff*[d1 d2]  batch_size], 'single', 'gpuArray');
row_shift_multi = zeros(1, batch_size, 'single', 'gpuArray');
col_shift_multi = zeros(1, batch_size, 'single', 'gpuArray');
CCmax_multi = zeros(1, batch_size, 'single', 'gpuArray');
rloc_multi = zeros(1, batch_size, 'single', 'gpuArray');
cloc_multi = zeros(1, batch_size, 'single', 'gpuArray');


if isscalar(min_shift); min_shift = min_shift*[1,1]; end
if isscalar(max_shift); max_shift = max_shift*[1,1]; end

conj_buf1ft = conj(buf1ft);
conj_buf2ft_multi = conj(buf2ft_multi);

if usfac == 0
elseif usfac == 1
elseif usfac > 1
    
    buf_prod_multi = buf1ft.*conj_buf2ft_multi;    
%     for tempi=1:batch_size
%        buf_pad_multi(:,:,tempi) = FTpad_jjb_v2(buf_prod_multi(:,:,tempi),pad_size);         
%     end          
    buf_pad_multi = FTpad_gpuBatch_jjb_v4(buf_prod_multi,[pad_size batch_size]);           
    
    if phase_flag
        buf_pad_multi = sign(buf_pad_multi); %#ok<*UNRCH>
    end    
    
    CC_multi= ifft2(buf_pad_multi);    
    
    CCabs_multi = abs(CC_multi);    
    
    [~, III] = max(CCabs_multi, [], [1 2],'linear');
    III = reshape(III,1,batch_size);
    [row_shift_multi,col_shift_multi,~] = ind2sub([pad_size batch_size],III);    
    
    % Now change shifts so that they represent relative shifts and not indices
    Nr2 = ifftshift(gpuArray(-fix(nr):ceil(nr)-1));
    Nc2 = ifftshift(gpuArray(-fix(nc):ceil(nc)-1));       
    
    for tempi=1:batch_size
        if Nr2(row_shift_multi(tempi))/2 > max_shift(1) || ...
                Nc2(col_shift_multi(tempi))/2 > max_shift(2) ||...
                Nr2(row_shift_multi(tempi))/2 < min_shift(1) ||...
                Nc2(col_shift_multi(tempi))/2 < min_shift(2)
            CCabs2 = CCabs_multi(:,:,tempi);
            CCabs2(Nr2/2>max_shift(1),:) = 0;
            CCabs2(:,Nc2/2>max_shift(2)) = 0;
            CCabs2(Nr2/2<min_shift(1),:) = 0;
            CCabs2(:,Nc2/2<min_shift(2)) = 0;
            [row_shift_multi(tempi), col_shift_multi(tempi)] = ...
                find(CCabs_multi(:,:,tempi) == max(CCabs2(:)),1,'first');
        end
    end
    
    for tempi=1:batch_size
        CCmax_multi(tempi) = CC_multi(row_shift_multi(tempi),col_shift_multi(tempi),tempi)*nr*nc;
        row_shift_multi(tempi) = single(Nr2(row_shift_multi(tempi))/2);
        col_shift_multi(tempi) = single(Nc2(col_shift_multi(tempi))/2);
    end         
    
    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid               
        row_shift_multi = round(row_shift_multi*usfac)/usfac;
        col_shift_multi = round(col_shift_multi*usfac)/usfac;        
                
        size_us = ceil(usfac*1.5);        
        dftshift = fix(size_us/2); %% Center of output array at dftshift+1
        
        % Matrix multiply DFT around the current shift estimate
        CC_multi = zeros([size_us size_us batch_size], 'single', 'gpuArray');        
        for tempi=1:batch_size
            %temp_range = (1:d2)+(tempi-1)*d2;
            CC_multi(:,:,tempi) = ...
                conj(...
                    dftups(...
                        buf2ft_multi(:,:,tempi).*conj_buf1ft,...
                        size_us,...
                        size_us,...
                        usfac,...
                        dftshift-row_shift_multi(tempi)*usfac,...
                        dftshift-col_shift_multi(tempi)*usfac...
                    )...
                ); 
        end                
        
        % Locate maximum and map back to original pixel grid         
        CCabs_multi = abs(CC_multi);
        
        [~, III] = max(CCabs_multi, [], [1 2],'linear');
        III = reshape(III,1,batch_size);
        [rloc_multi,cloc_multi,~] = ind2sub([size_us size_us batch_size],III);        
        for tempi=1:batch_size
            CCmax_multi(tempi) = CC_multi(rloc_multi(tempi),cloc_multi(tempi),tempi);
        end        
        rloc_multi = rloc_multi - dftshift - 1;
        cloc_multi = cloc_multi - dftshift - 1;
        row_shift_multi = row_shift_multi + rloc_multi/usfac;
        col_shift_multi = col_shift_multi + cloc_multi/usfac;           
 
    end       
    
    % If its only one row or column the shift along that dimension has no
    % effect. Set to zero.
    if nr == 1
        row_shift_multi(:) = 0;
    end
    if nc == 1
        col_shift_multi(:) = 0;
    end
    
end  


error = 1;
diffphase_multi = angle(CCmax_multi);
output=[error,diffphase_multi,row_shift_multi,col_shift_multi];

output_batch = cell(1, batch_size);
for tempi = 1:batch_size
    output_batch{tempi} = zeros(1, 4, 'single', 'gpuArray');
    output_batch{tempi}(1) = output(1);
    output_batch{tempi}(2) = output(1+tempi);
    output_batch{tempi}(3) = output(1+batch_size+tempi);
    output_batch{tempi}(4) = output(1+batch_size*2+tempi);
end
                
a = 1;

return


