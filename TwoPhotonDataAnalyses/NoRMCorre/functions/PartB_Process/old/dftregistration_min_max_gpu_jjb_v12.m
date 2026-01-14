function [output, Greg] = dftregistration_min_max_gpu_jjb_v12(buf1ft,buf2ft,usfac,min_shift,max_shift,phase_flag,if_FTpad,Greg_flag,Nr2,Nc2)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
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

% if ~exist('phase_flag','var')
%     phase_flag = true;
% end
% 
% if ~exist('usfac','var')
%     usfac = 1;
% end
% 
% if ~exist('max_shift','var')
%     min_shift = -Inf(1,2);
% end
% 
% if ~exist('max_shift','var')
%     max_shift = Inf(1,2);
% end


% [512,512]=size(buf2ft);

buf_prod = buf1ft.*conj(buf2ft);% To compute image cross correlation in complex domain

% toc(t0)
% a = 1;
% t1=tic;

% if usfac == 0
%     % Simple computation of error and phase difference without registration
% elseif usfac == 1
%     % Single pixel registration
% elseif usfac > 1
    %t01=tic;

    
    % Start with usfac == 2
    if if_FTpad == 1
        FTpad_coeff = 2;%2
        %buf_pad = FTpad_gpu_jjb_v3(buf_prod,FTpad_coeff*[512,512]); 
        buf_pad = FTpad_gpu_jjb_v4(buf_prod); 
    else
        FTpad_coeff = 1;
        buf_pad = buf_prod;
    end
    
    
    %toc(t01)
    %a = 1;
    %t02=tic;
    if phase_flag
        %buf_pad = buf_pad./(abs(buf_pad)+1e-10);
        buf_pad = sign(buf_pad);
    end
    
    %toc(t02)
    %a = 1;
    %t03=tic;
    
    CC = ifft2(buf_pad);

    
    %toc(t03)
    %a = 1;
    
    CCabs = abs(CC);
    
    %toc(t01)
    %a = 1;     
    %t04=tic;
    
    %[row_shift, col_shift] = find(CCabs == max(CCabs(:)),1,'first');   
    [MMM, III] = max(CCabs, [], 'all','linear'); %#ok<*ASGLU>
    %[row_shift,col_shift] = ind2sub(size(CCabs),III);
    [row_shift,col_shift] = ind2sub([1024 1024],III);

    
    %toc(t04)
    %a = 1;    
    %t05=tic;
    
    % Now change shifts so that they represent relative shifts and not indices
    %Nr2 = ifftshift(gpuArray(-fix(512):ceil(512)-1));
    %Nc2 = ifftshift(gpuArray(-fix(512):ceil(512)-1));
    
    %toc(t05)
    %a = 1;
    %t06=tic;
    
%     if Nr2(row_shift)/FTpad_coeff > max_shift(1) || Nc2(col_shift)/FTpad_coeff > max_shift(2) || Nr2(row_shift)/FTpad_coeff < min_shift(1) || Nc2(col_shift)/FTpad_coeff < min_shift(2)
%         CCabs2 = CCabs;
%         CCabs2(Nr2/FTpad_coeff>max_shift(1),:) = 0;
%         CCabs2(:,Nc2/FTpad_coeff>max_shift(2)) = 0;
%         CCabs2(Nr2/FTpad_coeff<min_shift(1),:) = 0;
%         CCabs2(:,Nc2/FTpad_coeff<min_shift(2)) = 0;
%         [row_shift, col_shift] = find(CCabs == max(CCabs2(:)),1,'first');
%     end     
    %CCmax = CC(row_shift,col_shift)*512*512;
    CCmax = CC(row_shift,col_shift)*262144;
    row_shift = single(Nr2(row_shift)/FTpad_coeff);
    col_shift = single(Nc2(col_shift)/FTpad_coeff);
    
    
    %toc(t06)
    %a = 1;
    %t07=tic;    
    
    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac; 
        col_shift = round(col_shift*usfac)/usfac;   
        size_us = ceil(usfac*1.5);
        dftshift = fix(size_us/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups_jjb_v1(buf2ft.*conj(buf1ft),size_us,size_us,usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac));
        % Locate maximum and map back to original pixel grid 
        CCabs = abs(CC);
        %[rloc, cloc] = find(CCabs == max(CCabs(:)),1,'first');        
        [MMM, III] = max(CCabs, [], 'all','linear');
        [rloc,cloc] = ind2sub(size(CCabs),III);
        
        
        CCmax = CC(rloc,cloc);
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;   

    end
    %toc(t07)
    %a = 1;    
                
%     toc(t03)
%     a = 1;
% end  

% toc(t1)
% a = 1;
% t2=tic;

error = 1;
diffphase = angle(CCmax);

output=[error,diffphase,row_shift,col_shift];

% Compute registered version of buf2ft
if Greg_flag == 1
    if (nargout > 1)&&(usfac > 0)
    %if usfac > 0
        %t002=tic;
        % Little slow
        Nr = ifftshift(-fix(512/2):ceil(512/2)-1);
        Nc = ifftshift(-fix(512/2):ceil(512/2)-1);
        
        [Nc,Nr] = meshgrid(Nc,Nr);
        Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/512-col_shift*Nc/512));
        Greg = Greg*exp(1i*diffphase);
        %toc(t002)
        %a = 1;
    elseif (nargout > 1)&&(usfac == 0)
    %elseif usfac == 0
        Greg = buf2ft*exp(1i*diffphase);
    end
else
    Greg = 0;
end

% return


