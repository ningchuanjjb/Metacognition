function [ imFTout ] = FTpad_gpu_jjb_v3(imFT,outsize)
% imFTout = FTpad(imFT,outsize)
% Pads or crops the Fourier transform to the desired ouput size. Taking 
% care that the zero frequency is put in the correct place for the output
% for subsequent FT or IFT. Can be used for Fourier transform based
% interpolation, i.e. dirichlet kernel interpolation. 
%
%   Inputs
% imFT      - Input complex array with DC in [1,1]
% outsize   - Output size of array [ny nx] 
%
%   Outputs
% imout   - Output complex image with DC in [1,1]
% Manuel Guizar - 2014.06.02

if ~ismatrix(imFT)
    error('Maximum number of array dimensions is 2')
end

% data_type = class(imFT);
% fun_data_type = str2func(data_type);
% outsize = fun_data_type(outsize);

% size_imFT = fun_data_type(size(imFT));

% Nout = outsize;
% Nin = size_imFT;
imFT = fftshift(imFT);
% center = floor(size_imFT/2)+1;

imFTout = zeros(outsize, 'single','gpuArray');

% centerout = floor(outsize/2)+1;


% t0 = tic;
% cenout_cen = centerout - center;
% imFTout(max(cenout_cen(1)+1,1):min(cenout_cen(1)+Nin(1),Nout(1)),max(cenout_cen(2)+1,1):min(cenout_cen(2)+Nin(2),Nout(2))) ...
%     = imFT(max(-cenout_cen(1)+1,1):min(-cenout_cen(1)+Nout(1),Nin(1)),max(-cenout_cen(2)+1,1):min(-cenout_cen(2)+Nout(2),Nin(2)));
%imFTout(257:768,257:768) = imFT(1:512,1:512);
imFTout(257:768,257:768) = imFT;


% fprintf('FTpad time = %.4f seconds\n',toc(t0));


%imFTout = ifftshift(imFTout)*Nout(1)*Nout(2)/(Nin(1)*Nin(2));
imFTout = ifftshift(imFTout)*4;
%a = 1;


% fprintf('FTpad time = %.4f seconds\n',toc(t0));
% a = 1;

% return