function [cY,mY,ng] = motion_metrics_gpuBatch_jjb_v1(Y,bnd,batch_size,var_name)

% computes several metrics for the quantitative assessment of the registration

% INPUTS
% Y:            registered time series in a loaded or memory mapped array
% bnd:          number of pixels to be exluded to avoid NaN effects 
%                   [x_beg,x_end,y_be,y_end,z_beg,z_end]
% batch_size:   size of batch to be read for memory mapped files
% var_name:     in case of memory mapped files use this variable

% OUTPUTS
% cY:           correlation coefficient of each frame with the mean
% mY:           mean image
% ng:           norm of gradient of mean image

if nargin == 1 || isempty(bnd); bnd = zeros(6,1); end
if nargin < 3|| isempty(batch_size); batch_size = 1000; end


sizY = size(Y);


dimsY = length(sizY);
nd = dimsY - 1;
d = prod(sizY(1:end-1));
T = sizY(end);

if isscalar(bnd); bnd = ones(2*(dimsY-1),1)*bnd; end
if dimsY == 3; sizY(3) = 1; bnd(5:6) = 0; end


Y = single(Y);
nd = ndims(Y)-1;
mY = nanmean(Y,nd+1);
if nd == 2
    Yr = Y(bnd(1)+1:sizY(1)-bnd(2),bnd(3)+1:sizY(2)-bnd(4),:);
    mYr = mY(bnd(1)+1:sizY(1)-bnd(2),bnd(3)+1:sizY(2)-bnd(4));
else
    Yr = Y(bnd(1)+1:sizY(1)-bnd(2),bnd(3)+1:sizY(2)-bnd(4),bnd(5)+1:sizY(3)-bnd(6),:);
    mYr = mY(bnd(1)+1:sizY(1)-bnd(2),bnd(3)+1:sizY(2)-bnd(4),bnd(5)+1:sizY(3)-bnd(6));
end
Yr = reshape(Yr,[],T);
mYr = mYr(:);
if any(any(isnan(Yr))) || any(isnan(mYr));
    cY = corr(Yr,mYr,'rows','p');
else
    cY = corr(Yr,mYr);
end
   

if ismatrix(mY);
    [gx,gy] = gradient(mY); 
    ng = norm(sqrt(gx.^2+gy.^2),'fro');
elseif ndims(Y)-1 == 3;
    [gx,gy,gz] = gradient(mY); 
    ng = sqrt(sum(gx(:).^2+gy(:).^2+gz(:).^2));
end

