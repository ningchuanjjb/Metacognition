clear

% g = gpuDevice(1);
% g = gpuDevice(2);
a = zeros(1024,1024,'single','gpuArray')+1;
tic
for tempi=1:500000
   b = fft(a); 
end
toc