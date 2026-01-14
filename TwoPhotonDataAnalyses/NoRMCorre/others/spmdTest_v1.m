if false
    delete(gcp); %#ok<*UNRCH>
    gcp;
end

b = zeros(512,512, 'single');
c = fft2(b);
spmd(0, 16)
    %for tempi=1:2500
        
        %a = zeros(1,1, 'single', 'gpuArray');
        %b = gpuArray(zeros(512,512, 'single'));
        %d = ifft2(gpuArray(c));
        %d = gpuArray(c);
    %end
end

% parfor tempi=1:2500
%    b = zeros(512,512,16, 'single', 'gpuArray');             
% end
    
% parfor tempi=1:2500
%     b = zeros(512,512,16, 'single', 'gpuArray');
% end

