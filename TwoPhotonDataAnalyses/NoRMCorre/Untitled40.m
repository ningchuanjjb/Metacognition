nhood = [1 1];
se = strel(nhood);

BW_raw = input_image > 200;
% BW_raw = gather(input_image > 200);

a = 1;

t0 = tic;
for tempi=1:1000
% BW1 = imerode(gpuArray(BW),nhood,'ispacked',512);
BW1 = imerode(BW,se,'ispacked',512);
end
toc(t0);

t1 = tic;
for tempi=1:1000
BW1 = imerode(BW_raw,nhood);
end
toc(t1);
a = 1;

%Conclusion: speed of two methods are similar.