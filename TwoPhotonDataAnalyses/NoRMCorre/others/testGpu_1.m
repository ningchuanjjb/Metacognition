profile on

X = rand(20000);

%G1 = gpuArray(single(X));
G2 = single(gpuArray(X));
G1 = gpuArray(single(X));
G2 = single(gpuArray(X));
G1 = gpuArray(single(X));
G2 = single(gpuArray(X));
G1 = gpuArray(single(X));
G2 = single(gpuArray(X));

profile viewer