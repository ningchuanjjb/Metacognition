tic
fA = parfeval(backgroundPool,@rand,1,10000);
% fA = parfeval(@rand,1,10000);
% fA = parfevalOnAll(@rand,1,10000);
% fA = parfevalOnAll(backgroundPool,@rand,1,10000);

% B = rand(10000);
fB = parfeval(backgroundPool,@rand,1,10000);

wait(fA)
wait(fB)
C = max(fetchOutputs(fA),B);
%  C = max(fetchOutputs(fA),fetchOutputs(fB));
toc

% tic
% % fA = parfeval(backgroundPool,@rand,1,10000);
% % B = ones(10000);
% % wait(fA)
% % C = max(rand(10000),B);
% C = max(rand(10000),rand(10000));
% toc