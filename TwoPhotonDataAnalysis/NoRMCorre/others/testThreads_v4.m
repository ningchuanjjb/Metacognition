%parpool("threads");
clear
% tic
% a = 1;
% parfor i = 1:100
%     %A{i} = rand(100);
%     a = a + i;
% end
% toc

threads0_processes1 = 1;
if_restartPool = 0;

if if_restartPool == 1
    delete(gcp);
    if threads0_processes1 == 1
        parpool('Processes');
    elseif threads0_processes1 == 0
        parpool('threads');
    end
end
pool = gcp;

% X = rand(20000, 20000);
% ticBytes(pool);
% tThreads = timeit(@() fetchOutputs(parfeval(@sum,1,X,'all')));
% f = parfeval(@sum,1,X,'all');
% tThreads = timeit(@() sum(X,'all'));
% tocBytes(pool)
% fprintf('\n%.3f seconds\n',tThreads);

% ticBytes(pool);
% tic
if_parfeval = 1;



tic
N = 160;
if if_parfeval == 1
    F1(1:N) = parallel.FevalFuture;
%     F2(1:N) = parallel.FevalFuture;
    for idx = 1:N
%         F1(idx) = parfeval(@tempASD2,2,2000,2000,idx);
        F1(idx) = parfeval(backgroundPool,@tempASD2,2,2000,2000,idx);

%         F2(idx) = parfeval(@tempASD2,2,2000,2000,idx);
%         F2(idx) = parfeval(backgroundPool,@tempASD2,2,2000,2000,idx);

%         F1(idx) = parfevalOnAll(@tempASD1,2,2000,2000,idx,aa);
        
        %     F1(idx) = parfeval(@rand,1,1000,1000);
        %     F1(idx) = parfevalOnAll(@rand,1,1000,1000);
    end
    fprintf('asd111\n');
    %F1 = parfeval(@tempASD1,2,2000,2000,5,aa,N);
    %F1 = parfevalOnAll(@tempASD1,2,3000,3000,5,N);

    %F1 = parfeval(backgroundPool,@tempASD1,2,3000,3000,5,N);    
    %F1 = parfevalOnAll(backgroundPool,@tempASD1,2,3000,3000,5,N);

    % toc
    % F2(1:N) = parallel.FevalFuture;
    % for idx = 1:N
    %     F2(idx) = parfeval(@tempASD2,2,10000,10000,idx);
    % end
    % toc
    %fprintf('asd111\n');
    %[out1,out12] = fetchOutputs(F1); % 10-by-10 concatenated output
    % toc


    completeNum = 0;
    completedIdx1 = zeros(1, N);
    completedIdx2 = zeros(1, N);
    for idx = 1:N
        [completedIdx1(idx),temp_out1,temp_out2] = fetchNext(F1); %#ok<*ASGLU>
%         [completedIdx2(idx),temp_out1,temp_out2] = fetchNext(F2);

        %completeNum = completeNum + 1;
        %     [temp_out1,temp_out2] = fetchOutputs(F(idx));
        %     temp_out2
        %     toc
    end
    %completedIdx1
else
    for idx = 1:N
       [out1,out12] = tempASD2(2000,2000,idx);
    end
    %[out1,out12] = tempASD1(3000,3000,5,N);
end


toc
% completeNum = length(completedIdx1) + length(completedIdx2);
% fprintf('completeNum = %d\n',completeNum);
% toc
% tocBytes(pool)
% [out1,out2] = tempASD(20,20,10,aa);

function [out1,out2] = tempASD1(in1,in2,in3,n)
for tempi=1:n
    out1 = ones(in1,in2)*in3;
    out2 = in3;
end
end

function [out1,out2] = tempASD2(in1,in2,in3)
    out1 = ones(in1,in2)*in3;
    out2 = in3;
end
