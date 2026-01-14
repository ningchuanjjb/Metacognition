global_tempi_reader = 1;

t_spmd = zeros(1,workerNum);
bin_count_spmd = zeros(1,workerNum);
F_reader(1:workerNum) = parallel.FevalFuture;


bin_count = i_worker_start(1);
t = (bin_count-1)*bin_width_eachWorker + 1;
F_reader(1) = parfeval(backgroundPool,@read_raw_file,1,...
    Y,t,min(t+bin_width_eachWorker-1,T)-t+1,[512,512],2);
[completedId_reader, Ytm] = fetchNext(F_reader(1));
%Ytm = read_raw_file(Y,t,min(t+bin_width_eachWorker-1,T)-t+1,[512,512],2);
size_3_Ytm = size(Ytm, 3);

for spmdindex=2:workerNum
    tempi_worker = i_worker_start(spmdindex):i_worker_end(spmdindex);
    if global_tempi_reader <= length(tempi_worker)
        bin_count_spmd(spmdindex) = tempi_worker(global_tempi_reader);
        t_spmd(spmdindex) = (bin_count_spmd(spmdindex)-1)*bin_width_eachWorker + 1;
        F_reader(spmdindex) = parfeval(backgroundPool,@read_raw_file,1,...
            Y,t_spmd(spmdindex),min(t_spmd(spmdindex)+bin_width_eachWorker-1,T)-t_spmd(spmdindex)+1,[512,512],2);
    end
end
F_reader;
% [completedId_reader, Ytm] = fetchNext(F_reader(2:end));








