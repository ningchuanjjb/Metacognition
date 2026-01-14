

template_n11n = rescale(template,min_template_in_raw,max_template_in_raw);
template_in;

%template_n11n = rescale(template,0,1);% scale to [0,1]
%template_n11n = adapthisteq(template_n11n,'NBins',4096,'Range','original','Distribution','rayleigh','Alpha',0.8);% histogram equalization
%template_n11n = rescale(template_n11n,min_template_in_raw,max_template_in_raw);

%template_in = rescale(template_in,0,1);% scale to [0,1]
%template_in = adapthisteq(gather(template_in),'NBins',4096,'Range','original','Distribution','rayleigh','Alpha',0.8);% histogram equalization
%template_in = rescale(template_in,min_template_in_raw,max_template_in_raw);


initial_template_n11n = rescale(gather(template_n11n),0,1);
initial_template_in = rescale(gather(template_in),0,1);

temp_template_n11n = initial_template_n11n;
temp_template_in = initial_template_in;


% To fit optimal alpha_histeq by fminsearch
range_alpha_histeq = [0 100];
params_range = range_alpha_histeq;

ModelName = 'adapthisteq';
Model = str2func(ModelName);
% func = @(params) Model(params,gather_template_in,options_nonrigid,gather_template_n11n,params_range);
func = @(params) Model(gather(template_in),'NBins',4096,'Range','original','Distribution','rayleigh','Alpha',params);

% Use fminsearch
optimset_options = optimset('fminsearch');
optimset_options.TolFun = 1e-7;%1e-6
optimset_options.TolX = 1e-3;%1e-1
%optimset_options.MaxFunEvals = 40;
optimset_options.MaxIter = 40;
optimset_options.Display = 'off'; % off, iter, final
optimset_options.FunValCheck = 'on';
optimset_options.PlotFcns = [];%@optimplotfval

t0 = tic;

initial_alpha_histeq = 0.4;%0.4
alpha_histeq_raw = ...
    -log(1./((initial_alpha_histeq-range_alpha_histeq(1))./(range_alpha_histeq(2)-range_alpha_histeq(1)))-1);
initialParams = alpha_histeq_raw;
tempN = length(initialParams);
rho = zeros(1,tempN);
temp_modelParams = zeros(1,tempN);
parfor tempi=1:tempN
    temp_optimset_options = optimset_options;
    if tempi==1
        temp_optimset_options.Display = 'iter'; % off, iter, final
    end
    temp_modelParams(tempi) = fminsearch(func, initialParams(tempi), temp_optimset_options);
    [~,~,rho(tempi)] = Model(temp_modelParams(tempi),gather_template_in,options_nonrigid,gather_template_n11n,params_range);
end
a = 1;
[rho2,I]=max(rho);
modelParams = temp_modelParams(I);


alpha_histeq_raw = modelParams(1);
alpha_histeq = (1/(1+exp(-alpha_histeq_raw)))*(range_alpha_histeq(2)-range_alpha_histeq(1))+range_alpha_histeq(1);% constrain in (1, 6)
fprintf('fimsearch time=%.1f seconds, rho=%.5f, alpha_histeq=%.3f(from%.3f)\n',...
    toc(t0),rho2,alpha_histeq,initial_alpha_histeq(I));

gather_template_in = gather_template_in.^(1/alpha_histeq);
gather_template_n11n = gather_template_n11n.^(1/alpha_histeq);
[min_gather_template_in,max_gather_template_in] = bounds(gather_template_in,'all');