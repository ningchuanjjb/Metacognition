

% template_n11n = rescale(template,min_template_in_raw,max_template_in_raw);
% template_in;

%template_n11n = rescale(template,0,1);% scale to [0,1]
%template_n11n = adapthisteq(template_n11n,'NBins',4096,'Range','original','Distribution','rayleigh','Alpha',0.8);% histogram equalization
%template_n11n = rescale(template_n11n,min_template_in_raw,max_template_in_raw);

%template_in = rescale(template_in,0,1);% scale to [0,1]
%template_in = adapthisteq(gather(template_in),'NBins',4096,'Range','original','Distribution','rayleigh','Alpha',0.8);% histogram equalization
%template_in = rescale(template_in,min_template_in_raw,max_template_in_raw);


% To fit nBins_histeq
gather_template_in = gather(template_in);
gather_template_n11n = gather(template_n11n);       
initial_template_n11n = rescale(gather_template_n11n,0,1);
initial_template_in = rescale(gather_template_in,0,1);

temp_start = 8;
temp_end = 13;
temp_step = 0.05;
tempN = floor((temp_end-temp_start)/temp_step);
temp_nBins_histeq_power = temp_start:temp_step:temp_end;
rho = zeros(1,tempN);

temp_options_nonrigid = options_nonrigid; 
t0_temp = tic;
parfor tempi=1:tempN
    temp_template_n11n = adapthisteq(initial_template_n11n,'NBins',floor(2^(temp_nBins_histeq_power(tempi))),'Range','original','Distribution','rayleigh','Alpha',0.4);
    temp_template_in = adapthisteq(initial_template_in,'NBins',floor(2^(temp_nBins_histeq_power(tempi))),'Range','original','Distribution','rayleigh','Alpha',0.4);       
    temp_output = normcorre_jjb_v2(temp_template_in,temp_options_nonrigid,temp_template_n11n);
    rho(tempi) = corr2(temp_output,temp_template_n11n);
end
[~,I]=max(rho);
optimal_nBins_histeq_power = temp_nBins_histeq_power(I);
optimal_nBins_histeq = floor(2^optimal_nBins_histeq_power);
fprintf(' optimal_nBins_histeq=%d, time=%.1f secs.',optimal_nBins_histeq,toc(t0_temp));

temp_template_n11n = adapthisteq(initial_template_n11n,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',0.4);
temp_template_in = adapthisteq(initial_template_in,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',0.4);
gather_template_n11n = rescale(temp_template_n11n,min_template_in_raw,max_template_in_raw); %#ok<*NASGU>
gather_template_in = rescale(temp_template_in,min_template_in_raw,max_template_in_raw);


% To fit alpha_histeq
gather_template_in = gather(template_in);
gather_template_n11n = gather(template_n11n);
initial_template_n11n = rescale(gather_template_n11n,0,1);
initial_template_in = rescale(gather_template_in,0,1);

temp_start = 0.1;%0.1
temp_end = 1;%10-->5-->3-->1
temp_step = 0.01;%0.1-->0.025-->0.01
tempN = floor((temp_end-temp_start)/temp_step);
temp_alpha_histeq = temp_start:temp_step:temp_end;
rho = zeros(1,tempN);

temp_options_nonrigid = options_nonrigid; 
t0_temp = tic;
parfor tempi=1:tempN
    temp_template_n11n = adapthisteq(initial_template_n11n,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',temp_alpha_histeq(tempi));
    temp_template_in = adapthisteq(initial_template_in,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',temp_alpha_histeq(tempi));       
    temp_output = normcorre_jjb_v2(temp_template_in,temp_options_nonrigid,temp_template_n11n);
    rho(tempi) = corr2(temp_output,temp_template_n11n);
end
[~,I]=max(rho);
optimal_alpha_histeq = temp_alpha_histeq(I);
fprintf(' optimal_alpha_histeq=%.3f, time=%.1f secs.',optimal_alpha_histeq,toc(t0_temp));

temp_template_n11n = adapthisteq(initial_template_n11n,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);
temp_template_in = adapthisteq(initial_template_in,'NBins',optimal_nBins_histeq,'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq);
gather_template_n11n = rescale(temp_template_n11n,min_template_in_raw,max_template_in_raw);
gather_template_in = rescale(temp_template_in,min_template_in_raw,max_template_in_raw);

