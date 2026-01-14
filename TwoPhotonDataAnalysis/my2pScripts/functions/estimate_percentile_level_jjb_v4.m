function [level,cdf_val] = estimate_percentile_level_jjb_v3(data,window,shift)

% Estimates the level of percentile filtering to be used in DF/F extraction
% using a kernel density estimator. 

if ~exist('window','var'); window = 4000; end
if ~exist('shift','var'); shift = 2000; end

data = data(:);
T = length(data);
start_point = 1:shift:T;
end_point = [window:shift:T,T];

min_ln = min(length(start_point),length(end_point));
start_point = start_point(1:min_ln);
end_point = end_point(1:min_ln);
ln_seg = end_point - start_point + 1;

if ln_seg(end) < 2*window/3
    start_point(end) = [];
    end_point(end-1) = [];
end

level = zeros(length(start_point),1);
cdf_val = zeros(length(start_point),1);

a = 1;
for i = 1:length(start_point)
    temp_data = data(start_point(i):end_point(i));
    %if sum(temp_data) == 0
    %   a = 1;
    %end
    if length(unique(temp_data)) == 1
       level(i) = nan;
       cdf_val(i) = nan;
       a = 1;
       continue
    end
    %[~,density,xmesh,cdf] = kde_jjb_v1(data(start_point(i):end_point(i)));
    %[~,density,xmesh,cdf] = kde_jjb_v5(data(start_point(i):end_point(i)),100); %#ok<*ASGLU>
    [~,density,xmesh,cdf] = kde_jjb_v7(data(start_point(i):end_point(i)),1000); %#ok<*ASGLU>
    [~,ind] = max(density);
    
    %ind = find(cdf>percentile,1);
    
    
    %[cdf,xmesh,~] = ksdensity(data(start_point(i):end_point(i)),'NumPoints',16384,'Function','cdf');
    %temp_data = data(start_point(i):end_point(i));    
    %[cdf,xmesh,~] = ksdensity(temp_data,'NumPoints',length(temp_data),'Function','cdf');
    %ind = find(cdf>0.5,1);

    
    
    level(i) = xmesh(ind);
    cdf_val(i) = cdf(ind)*100;
    a = 1;
end
a = 1;

level(isnan(level)) = median(level,'omitnan');
cdf_val(isnan(cdf_val)) = median(cdf_val,'omitnan');

