% clear all
close all

% temp_data = data(start_point(i):end_point(i));
% load temp_data

N = 500;
temp_data = [0.5*randn(1,N)+1,1.5*randn(1,N)+3,0.5*randn(1,N)+2];


% temp_data = x';

n=100;
n=2^ceil(log2(n)); % round up n to the next power of 2;

% [bandwidth,density,xmesh,cdf] = kde_jjb_v1(temp_data);
% [bandwidth,density,xmesh,cdf] = kde_jjb_v1(temp_data,n);
% [bandwidth,density,xmesh,cdf] = kde_jjb_v5(temp_data,n);
[bandwidth,density,xmesh,cdf] = kde_jjb_v7(temp_data,n);
% pdf = diff(cdf);
pdf = density;
% pdf = pdf./sum(pdf);
% pdf = pdf./max(pdf);

[pdf2,xmesh2,bandwidth2] = ksdensity(temp_data,'NumPoints',n,'Function','pdf');
% pdf2 = pdf2./sum(pdf2);
% [cdf2,xmesh2,bandwidth2] = ksdensity(temp_data,'NumPoints',n,'Function','cdf','Bandwidth',bandwidth);
% pdf2 = diff(cdf2);
% pdf2 = pdf2./sum(pdf2);
% pdf2 = pdf2./max(pdf2);

[~,temp_mode_index] = max(pdf2);
% [~,temp_median_index] = min(abs(cdf2-median(cdf2)));
% temp_median_index = find(cdf2>0.5,1);


temp_N = 100;
figure(1)
h = histogram(temp_data,temp_N,'Normalization','pdf');
% h = histogram(temp_data,100);
h.BinEdges;
h.Values;
temph = zeros(1,n);
for tempi=1:length(temph)
    temp_bool = xmesh2(tempi) < h.BinEdges;
    temp_index = min(find(temp_bool>0,1),temp_N);
    if isempty(temp_index) == 0
        temph(tempi) = h.Values(temp_index);
    end    
end
% temph2 = temph./sum(temph);
temph2 = temph;
temp_max = max(max(pdf,[],'all'),max(pdf2,[],'all'));

figure(2)
% plot(pdf,'r');
plot(xmesh,pdf,'r');
hold on 
% plot(pdf2,'g');
plot(xmesh2,pdf2,'g');
hold on
plot(xmesh2,temph2,'b');
% hold on
% plot([1,1]*temp_mode_index,[0,temp_max],'k');
% % plot([1,1]*temp_mode_index,[0,1],'k');
% hold on
% plot([1,1]*temp_median_index,[0,temp_max],'c');
% % plot([1,1]*temp_median_index,[0,1],'c');
% hold on



a = 1;

