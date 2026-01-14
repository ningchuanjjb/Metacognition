close all

temptempi = 1;
temp1 = meta_trialLevel_crossTime_length1(temptempi,:);
temp_interp_factor = 3;

% temptempi = 3;
% temp1 = meta_trialLevel_crossTime_length2(temptempi,:);
% temp_interp_factor = 1.5;


x = 1:length(temp1);
x_interp = (1:temp_interp_factor:length(temp1)*temp_interp_factor)+(temp_interp_factor-1);

temp1_interp = interp1(x_interp,temp1,1:length(temp1)*temp_interp_factor,'makima');

for temptempj=1:ceil((temp_interp_factor-1))
   temp1_interp(temptempj) = temp1_interp(ceil(temp_interp_factor));
end

plot(x_interp,temp1,'color',[0 0 0])
hold on
plot(temp1_interp+0.02,'color',[1 0 0])
hold on