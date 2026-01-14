function beta_shuffled = shuffleBeta_v2(beta,num_roi,numFrames,num_betaGroup)

a = 1;

beta_mean = zeros(num_roi,numFrames,num_betaGroup);
for temptempj=1:num_roi
    temp_beta_roi = beta(temptempj,:); %#ok<*PFBNS>
    for temptempk=1:num_betaGroup
        temp_range = (1:6)+(temptempk-1)*numFrames;        
        temptemp_beta = temp_beta_roi(temp_range);
        temptemp_beta_mean = mean(temptemp_beta);
        %temptemp_beta_mean = median(temptemp_beta);
        beta_mean(temptempj,:,temptempk) = temptemp_beta_mean;
    end
end
temp_beta_mean_raw = beta_mean; %#ok<*NASGU>
beta_mean = reshape(beta_mean,[num_roi,num_betaGroup*numFrames]);


beta_shuffled = beta_mean;
