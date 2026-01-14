function score = fun_seqSimilarity_v2(boolIndex_location_A,boolIndex_location_B,spatialInterferenceWeight) %#ok<*INUSD>

% spatialInterferenceWeight = 0.5;

% boolIndex_location_A = [1 1 0 1 0 1] == 1;
% % boolIndex_location_B = [1 0 1 1 1 1] == 1;
% boolIndex_location_B = [1 1 0 1 0 1] == 1;


if_16_near = 1;

numFrames = length(boolIndex_location_A);

temp_score1 = sum(boolIndex_location_A == boolIndex_location_B);

temptempBoolIndexA = boolIndex_location_A == 1 & boolIndex_location_B == 0;
temptempBoolIndexB = boolIndex_location_A == 0 & boolIndex_location_B == 1;

temptempIndexA = find(temptempBoolIndexA==true);
temptempIndexB = find(temptempBoolIndexB==true);

temp_spatialInterferenceNum = 0;
if length(temptempIndexA) > 0 %#ok<*ISMT>
    for tempi=1:length(temptempIndexA)
        if temptempIndexA(tempi)-1 >= 1
            if temptempBoolIndexB(temptempIndexA(tempi)-1) == true
                temp_spatialInterferenceNum = temp_spatialInterferenceNum + 1;
            end
        end
        if if_16_near == 1
            if temptempIndexA(tempi)-1 == 0
                if temptempBoolIndexB(numFrames) == true
                    temp_spatialInterferenceNum = temp_spatialInterferenceNum + 1;
                end
            end
        end
        if temptempIndexA(tempi)+1 <= numFrames
            if temptempBoolIndexB(temptempIndexA(tempi)+1) == true
                temp_spatialInterferenceNum = temp_spatialInterferenceNum + 1;
            end
        end
        if if_16_near == 1
            if temptempIndexA(tempi)+1 == (numFrames+1)
                if temptempBoolIndexB(1) == true
                    temp_spatialInterferenceNum = temp_spatialInterferenceNum + 1;
                end
            end
        end
    end
end
temp_spatialInterferenceNum = min([temp_spatialInterferenceNum,length(temptempIndexA),length(temptempIndexB)]);

temp_score2 = spatialInterferenceWeight * temp_spatialInterferenceNum;

score = temp_score1 + temp_score2;
