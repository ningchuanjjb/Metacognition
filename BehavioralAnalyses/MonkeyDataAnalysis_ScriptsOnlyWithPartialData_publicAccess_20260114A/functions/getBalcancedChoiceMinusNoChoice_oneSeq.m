function [choiceMinusNoChoice_oneSeq, weight] = getBalcancedChoiceMinusNoChoice_oneSeq(...
    choice_correct_trialNum, choice_trialNum, noChoice_correct_trialNum, noChoice_trialNum)

% choice_correct_trialNum = 5;
% choice_trialNum = 20;
% noChoice_correct_trialNum = 20;
% noChoice_trialNum = 100;

choiceMinusNoChoice_oneSeq = ...
    choice_correct_trialNum/choice_trialNum - noChoice_correct_trialNum/noChoice_trialNum;

weight = min(choice_trialNum/noChoice_trialNum, noChoice_trialNum/choice_trialNum);


end