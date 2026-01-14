%% Begin
%Enable parallel controller ioport so that giveWater can work certainly.
% !setpci -d 9710:9900 command=101
%giveWater(160000);
%giveWater(80000);
%giveWater(30000);
%giveWater(20000);
%giveWater(15000);
%giveWater(10000);
%giveWater(5000);
%giveWater(1000);

%reCenterGazeBias_x;  
% Clear the workspace and the screen
sca;
close all;
clear all; %#ok<CLALL>

searchName_reCenterGazeBias = 'reCenterGazeBias_Ding.mat';

% preFixationTime = 1;

Screen('Preference', 'VisualDebuglevel', 3);% disable the startup screen, replace it by a black display until calibration is finished
% Bring the Command Window to the front if it is already open
if ~IsOctave; commandwindow; end
%% Basic adjustable parameters begin
%%%%%%%%%%%%%%%%%% basic adjustable parameters begin %%%%%%%%%%%%%%%%%%
[dynamicFixThreshold, basic_para.dynamicFixThreshold] = deal(0.45);
[dynamicFix, basic_para.dynamicFix] = deal(2000);



[taskMode_noChoice_Green0_Red1, basic_para.taskMode_noChoice_Green0_Red1] = deal(0);

[pointShowTime, basic_para.pointShowTime] = deal(600);%437.5-->432.5-->430-->433.5-->500-->470-->500-->400-->380
%450-->470-->550-->600-->725-->750-->600-->600
[pointShowPWM, basic_para.pointShowPWM] = deal(0.666);%0.8-->0.5-->0.43-->0.45-->0.5-->1-->0.5-->0.69-->0.733-->0.666-->0.6923
% here PWM is off/all

[trial_num, basic_para.trial_num] = deal(3000);
%mode flags begin
[if_221False_329True, basic_para.if_221False_329True] = deal(1);
[if_eyeResponse, basic_para.if_eyeResponse] = deal(1);
[if_mouse_mannual_eyeResponse, basic_para.if_mouse_mannual_eyeResponse] = deal(1);%0

[if_choiceToR, basic_para.if_choiceToR] = deal(0);%0


[if_length1_alone, basic_para.if_length1_alone] = deal(0);
[length1_runNum, basic_para.length1_runNum] = deal(5);%5

[if_justFix0_touch1, basic_para.if_justFix0_touch1] = deal(1);%1
[if_fixTolerence, basic_para.if_fixTolerence] = deal(1);
[fixTolerenceTime, basic_para.fixTolerenceTime] = deal(50);%300
[fixTolerenceTime_delay2, basic_para.fixTolerenceTime_delay2] = deal(50);%300-->50-->150-->80-->50
if if_justFix0_touch1 == 0
    [ifSmallWater, basic_para.ifSmallWater] = deal(1);
elseif if_justFix0_touch1 == 1
    [ifSmallWater, basic_para.ifSmallWater] = deal(0);%0    
end
[ifDelay1OffWater, basic_para.ifDelay1OffWater] = deal(0);
[ifLongWater, basic_para.ifLongWater] = deal(0);

[ifDrawDots_eye, basic_para.ifDrawDots_eye] = deal(0);
[ifDrawDots_touch, basic_para.ifDrawDots_touch] = deal(0);
[if_freeview0_fixation1, basic_para.if_freeview0_fixation1] = deal(1);%1
[if_fix_noHold0_hold1, basic_para.if_fix_noHold0_hold1] = deal(0);%1


[joyMode_JJBJoy0_JoyMex1, basic_para.joyMode_JJBJoy0_JoyMex1] = deal(1);

[ifLinkTriangleSquare, basic_para.ifLinkTriangleSquare] = deal(0);
[ifLinkTriangleSquare_inOneTrial, basic_para.ifLinkTriangleSquare_inOneTrial] = deal(1);
[ifCircle, basic_para.ifCircle] = deal(1);
[ifAllowRedMiddle, basic_para.ifAllowRedMiddle] = deal(0);
[ifGreenRedAutoSwitch, basic_para.ifGreenRedAutoSwitch] = deal(0);

[ifSelection, basic_para.ifSelection] = deal(1);
[ifSelectRandom, basic_para.ifSelectRandom] = deal(1);
[ifStep0, basic_para.ifStep0] = deal(0);

[ifSeqLength14_sparse, basic_para.ifSeqLength14_sparse] = deal(0);
[ifSeqLength1_sparse, basic_para.ifSeqLength1_sparse] = deal(0);
[ifSeq_length_3to1, basic_para.ifSeq_length_3to1] = deal(0);%please let seq_length_rangeHead == 3

[ifTurnOff, basic_para.ifTurnOff] = deal(1);
[ifBeginEnd, basic_para.ifBeginEnd] = deal(0);

[ifSpecializedTraining_R, basic_para.ifSpecializedTraining_R] = deal(0);
[ifSpecializedTraining_G, basic_para.ifSpecializedTraining_G] = deal(0);
[ifFormalTask, basic_para.ifFormalTask] = deal(1);

% [if_seqList, basic_para.if_seqList] = deal(1); % pls set 1 with formal test!
if ifFormalTask == 1
    [if_seqList, basic_para.if_seqList] = deal(1);%1
elseif ifFormalTask == 0
    [if_seqList, basic_para.if_seqList] = deal(0);
end
[if_seqList, basic_para.if_seqList] = deal(1);

[ifFreeTouch_origin, basic_para.ifFreeTouch_origin] = deal(1);
[ifNoBackTouch, basic_para.ifNoBackTouch] = deal(1);%1
[ifCorrectLock, basic_para.ifCorrectLock] = deal(0);
[ifCorrectLock_last_toErrorStop, basic_para.ifCorrectLock_last_toErrorStop] = deal(0);
[if_R_mustCorrect, basic_para.if_R_mustCorrect] = deal(0);
[if_R_mustErrorStop, basic_para.if_R_mustErrorStop] = deal(0);%0
[if_validHoldUntilCorrect, basic_para.if_validHoldUntilCorrect] = deal(0);


[if_length4_freetouchCCL0, basic_para.if_length4_freetouchCCL0] = deal(0);
[if_length3_lessCCL, basic_para.if_length3_lessCCL] = deal(0);
[if_length1_moreCCL, basic_para.if_length1_moreCCL] = deal(0);

[if_breakWithNewSeq, basic_para.if_breakWithNewSeq] = deal(1);%0
[if_cancelBreakForCCL, basic_para.if_cancelBreakForCCL] = deal(1);%0

%mode flags end


% Some eyelink paremeters
[dummymode, basic_para.dummymode] = deal(0);
[meanGazeNum, basic_para.meanGazeNum] = deal(4);
[fixWinSize_DisplayPC, basic_para.fixWinSize_DisplayPC] = deal(63);%100-->70-->140-->60-->45-->55
[fixWinSize_ELHost, basic_para.fixWinSize_ELHost] = deal(63);%100-->70-->140-->60-->45-->55

[fixWinSize_delay2, basic_para.fixWinSize_delay2] = deal(193);%%400-->187-->193


[smallWater, basic_para.smallWater] = deal(20);%9-->6-->5-->8
[smallWaterInterval, basic_para.smallWaterInterval] = deal(1500);%500-->200-->500-->750-->1000
[totalWater, basic_para.totalWater] = deal(7000/smallWaterInterval*smallWater);

% [initGiveWaterTime, basic_para.initGiveWaterTime] = deal(170);%160-->180-->190-->180-->200-->120-->180-->165
if if_221False_329True == 0
    [initGiveWaterTime, basic_para.initGiveWaterTime] = deal(245);%
elseif if_221False_329True == 1
    [initGiveWaterTime, basic_para.initGiveWaterTime] = deal(240);
    %240-->260-->300-->270-->285-->277
    %277-->346-->325385-->366-->379
    %379-->290-->230-->165-->174
end
ifDelay1OffWater;
[delay1Water, basic_para.delay1Water] = deal(40);%100
[saccadeWater, basic_para.saccadeWater] = deal(20);%20-->10-->5-->20

% [submitHoldingTime, basic_para.submitHoldingTime] = deal(15);%100-->75-->50-->35-->25-->20-->15
[submitHoldingTime, basic_para.submitHoldingTime] = deal(400);%600-->400
[eyeResponseTime, basic_para.eyeResponseTime] = deal(345);%300-->450-->420-->400-->360-->330-->345
[fixupdateTimeBound, basic_para.fixupdateTimeBound] = deal(270); %150-->175-->205-->270
[fixupdateTimeBound_first, basic_para.fixupdateTimeBound_first] = deal(270);%270-->220
[fixupdateTimeFoldNum, basic_para.fixupdateTimeFoldNum] = deal(10);%3-->10
[fixupdateTimeBound_lock, basic_para.fixupdateTimeBound_lock] = deal(180);%180-->210

fixupdateTimeBound_fold = round(linspace(fixupdateTimeBound/(fixupdateTimeFoldNum+1),fixupdateTimeBound,fixupdateTimeFoldNum+1));

[InternalMemoryWaterRatio, basic_para.InternalMemoryWaterRatio] = deal(1.3);%1.5-->1.4-->1.32-->1.22-->1.3
[waterDropNum, basic_para.waterDropNum] = deal(1);%2
[waterDropInterval, basic_para.waterDropInterval] = deal(600);%ms, 500-->600

[endBetRatio_correctFast, basic_para.endBetRatio_correctFast] = deal(1);
[endBetRatio_correctSlow, basic_para.endBetRatio_correctSlow] = deal(1);
[endBetRatio_errorFast, basic_para.endBetRatio_errorFast] = deal(0);
[endBetRatio_errorSlow, basic_para.endBetRatio_errorSlow] = deal(0);

[endBet_fastSlowBoundary, basic_para.endBet_fastSlowBoundary] = deal(700);



[waterStepPer100, basic_para.waterStepPer100] = deal(0);%10-->8-->10-->12-->15-->12-->10-->8-->0
[stimulusExtendTime, basic_para.stimulusExtendTime] = deal(0);

[seq_length_rangeHead, basic_para.seq_length_rangeHead] = deal(1);
[seq_length_rangeTail, basic_para.seq_length_rangeTail] = deal(4);%3-->4
[numFrames_rangeHead, basic_para.numFrames_rangeHead] = deal(6);
[numFrames_rangeTail, basic_para.numFrames_rangeTail] = deal(6);

% if seq_length_rangeHead ~= 1
%     if_seqList = 0;
% end

[probFreeTouch, basic_para.probFreeTouch] = deal(1);%0-->0.4

[choiceCondition_boardLine1, basic_para.choiceCondition_boardLine1] = deal(0.25);%0.35
[choiceCondition_boardLine2, basic_para.choiceCondition_boardLine2] = deal(0.50);%0.35
[freeTouch_RT_limit, basic_para.freeTouch_RT_limit] = deal(25000);%250000
[selecting_RT_limit, basic_para.selecting_RT_limit] = deal(7000);


[touchBufferSize, basic_para.touchBufferSize] = deal(4);
[touchBuffer_enable, basic_para.touchBuffer_enable] = deal(0);


%[InitAngle_arc, basic_para.InitAngle_arc] = deal(-45);%15-->-15-->-70-->-45
[InitAngle_arc, basic_para.InitAngle_arc] = deal(-59);%after correction:-22.5-->0-->-22.5-->-30-->-50-->-67.5-->-80-->-50
%-63-->-57.5-->-40(numFrames=6)-->-50(numFrames=7)-->-59
[cumulativeCorrectLimit_seqLengthSwitch, basic_para.cumulativeCorrectLimit_seqLengthSwitch] = deal(1);
[cumulativeCorrectLimit_greenRedSwitch, basic_para.cumulativeCorrectLimit_greenRedSwitch] = deal(1);
[consecutiveCorrectLimit, basic_para.consecutiveCorrectLimit] = deal(0);%0

[radius, basic_para.radius] = deal(283);%215-->250-->275-->300-->315-->300-->283
[selectionYpos_coefficient, basic_para.selectionYpos_coefficient] = deal(0.405);%0.435-->0.385-->0.405-->0.385-->0.445
%0.425-->0.405
if if_221False_329True == 0
    [screenXpixels_bias, basic_para.screenXpixels_bias] = deal(0);%60+150-75-->0-->-50-->-15-->0
    [screenYpixels_bias, basic_para.screenYpixels_bias] = deal(-435);%0-->50-->0-->-50-->-140-->-125-->-75-->-150-->-190-->-170-->-160-->-240-->-280-->-255
    %-275-->-285-->-295-->-310-->-325-->-310-->-330-->-380-->-410-->-450-->-510
    %590, %move screen up 160 pix,
    %-590-->-430(horizon)-->390-->350-->325-->300-->330-->320-->360-->420-->-435
elseif if_221False_329True == 1
    [screenXpixels_bias, basic_para.screenXpixels_bias] = deal(0);%35-->0
    [screenYpixels_bias, basic_para.screenYpixels_bias] = deal(-230);%-100-->-70-->-200-->-230
end
[enableRadiantDetect, basic_para.enableRadiantDetect] = deal(0);
[radiantDetectR1, basic_para.radiantDetectR1] = deal(radius-85);%150-->100-->85
[radiantDetectR2, basic_para.radiantDetectR2] = deal(radius+85);%150-->100-->85
    

[ifExplore, basic_para.ifExplore] = deal(0);
[propExplore, basic_para.propExplore] = deal(0);

%-----------


[shortShowTime, basic_para.shortShowTime] = deal(16.6667);%refresh period of screen /ms
%[pointGocueTime, basic_para.pointGocueTime] = deal(500);
% [ItiTime, basic_para.ItiTime] = deal(900);%750-->400

rng('default');%rng('shuffle');
[randArray1, basic_para.randArray1] = deal(rand(1, trial_num)); %for delay1
[randArray2, basic_para.randArray2] = deal(rand(1, trial_num)); %for noChoice and choice
[randArray3, basic_para.randArray3] = deal(rand(1, trial_num)); %for delay2
[randArray4, basic_para.randArray4] = deal(rand(1, trial_num)); %for free touch
[randArray5, basic_para.randArray5] = deal(rand(1, trial_num)); %for if_breakWithNewSeq
[randArray6, basic_para.randArray6] = deal(rand(1, trial_num)); %for ITI
[randArray7, basic_para.randArray7] = deal(rand(1, trial_num)); %
temp_rand = randi(5,1,2000);

% [ifFreeTouch, basic_para.ifFreeTouch] = deal(ifFreeTouch_origin * ones(1, trial_num) .* (randArray4<probFreeTouch));

[ItiTime, basic_para.ItiTime] = deal(0*randArray6 + 1000*ones(1, trial_num));%1100-->800~1200-->800~1000-->2800~3000-->1100~1300-->1000
[preFixationTime, basic_para.preFixationTime] = deal(700*ones(1, trial_num));%300-->700
% [preFixationTime, basic_para.preFixationTime] = deal(800*randArray1 + 100*ones(1, trial_num));
%[fixationTime, basic_para.fixationTime] = deal(400*randArray1 + 450*ones(1, trial_num));%500-->400~600-->600~1000 /2-->400~800-->800~1200-->950~1450
if if_justFix0_touch1 == 0
    %[fixationTime, basic_para.fixationTime] = deal(400*randArray1 + 6800*ones(1, trial_num));
    [fixationTime, basic_para.fixationTime] = deal(400*randArray1 + 20000*ones(1, trial_num));
elseif if_justFix0_touch1 == 1
    [fixationTime, basic_para.fixationTime] = deal(200*randArray1 + 1500*ones(1, trial_num));%1100~1500-->1300~1700-->1500~1900-->1300~1700-->1500~1700
end
[twoFixationTime, basic_para.twoFixationTime] = deal(200*randArray3 + 1500*ones(1, trial_num));
%1100~1500-->700~800-->1500~1900-->1500~1700



[numSquares, basic_para.numSquares] = deal(6);


[cumulativeErrorLimit, basic_para.cumulativeErrorLimit] = deal(4 * ones(1, trial_num));%6-->7-->5-->3-->2-->7
[encourageCumuErrorLimit, basic_para.encourageCumuErrorLimit] = deal(4);%5-->10-->5
[encourageStimET, basic_para.encourageStimET] = deal(20000);


% Set the line width for our fixation cross
[lineWidthPix, basic_para.lineWidthPix] = deal(10);

% Here we set the size of the arms of our fixation cross
[fixCrossDimPix, basic_para.fixCrossDimPix] = deal(40);
[crossLength, basic_para.crossLength] = deal(24);%100-->45-->25-->15
[crossWidth, basic_para.crossWidth] = deal(12);%20-->12-->7-->5
[fixSquareSize, basic_para.fixSquareSize] = deal(100);
[fixTriRadius, basic_para.fixTriRadius] = deal(66.7);

[frameWidth, basic_para.frameWidth] = deal(5);%3-->5


[showBaseRect, basic_para.showBaseRect] = deal([0 0 100 100]);%100-->75-->100-->135-->125-->115-->125-->120-->110
%110-->100-->85-->100-->85-->75-->85
[detectBaseRect, basic_para.detectBaseRect] = deal([0 0 180 180]);%145-->135-->125-->145-->155-->165-->190
[firstDetectRectModifed, basic_para.firstDetectRectModifed] = deal(0);%45-->80-->45-->20-->0
[showBaseTriRadius, basic_para.showBaseTriRadius] = deal(90);
[selectbaseRect, basic_para.selectbaseRectSize] = deal([0 0 110 110]);%125-->110-->95-->110
[selectbaseRect_detection, basic_para.selectbaseRect_detection] = deal([0 0 150 150]);%220-->250-->150
[selectBaseTriRadius, basic_para.selectBaseTriRadius] = deal(80);%90-->80-->67.5-->80
[submitBaseRect, basic_para.submitBaseRect] = deal([0 0 70 70]);%100-->90-->70
[submitBaseRect_detection, basic_para.submitBaseRect_detection] = deal([0 0 130 130]);%150-->95-->120-->130
[holdBaseRect, basic_para.holdBaseRect] = deal([0 0 250 120]);
[holdBaseRect_detection, basic_para.holdBaseRect_detection] = deal([0 0 250 120]);



[bgColor, basic_para.bgColor] = deal([1 1 1]*0.15);%0.25-->0.15
% [normalColor, basic_para.normalColor] = deal([0 0.5 0 0.5]);%[0 1 0 0.5]
[normalColor, basic_para.normalColor] = deal([0 170/255 85/255 0.5]);%[0 1 0 0.5]
[specialColor_decision, basic_para.normalColor] = deal([0.5 0 0 0.5]);%[1 0 0 0.5]
[specialColor, basic_para.specialColor] = deal([0.65 0 0]);%[0.5 0 0]-->[1 0 0]-->[0.8 0 0]
[frameColor, basic_para.frameColor] = deal([1 1 1]*0.45);%1-->0.4-->0.8
[frameColor2, basic_para.frameColor] = deal([1 1 1]*0.15);%0.25-->0.15
[frameColor3, basic_para.frameColor] = deal([1 1 1]*1);
[errorColor, basic_para.errorColor] = deal([0 0 1 0.75]);
[normalColor_offload, basic_para.normalColor_offload] = deal([1 1 1 0.75]);

[frameDisappearTime, basic_para.frameDisappearTime] = deal(200);

if ifSpecializedTraining_R == 1
    [taskMode_noChoice_Green0_Red1, basic_para.taskMode_noChoice_Green0_Red1] = deal(1);    
    [waterStepPer100, basic_para.waterStepPer100] = deal(0);%10-->8-->10-->12-->15-->12-->10-->8-->0
    if if_221False_329True == 0
        [initGiveWaterTime, basic_para.initGiveWaterTime] = deal(300);%140-->120
    elseif if_221False_329True == 1
        [initGiveWaterTime, basic_para.initGiveWaterTime] = deal(260);%120
    end

    if if_justFix0_touch1 == 0
        [initGiveWaterTime, basic_para.initGiveWaterTime] = deal(2);
        [seq_length_rangeHead, basic_para.seq_length_rangeHead] = deal(4);
        [seq_length_rangeTail, basic_para.seq_length_rangeTail] = deal(5);
    elseif if_justFix0_touch1 == 1
        [seq_length_rangeHead, basic_para.seq_length_rangeHead] = deal(1);
        [seq_length_rangeTail, basic_para.seq_length_rangeTail] = deal(2);%3
    end  
    
    %[numFrames_rangeHead, basic_para.numFrames_rangeHead] = deal(7);
    %[numFrames_rangeTail, basic_para.numFrames_rangeTail] = deal(7);
    [choiceCondition_boardLine1, basic_para.choiceCondition_boardLine1] = deal(0);
    [choiceCondition_boardLine2, basic_para.choiceCondition_boardLine2] = deal(1);
    [cumulativeCorrectLimit_seqLengthSwitch, basic_para.cumulativeCorrectLimit_seqLengthSwitch] = deal(1);%4-->6
    [consecutiveCorrectLimit, basic_para.consecutiveCorrectLimit] = deal(0);%2
    [cumulativeErrorLimit, basic_para.cumulativeErrorLimit] = deal(7 * ones(1, trial_num));%6-->7-->5-->3-->2
    [encourageCumuErrorLimit, basic_para.encourageCumuErrorLimit] = deal(5);%5-->10        
    
    [probFreeTouch, basic_para.probFreeTouch] = deal(0);%0.8, 0
    [ifFormalTask, basic_para.ifFormalTask] = deal(0);
    [if_seqList, basic_para.if_seqList] = deal(0);
    
    [if_length3_lessCCL, basic_para.if_length3_lessCCL] = deal(0);
    [if_length1_moreCCL, basic_para.if_length1_moreCCL] = deal(0);
     
    [if_validHoldUntilCorrect, basic_para.if_validHoldUntilCorrect] = deal(0);
    [if_R_mustCorrect, basic_para.if_R_mustCorrect] = deal(0);%1


end

if ifSpecializedTraining_G == 1
    [taskMode_noChoice_Green0_Red1, basic_para.taskMode_noChoice_Green0_Red1] = deal(0);
    [waterStepPer100, basic_para.waterStepPer100] = deal(0);%10-->8-->10-->12-->15-->12-->10-->8-->0
    if if_221False_329True == 0
        [initGiveWaterTime, basic_para.initGiveWaterTime] = deal(245);%240-->200
    elseif if_221False_329True == 1
        [initGiveWaterTime, basic_para.initGiveWaterTime] = deal(300);%340-->430-->420
    end
    [InternalMemoryWaterRatio, basic_para.InternalMemoryWaterRatio] = deal(1);
    %[seq_length_rangeHead, basic_para.seq_length_rangeHead] = deal(1);
    %[seq_length_rangeTail, basic_para.seq_length_rangeTail] = deal(1);
    
    if if_justFix0_touch1 == 0
        [initGiveWaterTime, basic_para.initGiveWaterTime] = deal(2);
        [seq_length_rangeHead, basic_para.seq_length_rangeHead] = deal(4);
        [seq_length_rangeTail, basic_para.seq_length_rangeTail] = deal(5);
    elseif if_justFix0_touch1 == 1
        [seq_length_rangeHead, basic_para.seq_length_rangeHead] = deal(1);
        [seq_length_rangeTail, basic_para.seq_length_rangeTail] = deal(1);%3
    end
    
    %[numFrames_rangeHead, basic_para.numFrames_rangeHead] = deal(6);
    %[numFrames_rangeTail, basic_para.numFrames_rangeTail] = deal(6);    
    [choiceCondition_boardLine1, basic_para.choiceCondition_boardLine1] = deal(1);
    [choiceCondition_boardLine2, basic_para.choiceCondition_boardLine2] = deal(1);
    [cumulativeCorrectLimit_seqLengthSwitch, basic_para.cumulativeCorrectLimit_seqLengthSwitch] = deal(1);%6
    [consecutiveCorrectLimit, basic_para.consecutiveCorrectLimit] = deal(0);%2
    [cumulativeErrorLimit, basic_para.cumulativeErrorLimit] = deal(6 * ones(1, trial_num));%6-->7-->5-->3-->2-->8
    [encourageCumuErrorLimit, basic_para.encourageCumuErrorLimit] = deal(5);%5-->10-->5-->4        
    
    [if_length4_freetouchCCL0, basic_para.if_length4_freetouchCCL0] = deal(0);
    
    [probFreeTouch, basic_para.probFreeTouch] = deal(0);%0.8, 0
    [ifFormalTask, basic_para.ifFormalTask] = deal(0);%0
    [if_seqList, basic_para.if_seqList] = deal(0);
    [if_breakWithNewSeq, basic_para.if_breakWithNewSeq] = deal(1);%0
    [if_cancelBreakForCCL, basic_para.if_cancelBreakForCCL] = deal(1);%0

    [if_length3_lessCCL, basic_para.if_length3_lessCCL] = deal(0);
    [if_length1_moreCCL, basic_para.if_length1_moreCCL] = deal(0);    
end

[ifFreeTouch, basic_para.ifFreeTouch] = deal(ifFreeTouch_origin * ones(1, trial_num) .* (randArray4<probFreeTouch));
ifFreeTouch_raw = ifFreeTouch;

trial_para.marker = struct;
trial_para.marker.content = [];
trial_para.marker.count = 0;
trial_para.marker.time = [];
%%%%%%%%%%%%%%%%%% basic adjustable parameters end %%%%%%%%%%%%%%%%%%
%% Some setup

% Setup for sending marker by mmep (20220710)
cd '/home/middle/STUDY/MatlabScript/SomeHardDriver/marker'
mSize = zeros(1,1);%这是你要用到的共享内存的大小，注意，必须大小实现敲定好
fileID =  fopen('m_marker.dat','w');   
fwrite(fileID, mSize, 'uint8');
fclose(fileID);
m_marker = memmapfile('m_marker.dat', 'Writable', true, 'Format', 'uint8'); %设置好共享内存

mSize = zeros(1,1);%这是你要用到的共享内存的大小，注意，必须大小实现敲定好
fileID =  fopen('m_markerCount.dat','w');   
fwrite(fileID, mSize, 'int32');
fclose(fileID);
m_markerCount = memmapfile('m_markerCount.dat', 'Writable', true, 'Format', 'int32'); %设置好共享内存
m_markerCount.Data = int32(trial_para.marker.count);



cd '/home/middle/STUDY/MatlabScript'

if if_eyeResponse == 0
    if joyMode_JJBJoy0_JoyMex1 == 0
        %use memmapfile to quickly get joystick state
        fileID =  fopen('mJoyState.dat','r+');
        % fileID =  fopen('mJoyState.dat','w');
        % mJoyState = zeros(1,1);%这是你要用到的共享内存的大小，注意，必须大小实现敲定好
        % fwrite(fileID, mJoyState, 'double');
        % fclose(fileID);
        m = memmapfile('mJoyState.dat', 'Writable', true, 'Format', 'double'); %设置好共享内存
    elseif joyMode_JJBJoy0_JoyMex1 == 1
        js = joyMex('/dev/input/js0');
        joyState = 0;
    end
end


if exist([pwd, '/MATDATA/'], 'dir') == 0
    mkdir MATDATA
end
filenum = 1;
while exist([pwd, '/MATDATA/', datestr(now, 'yyyy-mm-dd'), mfilename, '-', num2str(filenum), '.mat'], 'file') == 2
    filenum = filenum + 1;
end
filename = [pwd, '/MATDATA/', datestr(now, 'yyyy-mm-dd'), mfilename, '-', num2str(filenum), '.mat'];
shortfilename = [datestr(now, 'yyyy-mm-dd'), mfilename, '-', num2str(filenum), '.mat'];


%Save basic parameters
basic_para.mfilename = mfilename;
%save(filename,'basic_para');

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

if if_eyeResponse == 0
    dev = min(GetTouchDeviceIndices([], 1));
    fprintf('Touch device properties:\n');
    info = GetTouchDeviceInfo(dev);
    disp(info);
end

% Open an on screen window
% [window, windowRect] = PsychImaging('OpenWindow', screenNumber, [64/255 64/255 64/255]);
% [window, windowRect] = PsychImaging('OpenWindow', screenNumber, [0/255 0/255 0/255]);
% [window, windowRect] = PsychImaging('OpenWindow', screenNumber, [25.5/255 25.5/255 25.5/255]);
% [window, windowRect] = PsychImaging('OpenWindow', screenNumber, [0.15 0.15 0.15]);
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, bgColor);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% The avaliable keys to press
escapeKey = KbName('ESCAPE');
leftaltKey = KbName('LeftAlt');
applicationKey = KbName('Application');
% minusKey = KbName('-_');%giveWaterTime-
% equalKey = KbName('=+');%giveWaterTime+
% minusKey = KbName('-_');%preFixationTime-
% equalKey = KbName('=+');%preFixationTime+
minusKey = KbName('-_');%smallWaterInterval-
equalKey = KbName('=+');%smallWaterInterval+
leftBracketsKey = KbName('[{');%consecutiveCorrectLimit-
rightBracketsKey = KbName(']}');%consecutiveCorrectLimit+
orKey = KbName('\|');%giveWater function
% commaKey = KbName(',<');%cumulativeCorrectLimit_seqLengthSwitch-
% fullstopKey = KbName('.>');%cumulativeCorrectLimit_seqLengthSwitch+
% commaKey = KbName(',<');%twoFixationTime-
% fullstopKey = KbName('.>');%twoFixationTime+
commaKey = KbName(',<');%numFrames-
fullstopKey = KbName('.>');%numFrames+
% commaKey = KbName(',<');%freeTouch_RT_limit-
% fullstopKey = KbName('.>');%freeTouch_RT_limit+
% commaKey = KbName(',<');%probFreeTouch-
% fullstopKey = KbName('.>');%probFreeTouch+
% questionKey = KbName('/?');%mannul lever
questionKey = KbName('/?');%whether disable right selection
nineKey = KbName('9(');%fixWinSize-
zeroKey = KbName('0)');%fixWinSize+
% nineKey = KbName('9(');%stimulusExtendTime-
% zeroKey = KbName('0)');%stimulusExtendTime+
% nineKey = KbName('9(');%fixationTime--
% zeroKey = KbName('0)');%fixationTime-+
% nineKey = KbName('9(');%InternalMemoryWaterRatio-
% zeroKey = KbName('0)');%InternalMemoryWaterRatio+

%fixationTime
reCenterKey = KbName('r');%reCenter fix point according to current gaze position
tKey = KbName('t');%reset fix point
% oKey = KbName('o');%dynamicFixThreshold-
% pKey = KbName('p');%dynamicFixThreshold+
% oKey = KbName('o');%fixTolerenceTime-
% pKey = KbName('p');%fixTolerenceTime+
oKey = KbName('o');%submitHoldingTime-
pKey = KbName('p');%submitHoldingTime+
% oKey = KbName('o');%InternalMemoryWaterRatio-
% pKey = KbName('p');%InternalMemoryWaterRatio+
% oKey = KbName('o');%pointShowTime-
% pKey = KbName('p');%pointShowTime+
% oKey = KbName('o');%endBet_fastSlowBoundary-
% pKey = KbName('p');%endBet_fastSlowBoundary+
% zKey = KbName('z');%propExplore-
% xKey = KbName('x');%propExplore+
% zKey = KbName('z');%cumulativeCorrectLimit_greenRedSwitch-
% xKey = KbName('x');%cumulativeCorrectLimit_greenRedSwitch+
% zKey = KbName('z');%choiceCondition_boardLine1-
% xKey = KbName('x');%choiceCondition_boardLine1+
% cKey = KbName('c');%choiceCondition_boardLine2-
% vKey = KbName('v');%choiceCondition_boardLine2+
xKey = KbName('x');%adjust rFPW evenly
cKey = KbName('c');%adjust rFPW automatically
vKey = KbName('v');%display rFPW 
% qKey = KbName('q');%probability of first frame is redFill+
% qKey = KbName('q');%probability of last frame is redFill+
% wKey = KbName('w');%probability of second frame is redFill+
% eKey = KbName('e');%probability of third frame is redFill+
% rKey = KbName('r');%probability of fourth frame is redFill+
% tKey = KbName('t');%probability of fifth frame is redFill+
qKey = KbName('q');%reCenterGazeBias_x+
wKey = KbName('w');%reCenterGazeBias_x-
sKey = KbName('s');%reCenterGazeBias_y+
dKey = KbName('d');%reCenterGazeBias_y-
%sKey = KbName('s');%switch frame-showing mode
% sKey = KbName('s');%switch between force-choice-g, force-choice-r and free-choice-g/r
% dKey = KbName('d');%switch between internal memory task and offloading task
lKey = KbName('l');%switch ifLinkTriangleSquare mode
mKey = KbName('m');%switch if_mouse_mannual_eyeResponse


% Screen positions of our  rectangles
squareXpos = ones(1, numSquares); 
squareYpos = ones(1, numSquares); 

%Arc-like points metrix
theta = zeros(1, numSquares);
for i = 1:numSquares
    %theta(i) = (i-1)/(numSquares-1) * 330/360 * pi + 15/360 * pi;
    theta(i) = (i-1)/(numSquares-1) * (180-2*InitAngle_arc)/180 * pi + InitAngle_arc/180 * pi;
    squareXpos(i) = -1 * radius * cos(theta(i));
    squareYpos(i) = tan(theta(i)) * squareXpos(i) + screenYpixels* 0.7;
    squareXpos(i) = squareXpos(i) + screenXpixels/2;
end


squareXpos = squareXpos + screenXpixels_bias;
squareYpos = squareYpos + screenYpixels_bias;

originalPointXpos = screenXpixels/2 + screenXpixels_bias;
originalPointYpos = screenYpixels* 0.7 + screenYpixels_bias;

% Make our rectangle coordinates
showRects = nan(4, numSquares);
for i = 1:numSquares
    showRects(:, i) = CenterRectOnPointd( showBaseRect, squareXpos(i), squareYpos(i) );
end
fixupdateTimeBound_fold;
fixupdateTimeFoldNum;
showRects_fold = nan(fixupdateTimeFoldNum,4, numSquares);
for tempk=1:fixupdateTimeFoldNum
    temp_scale = tempk/(1+fixupdateTimeFoldNum);
    for i = 1:numSquares
        showRects_fold(tempk,:, i) = CenterRectOnPointd( (showBaseRect+10)*temp_scale, squareXpos(i), squareYpos(i) );
    end
end
% tempk = 1;

detectRects = nan(4, numSquares);
for i = 1:numSquares
    if i == 1
        tempDetectBaseRect = [0 0 detectBaseRect(3:4)+firstDetectRectModifed];
        detectRects(:, i) = CenterRectOnPointd( tempDetectBaseRect, squareXpos(i), squareYpos(i) );
    else
        detectRects(:, i) = CenterRectOnPointd( detectBaseRect, squareXpos(i), squareYpos(i) );
    end
end
%showTriangles = nan(8, numSquares);
showTriangles = cell(1, numSquares);
anglesDeg = linspace(-90, 270, 3 + 1);
anglesRad = anglesDeg * (pi / 180);
for i = 1:numSquares
    tri_yPosVector = sin(anglesRad) .* showBaseTriRadius + squareYpos(i);
    tri_xPosVector = cos(anglesRad) .* showBaseTriRadius + squareXpos(i);
    showTriangles(i) = {[tri_xPosVector; tri_yPosVector]'};
end


% Screen positions of selection rectangles
%selectionXpos = [screenXpixels*9/20 screenXpixels*11/20];
%selectionXpos = [screenXpixels*0.4625 screenXpixels*0.5375];
% selectionXpos = [screenXpixels*0.455 screenXpixels*0.545]+screenXpixels_bias;%0.44-->0.45-->0.455
selectionXpos = [screenXpixels*0.43 screenXpixels*0.57]+screenXpixels_bias;%0.44-->0.45-->0.455-->0.43


% selectionYpos = [screenYpixels*selectionYpos_coefficient screenYpixels*selectionYpos_coefficient];
% selectionYpos = [screenYpixels/2+screenYpixels_bias+205-20 screenYpixels/2+screenYpixels_bias+205-20];
selectionYpos = [screenYpixels/2+screenYpixels_bias+215 screenYpixels/2+screenYpixels_bias+215];

selectionTriangles = cell(1, 2);
anglesDeg = linspace(-90, 270, 3 + 1);
anglesRad = anglesDeg * (pi / 180);
for i = 1:2
    %select_tri_yPosVector = sin(anglesRad) .* selectBaseTriRadius + selectionYpos(i) + 20;
    select_tri_yPosVector = sin(anglesRad) .* selectBaseTriRadius + selectionYpos(i) + 20;
    select_tri_xPosVector = cos(anglesRad) .* selectBaseTriRadius + selectionXpos(i);
    selectionTriangles(i) = {[select_tri_xPosVector; select_tri_yPosVector]'};
end


% Make selection rectangle coordinates
selectionRects = nan(4, 2);
for i = 1:2
    selectionRects(:, i) = CenterRectOnPointd(selectbaseRect, selectionXpos(i), selectionYpos(i));
end

%selectbaseRect_detection
% selectionXpos_detection(1) = selectionXpos(1) + selectbaseRect(3)/2 - selectbaseRect_detection(3)/2;
% selectionXpos_detection(2) = selectionXpos(2) - selectbaseRect(3)/2 + selectbaseRect_detection(3)/2;
selectionXpos_detection(1) = selectionXpos(1);
selectionXpos_detection(2) = selectionXpos(2);
selectionYpos_detection = selectionYpos;

selectionRects_detection = nan(4, 2);
for i = 1:2
    selectionRects_detection(:, i) = CenterRectOnPointd(selectbaseRect_detection,...
        selectionXpos_detection(i), selectionYpos_detection(i));
end


if ifSelectRandom == 1
    isGreenRectLeft = [zeros(1, floor(trial_num/2)) ones(1, trial_num-floor(trial_num/2))];
    isGreenRectLeft = isGreenRectLeft(randperm(trial_num));
else
    isGreenRectLeft = ones(1, trial_num);    
end
trial_para.isGreenRectLeft = isGreenRectLeft;


submitXpos = screenXpixels*0.5+screenXpixels_bias;
% submitYpos = screenYpixels*0.5+screenYpixels_bias+205-20;
% submitYpos = screenYpixels*0.5+screenYpixels_bias+205-20+100;
% submitYpos = screenYpixels*0.5+screenYpixels_bias+205-20+150;
submitYpos = screenYpixels/2+screenYpixels_bias+215;

% Make submit rectangle coordinates
submitRects = CenterRectOnPointd(submitBaseRect, submitXpos, submitYpos);

%submitBaseRect_detection
submitRects_dection = CenterRectOnPointd(submitBaseRect_detection, submitXpos, submitYpos);


holdXpos = screenXpixels*0.5+screenXpixels_bias;
holdYpos = screenYpixels*0.5+screenYpixels_bias+205-20+400;

% Make hold rectangle coordinates
holdRects = CenterRectOnPointd(holdBaseRect, holdXpos, holdYpos);

%holdBaseRect_detection
holdRects_dection = CenterRectOnPointd(holdBaseRect_detection, holdXpos, holdYpos);



% Now we set the coordinates (these are all relative to zero we will let
% the drawing routine center the cross in the center of our monitor for us)
xCross = [-fixCrossDimPix fixCrossDimPix 0 0];
yCross = [0 0 -fixCrossDimPix fixCrossDimPix];
xyCross = [xCross; yCross];
%crossCenter = [screenXpixels/2+screenXpixels_bias screenYpixels/2+screenYpixels_bias+75+80];%+75 is new, +80
% crossCenter = [screenXpixels/2 screenYpixels/2-105];%-75-->-105
% crossCenter = [screenXpixels/2+screenXpixels_bias screenYpixels/2+screenYpixels_bias+205-20-25];
% crossCenter = [screenXpixels/2+screenXpixels_bias screenYpixels/2+screenYpixels_bias+205-20-25+30];
crossCenter = [screenXpixels/2+screenXpixels_bias screenYpixels/2+screenYpixels_bias+215];%+190-->+200

crossVerticalRect = [0 0 crossWidth crossLength];%20-->12
crossHorizontalRect = [0 0 crossLength crossWidth];

centeredVerticalCrossRect = CenterRectOnPointd(crossVerticalRect, crossCenter(1), crossCenter(2));
centeredHorizontalCrossRect = CenterRectOnPointd(crossHorizontalRect, crossCenter(1), crossCenter(2));

fixSquare = [0 0 fixSquareSize fixSquareSize];
centeredFixSquare = CenterRectOnPointd(fixSquare, crossCenter(1), crossCenter(2));

fixTri_yPosVector = sin(anglesRad) .* fixTriRadius + crossCenter(2);
fixTri_xPosVector = cos(anglesRad) .* fixTriRadius + crossCenter(1);
%Some flags initialization
exitProgram = 0;
isManual = 0;
isManual_last = 0;
oldState = 0;
isReleaseToHold = 0;
isHoldToRelease = 0;
isManual_rFPW = 1;
frameMode = 0;
choiceMode = 2;
choiceMode_backup = 2;
choiceCondition_flag = zeros(1, trial_num);
%ifLinkTriangleSquare = 0;
isTriangletrialMiddle = 0;

trial_count = 1;
uniqueSeq_trial_count = 0;
internal_trial_count = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);
offloading_trial_count = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);
internal_trial_count_newSeq = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);
offloading_trial_count_newSeq = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);
internal_trial_count_last20NewSeq = zeros(1, 20+1);
offloading_trial_count_last20NewSeq = zeros(1, 20+1);
internal_correct_count = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);
offloading_correct_count = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);
internalAccuracy = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);
offloadingAccuracy = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);
offloadingProb = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);

internal_noChoice_trial_count = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);
internal_noChoice_correct_count = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);
internalNoChoiceAccuracy = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);


break_sumCount = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);
seqLength_trial_count = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);
breakRate = zeros(1,seq_length_rangeTail-seq_length_rangeHead+1+1);

%break_count = 0;
break_count = zeros(1, trial_num);
showing_count = 0;
isCorrect = zeros(1, trial_num);
isIppatsuCorrect = zeros(1, trial_num);
sort_correctCount = zeros(1, numFrames_rangeTail);
sort_allCount = zeros(1, numFrames_rangeTail);
sort_isCorrect = zeros(1, numFrames_rangeTail);
sort_correctCount_last40 = zeros(40, numFrames_rangeTail);
sort_allCount_last40 = zeros(40, numFrames_rangeTail);
isWaitGocueBreak = -1*ones(1, trial_num);
isTouchFrameCorrect = zeros(1, trial_num);
isTouchingTimeOut = zeros(1, trial_num);
isBackTouch = zeros(1, trial_num);
consecutiveCorrect = 0;
cumulativeError = 0;
isProbWeightChange = 0;
isProbWeightEven = 0;
isRelaSeqChange = 0;
tempWaterTime = initGiveWaterTime;
cumulativeCorrect_greenRedSwitch = 0;
mytouchQueueBuffer = cell(1, touchBufferSize);
mytouchQueueBuffer_selection = cell(1, touchBufferSize);
mytouchQueueBuffer_submit = cell(1, touchBufferSize);

newTouchEnableFlag = 1;
selectFlag_newSeq = zeros(1, trial_num);

preFixationEndFlag = 0;
last_preFixationEndFlag = 0;
initialFlag = 0;

t0_showReady = tic;

if taskMode_noChoice_Green0_Red1 == 0
    last_ifSelectOffloading = 0;        
elseif taskMode_noChoice_Green0_Red1 == 1
    last_ifSelectOffloading = 1;        
end

fixOn_tolerenceFlag = 0;
currentTrial_fixPos = struct;
currentTrial_fixPos.x = [];
currentTrial_fixPos.y = [];
gazeBias_offset_count = 0;

fake_isHold = 0;

if_disableRightSelection = 0;

isShowingInitial = 1;
isShowing = 0;
isDelay = 0;
isReleaseForTouch = 0;
isSelecting = 0;
isTwoHolding = 0;
ifTwoHoldFlag = 0;
isTwoDelay = 0;
isTouching = 0;
isItiInitial = 0;
isIti = 0;
isItiInitial_break = 0;

%% Get seq_length and numFrames
%cumulativeCorrectLimit_seqLengthSwitch
if ifSeq_length_3to1 == 1
    temp_seq_length = [1 4];
    temp_seq_length2 = cell(1, seq_length_rangeTail-seq_length_rangeHead+1);
    for i=1:seq_length_rangeTail-seq_length_rangeHead+1
        temp_seq_length2{i} = repmat(temp_seq_length(i), 1, cumulativeCorrectLimit_seqLengthSwitch);
    end
elseif ifSeq_length_3to1 == 0 
    temp_seq_length = seq_length_rangeHead:seq_length_rangeTail;
    temp_seq_length2 = cell(1, seq_length_rangeTail-seq_length_rangeHead+1);
    for i=1:seq_length_rangeTail-seq_length_rangeHead+1
        temp_seq_length2{i} = repmat(temp_seq_length(i), 1, cumulativeCorrectLimit_seqLengthSwitch);
    end
end

% random the blockIndex every X blocks
period = 4;
period_seq = period*length(temp_seq_length);

if ifSeqLength14_sparse == 1 && seq_length_rangeHead == 1 && seq_length_rangeTail == 4
    temp_seq_length_sparse = temp_seq_length(2:end-1);
    temp_sparse = [temp_seq_length temp_seq_length_sparse];
    period_seq = period_seq - 2*(period/2);
    temp_blockIndex = repmat(temp_sparse, 1, ceil(trial_num/cumulativeCorrectLimit_seqLengthSwitch/length(temp_sparse)));   
elseif ifSeqLength1_sparse == 1 && seq_length_rangeHead == 1
    temp_seq_length_sparse = temp_seq_length(2:end);
    temp_sparse = [temp_seq_length temp_seq_length_sparse];
    period_seq = period_seq - period/2;
    temp_blockIndex = repmat(temp_sparse, 1, ceil(trial_num/cumulativeCorrectLimit_seqLengthSwitch/length(temp_sparse)));       
else
    temp_blockIndex = repmat(temp_seq_length, 1, ceil(trial_num/cumulativeCorrectLimit_seqLengthSwitch/length(temp_seq_length)));
end



for i=1:length(temp_blockIndex)/period_seq
    temp_blockIndex(1+(i-1)*period_seq:i*period_seq) = temp_blockIndex(randperm(period_seq)+(i-1)*period_seq);
end

seq_length = zeros(1, trial_num);
for i=1:length(temp_blockIndex)
    for j=1:cumulativeCorrectLimit_seqLengthSwitch
        seq_length((i-1)*cumulativeCorrectLimit_seqLengthSwitch+j) = temp_seq_length2{find(temp_seq_length==temp_blockIndex(i))}(j); %#ok<*FNDSB>
    end
end


temp_numFrames = numFrames_rangeHead:numFrames_rangeTail;
temp_numFrames2 = repmat(temp_numFrames,1,ceil((trial_num)/(length(temp_numFrames))));
numFrames = temp_numFrames2(randperm((trial_num)));


trial_para.seq_length = seq_length;
trial_para.numFrames = numFrames;

%% Initialization for redFillProbWeight & relativeSequence
temprand = rand(max(seq_length), trial_num);
redFillProbWeight = cell(1, trial_num);
for temp_trial_count=1:trial_num
    redFillProbWeight{temp_trial_count} = ones(1, numFrames(temp_trial_count)) / numFrames(temp_trial_count);
end
relativeSequence = cell(1, trial_num);

if if_seqList == 1 && ifFormalTask == 1
    relativeSequence = [];
    choiceCondition_flag_list = [];
    temp_seq_length_rangeHead = seq_length_rangeHead;
    
    randForLength4 = 1:nchoosek(numFrames_rangeHead,4);
    
    if if_length1_alone == 1
        temp_seq_length_rangeHead = 2;
        [seq_list_oneRun, choiceCondition_flag_oneRun, randForLength4] = ...
            get_seq_list_D9(temp_seq_length_rangeHead, seq_length_rangeTail, numFrames_rangeTail, randForLength4);
    else
        [seq_list_oneRun, choiceCondition_flag_oneRun, randForLength4] = ...
            get_seq_list_D9(seq_length_rangeHead, seq_length_rangeTail, numFrames_rangeTail, randForLength4);
    end
    temp_trialNum_oneRun = length(choiceCondition_flag_oneRun);
    %temp_runNum = ceil(trial_num/temp_trialNum_oneRun);
    
    temp_ok_flag = 0;
    temp_length1_flag = 0;
    temp_length1_runCount = 0;
    trialNum_length1Run = 0;
    while temp_ok_flag == 0
        if if_length1_alone == 1 && temp_length1_flag == 0
            temp_length1_runCount = temp_length1_runCount + 1;
            if temp_length1_runCount > length1_runNum
                temp_length1_flag = 1;
                trialNum_length1Run = length(choiceCondition_flag_list);
                continue
            else
                [seq_list_oneRun, choiceCondition_flag_oneRun] = ...
                    get_seq_list_D9(1, 1, numFrames_rangeTail); %#ok<*ASGLU>
                a = 1;
            end
        else
            [seq_list_oneRun, choiceCondition_flag_oneRun, randForLength4] = ...
                get_seq_list_D9(temp_seq_length_rangeHead, seq_length_rangeTail, numFrames_rangeTail, randForLength4);
        end

        relativeSequence = [relativeSequence seq_list_oneRun]; %#ok<AGROW>
        choiceCondition_flag_list = [choiceCondition_flag_list choiceCondition_flag_oneRun]; %#ok<AGROW>
        if length(choiceCondition_flag_list) >= trial_num
            temp_ok_flag = 1;
        end
    end
    relativeSequence = relativeSequence(1, 1:trial_num);
    choiceCondition_flag_list = choiceCondition_flag_list(1, 1:trial_num);
    
    for tempi=1:trial_num
        seq_length(tempi) = length(relativeSequence{tempi});
    end
    trial_para.seq_length = seq_length;
else
    if max(seq_length) <= 2
        
        for temp_trial_count=1:trial_num
            for seq_index=1:seq_length(temp_trial_count)
                for i=1:numFrames(temp_trial_count)
                    if i==1
                        if temprand(seq_index, temp_trial_count) < redFillProbWeight{temp_trial_count}(1)
                            relativeSequence{temp_trial_count}(seq_index) = 1;
                        end
                    else
                        if temprand(seq_index, temp_trial_count) > sum(redFillProbWeight{temp_trial_count}(1:i-1))...
                                && temprand(seq_index, temp_trial_count) < sum(redFillProbWeight{temp_trial_count}(1:i))
                            relativeSequence{temp_trial_count}(seq_index) = i;
                        end
                    end
                end
                %Change same red fill so that they are unique
                while length(relativeSequence{temp_trial_count})-length(unique(relativeSequence{temp_trial_count})) > 0
                    %if length(relativeSequence{temp_trial_count})-length(unique(relativeSequence{temp_trial_count})) > 0
                    tempValue = relativeSequence{temp_trial_count}(end);
                    seg_temprand = temprand(seq_index, temp_trial_count);
                    if tempValue > 1
                        seg_temprand = seg_temprand - sum(redFillProbWeight{temp_trial_count}(1:tempValue-1));
                    end
                    tempValueWeight = redFillProbWeight{temp_trial_count}(tempValue);
                    for i=1:numFrames(temp_trial_count)
                        if i < tempValue
                            if i == 1
                                if seg_temprand/tempValueWeight < ...
                                        redFillProbWeight{temp_trial_count}(1)/(1-tempValueWeight)
                                    relativeSequence{temp_trial_count}(end) = 1;
                                end
                            else
                                if seg_temprand/tempValueWeight > ...
                                        sum(redFillProbWeight{temp_trial_count}(1:i-1))/(1-tempValueWeight)...
                                        && seg_temprand/tempValueWeight <...
                                        sum(redFillProbWeight{temp_trial_count}(1:i))/(1-tempValueWeight)
                                    relativeSequence{temp_trial_count}(end) = i;
                                end
                            end
                        elseif i > tempValue
                            if seg_temprand/tempValueWeight > ...
                                    (sum(redFillProbWeight{temp_trial_count}(1:i-1))-tempValueWeight)/(1-tempValueWeight)...
                                    && seg_temprand/tempValueWeight <...
                                    (sum(redFillProbWeight{temp_trial_count}(1:i))-tempValueWeight)/(1-tempValueWeight)
                                relativeSequence{temp_trial_count}(end) = i;
                            end
                            
                        end
                    end
                end
            end
        end
        
    elseif max(seq_length) >= 3
        for temp_trial_count=1:trial_num
            relativeSequence{temp_trial_count} = randperm(numFrames(temp_trial_count), seq_length(temp_trial_count));
        end
    end
    
    %temporary code, transform seq_length == 3 to seq_length == 1
    if ifSeq_length_3to1 == 1
        for i=1:trial_num
            if length(relativeSequence{i}) == 3
                relativeSequence{i} = relativeSequence{i}(1);
                seq_length(i) = 1;
            end
        end
    end
    
    
    for i=1:trial_num
        a(i) = length(relativeSequence{i})-length(unique(relativeSequence{i})); %#ok<SAGROW>
    end
    
    
end

randExplore = rand(1, trial_num);

%% Eyelink preparation
if true
    %% Eyelink STEP 1: INITIALIZE EYELINK CONNECTION; OPEN EDF FILE; GET EYELINK TRACKER VERSION
    %commandwindow;
    % ListenChar(-1);
    % Initialize EyeLink connection (dummymode = 0) or run in "Dummy Mode" without an EyeLink connection (dummymode = 1);
    EyelinkInit(dummymode); % Initialize EyeLink connection
    status = Eyelink('IsConnected');
    if status < 1 % If EyeLink not connected
        dummymode = 1;
    end
    
    edfFile = 'demo3';
    % edfFile = [datestr(now, 'yyyy-mm-dd'), mfilename, '-', num2str(filenum)];
    
    % Open an EDF file and name it
    failOpen = Eyelink('OpenFile', edfFile);
    if failOpen ~= 0 % Abort if it fails to open
        fprintf('Cannot create EDF file %s', edfFile); % Print some text in Matlab's Command Window
        % Cleanup function used throughout the script above
        try %#ok<TRYNC>
            Screen('CloseAll'); % Close window if it is open
        end
        Eyelink('Shutdown'); % Close EyeLink connection
        % ListenChar(0); % Restore keyboard output to Matlab
        return
    end
    
    % Get EyeLink tracker and software version
    % <ver> returns 0 if not connected
    % <versionstring> returns 'EYELINK I', 'EYELINK II x.xx', 'EYELINK CL x.xx' where 'x.xx' is the software version
    ELsoftwareVersion = 0; % Default EyeLink version in dummy mode
    [ver, versionstring] = Eyelink('GetTrackerVersion');
    if dummymode == 0 % If connected to EyeLink
        % Extract software version number.
        [r1, vnumcell] = regexp(versionstring,'.*?(\d)\.\d*?','Match','Tokens'); % Extract EL version before decimal point
        ELsoftwareVersion = str2double(vnumcell{1}{1}); % Returns 1 for EyeLink I, 2 for EyeLink II, 3/4 for EyeLink 1K, 5 for EyeLink 1KPlus, 6 for Portable Duo
        % Print some text in Matlab's Command Window
        fprintf('Running experiment on %s version %d\n', versionstring, ver );
    end
    % Add a line of text in the EDF file to identify the current experimemt name and session. This is optional.
    % If your text starts with "RECORDED BY " it will be available in DataViewer's Inspector window by clicking
    % the EDF session node in the top panel and looking for the "Recorded By:" field in the bottom panel of the Inspector.
    preambleText = sprintf('RECORDED BY Psychtoolbox demo %s session name: %s', mfilename, edfFile);
    Eyelink('Command', 'add_file_preamble_text "%s"', preambleText);
    
    %% Eyelink STEP 2: SELECT AVAILABLE SAMPLE/EVENT DATA
    % See EyeLinkProgrammers Guide manual > Useful EyeLink Commands > File Data Control & Link Data Control
    
    % Select which events are saved in the EDF file. Include everything just in case
    Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    % Select which events are available online for gaze-contingent experiments. Include everything just in case
    Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,BUTTON,FIXUPDATE,INPUT');
    % Select which sample data is saved in EDF file or available online. Include everything just in case
    if ELsoftwareVersion > 3  % Check tracker version and include 'HTARGET' to save head target sticker data for supported eye trackers
        Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,RAW,AREA,HTARGET,GAZERES,BUTTON,STATUS,INPUT');
        Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT');
    else
        Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,RAW,AREA,GAZERES,BUTTON,STATUS,INPUT');
        Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    end
    
    %% Eyelink STEP 3: Create central square fixation window
    fixationBase_ELHost = [-fixWinSize_ELHost -fixWinSize_ELHost ...
        fixWinSize_ELHost fixWinSize_ELHost];
    fixationRect_ELHost = CenterRectOnPointd(fixationBase_ELHost, crossCenter(1), crossCenter(2));
    
    %% Eyelink STEP 4: SET CALIBRATION SCREEN COLOURS; PROVIDE WINDOW SIZE TO EYELINK HOST & DATAVIEWER; SET CALIBRATION PARAMETERS; CALIBRATE
    
    % Provide EyeLink with some defaults, which are returned in the structure "el".
    el = EyelinkInitDefaults(window);
    % set calibration/validation/drift-check(or drift-correct) size as well as background and target colors.
    % It is important that this background colour is similar to that of the stimuli to prevent large luminance-based
    % pupil size changes (which can cause a drift in the eye movement data)
    el.calibrationtargetsize = 3;% Outer target size as percentage of the screen
    el.calibrationtargetwidth = 0.7;% Inner target size as percentage of the screen
    el.backgroundcolour = [0 0 0];% [128 128 128]
    el.calibrationtargetcolour = [256 256 256];% [0 0 0]
    % set "Camera Setup" instructions text colour so it is different from background colour
    el.msgfontcolour = [128 128 128];% [0 0 0]
    
    % You must call this function to apply the changes made to the el structure above
    EyelinkUpdateDefaults(el);
    
    %% Eyelink STEP 5: Some other preparations
    % STEP 5.2: START TRIAL; SHOW TRIAL INFO ON HOST PC; SHOW BACKDROP IMAGE AND/OR DRAW FEEDBACK GRAPHICS ON HOST PC; DRIFT-CHECK/CORRECTION
    % Supply the trial number as a line of text on Host PC screen
    %Eyelink('Command', 'record_status_message "TRIAL %d/%d"', i, length(imgList));
    
    % Draw graphics on the EyeLink Host PC display. See COMMANDS.INI in the Host PC's exe folder for a list of commands
    Eyelink('SetOfflineMode'); % Put tracker in idle/offline mode before drawing Host PC graphics and before recording
    Eyelink('Command', 'clear_screen 0'); % Clear Host PC display from any previus drawing
    
    % Optional: draw feedback box and lines on Host PC interface instead of (or on top of) backdrop image.
    % See section 25.7 'Drawing Commands' in the EyeLink Programmers Guide manual
    Eyelink('Command', 'draw_filled_box %d %d %d %d 15', fixationRect_ELHost(1), fixationRect_ELHost(2), fixationRect_ELHost(3), fixationRect_ELHost(4)); % Fixation window
    
    tempSquareXpos = round(squareXpos);
    tempSquareYpos = round(squareYpos);
    
    for tempi=1:numSquares
        Eyelink('Command', 'draw_filled_box %d %d %d %d 15',...
            tempSquareXpos(tempi)-20, tempSquareYpos(tempi)-20,...
            tempSquareXpos(tempi)+20, tempSquareYpos(tempi)+20);
    end
    Eyelink('Command', 'draw_cross %d %d 0', crossCenter(1), crossCenter(2)); % Central crosshairs
    
    %     % Perform a drift check/correction.
    %     % Optionally provide x y target location, otherwise target is presented on screen centre
    %     EyelinkDoDriftCorrection(el, round(screenXpixels/2), round(screenYpixels/2));
    
    %STEP 5.3: START RECORDING
    
    Eyelink('SetOfflineMode');% Put tracker in idle/offline mode before recording
    Eyelink('StartRecording'); % Start tracker recording
    WaitSecs(0.1); % Allow some time to record a few samples before presenting first stimulus
    
    % STEP 5.4: PRESENT STIMULUS; CREATE DATAVIEWER BACKDROP AND INTEREST AREA
    
    % Check which eye is available for gaze-contingent drawing. Returns 0 (left), 1 (right) or 2 (binocular)
    eyeUsed = Eyelink('EyeAvailable');
    % Get events from right eye if binocular
    if eyeUsed == 2
        eyeUsed = 1;
    end
%     reCenterGazeBias_x = 0;
%     reCenterGazeBias_y = 0;
%     ELHost_changeFLag = 0;    
    %reCenterGazeBias_x = -44.3286;
    %reCenterGazeBias_y = -222.6857;
    
    MAT_file=dir([searchName_reCenterGazeBias, '*']);
    if isempty(MAT_file) == 0
        path = pwd;
        load_reCenterGazeBias = loadMat_single(searchName_reCenterGazeBias, path);
        reCenterGazeBias_x = load_reCenterGazeBias.reCenterGazeBias_x;
        reCenterGazeBias_y = load_reCenterGazeBias.reCenterGazeBias_y;
    else
        path = pwd;
        reCenterGazeBias_x = -44.3286;
        reCenterGazeBias_y = -222.6857;
    end    
    
    ELHost_changeFLag = 1;
    
end
%%
% Sync us and get a time stamp
vbl = Screen('Flip', window);
shortWaitFrames = shortShowTime * 60 /1000;

if if_freeview0_fixation1 == 1
    t0_fixOn = tic;
end

temp_consecutiveCorrectLimit = consecutiveCorrectLimit;

if if_eyeResponse == 0
    TouchQueueCreate(window, dev);
    TouchQueueStart(dev);
end

% Maximum priority level
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);
%% While loop
% Loop the animation until the escape key is pressed
while exitProgram == 0  
    %% Keyboard check
    % Check the keyboard to see if a button has been pressed
    [keyIsDown,secs, keyCode] = KbCheck;          
    
    if keyCode(leftaltKey) == 1 && keyCode(qKey) == 1 && lastKeyCode(qKey) == 0
        reCenterGazeBias_x = reCenterGazeBias_x + 5;
        fprintf("reCenterGazeBias_x = %d\n",reCenterGazeBias_x);
        ELHost_changeFLag = 1;
    end    
    if keyCode(leftaltKey) == 1 && keyCode(wKey) == 1 && lastKeyCode(wKey) == 0
        reCenterGazeBias_x = reCenterGazeBias_x - 5;
        fprintf("reCenterGazeBias_x = %d\n",reCenterGazeBias_x);
        ELHost_changeFLag = 1;
    end     
    if keyCode(leftaltKey) == 1 && keyCode(sKey) == 1 && lastKeyCode(sKey) == 0
        reCenterGazeBias_y = reCenterGazeBias_y + 5;
        fprintf("reCenterGazeBias_y = %d\n",reCenterGazeBias_y);
        ELHost_changeFLag = 1;
    end    
    if keyCode(leftaltKey) == 1 && keyCode(dKey) == 1 && lastKeyCode(dKey) == 0
        reCenterGazeBias_y = reCenterGazeBias_y - 5;
        fprintf("reCenterGazeBias_y = %d\n",reCenterGazeBias_y);
        ELHost_changeFLag = 1;
    end      
    
    if keyCode(leftaltKey) == 1 && keyCode(tKey) == 1 && lastKeyCode(tKey) == 0
        reCenterGazeBias_x = 0;
        reCenterGazeBias_y = 0;
        fprintf('Now is reset!\n');
        ELHost_changeFLag = 1;
    end    
    if keyCode(leftaltKey) == 1 && keyCode(reCenterKey) == 1 && lastKeyCode(reCenterKey) == 0
%         if if_freeview0_fixation1 == 1
            reCenterGazeBias_x = crossCenter(1) - elEvt.gx;
            reCenterGazeBias_y = crossCenter(2) - elEvt.gy;
%         else
%             reCenterGazeBias_x = 0;
%             reCenterGazeBias_y = 0;
%         end
        fprintf('Now is reCenter!\n');
        ELHost_changeFLag = 1;
    end
    if keyCode(leftaltKey) == 1 && keyCode(nineKey) == 1 && lastKeyCode(nineKey) == 0
        if fixWinSize_DisplayPC > 30
            fixWinSize_DisplayPC = fixWinSize_DisplayPC - 10;
        else
            fixWinSize_DisplayPC = fixWinSize_DisplayPC - 3;
        end
        if fixWinSize_DisplayPC < 0
            fixWinSize_DisplayPC = 0;
        end
        fprintf("fixWinSize = %d\n",fixWinSize_DisplayPC);
        fixWinSize_ELHost = fixWinSize_DisplayPC;
        fixationBase_ELHost = [-fixWinSize_ELHost -fixWinSize_ELHost ...
            fixWinSize_ELHost fixWinSize_ELHost];
        ELHost_changeFLag = 1;
    end
    if keyCode(leftaltKey) == 1 && keyCode(zeroKey) == 1 && lastKeyCode(zeroKey) == 0
        if fixWinSize_DisplayPC > 30
            fixWinSize_DisplayPC = fixWinSize_DisplayPC + 10;
        else
            fixWinSize_DisplayPC = fixWinSize_DisplayPC + 3;
        end
        fprintf("fixWinSize = %d\n",fixWinSize_DisplayPC);
        fixWinSize_ELHost = fixWinSize_DisplayPC;
        fixationBase_ELHost = [-fixWinSize_ELHost -fixWinSize_ELHost ...
            fixWinSize_ELHost fixWinSize_ELHost];
        ELHost_changeFLag = 1;
    end
    
    if ELHost_changeFLag == 1
        tempWinPos = [crossCenter(1)-reCenterGazeBias_x ...
            crossCenter(2)-reCenterGazeBias_y];
        tempWinPos = round(tempWinPos);        
        fixationRect_ELHost = CenterRectOnPointd(fixationBase_ELHost, tempWinPos(1), tempWinPos(2));
        
        for tempi=1:4
            if fixationRect_ELHost(tempi) < 0
                fixationRect_ELHost(tempi) = 0;
            end
        end
        if fixationRect_ELHost(1) > screenXpixels
            fixationRect_ELHost(1) = screenXpixels;
        end
        if fixationRect_ELHost(3) > screenXpixels
            fixationRect_ELHost(3) = screenXpixels;
        end
        if fixationRect_ELHost(2) > screenYpixels
            fixationRect_ELHost(2) = screenYpixels;
        end
        if fixationRect_ELHost(4) > screenYpixels
            fixationRect_ELHost(4) = screenYpixels;
        end
        fixationRect_ELHost = round(fixationRect_ELHost);
        Eyelink('Command', 'clear_screen 0'); % Clear Host PC display from any previus drawing
        Eyelink('Command', 'draw_filled_box %d %d %d %d 15', fixationRect_ELHost(1), fixationRect_ELHost(2), fixationRect_ELHost(3), fixationRect_ELHost(4)); % Fixation window
        
        tempSquareXpos = round(squareXpos - reCenterGazeBias_x);
        tempSquareYpos = round(squareYpos - reCenterGazeBias_y);
        tempSquareXpos(tempSquareXpos < 0) = 0;
        tempSquareXpos(tempSquareXpos > screenXpixels) = screenXpixels;
        tempSquareYpos(tempSquareYpos < 0) = 0;
        tempSquareYpos(tempSquareYpos > screenYpixels) = screenYpixels;
        
        for tempi=1:numSquares
            Eyelink('Command', 'draw_filled_box %d %d %d %d 15',...
                tempSquareXpos(tempi)-20, tempSquareYpos(tempi)-20,...
                tempSquareXpos(tempi)+20, tempSquareYpos(tempi)+20);
        end
        
        
        Eyelink('Command', 'draw_cross %d %d 0', tempWinPos(1), tempWinPos(2)); % Central crosshairs
        
        ELHost_changeFLag = 0;
    end
    % Check that eye tracker is  still recording. Otherwise close and transfer copy of EDF file to Display PC
    err = Eyelink('CheckRecording');
    if(err ~= 0)
        fprintf('EyeLink Recording stopped!\n');
        exitProgram = 1;
    end
    
    %elEvt_total = cell(1, meanGazeNum);
    elEvt_total = cell(1, 1);
    tempMeanGaze_x = 0;
    tempMeanGaze_y = 0;
    tempCount = 0;
    
    t0_getMeanGaze = tic;
    while toc(t0_getMeanGaze) < (meanGazeNum+1)/1000
        if tempCount == meanGazeNum
            break
        else
            if Eyelink('NewFloatSampleAvailable') > 0
                elEvt_total{1} = Eyelink('NewestFloatSample');
                tempMeanGaze_x = tempMeanGaze_x + elEvt_total{1}.gx(eyeUsed+1);
                tempMeanGaze_y = tempMeanGaze_y + elEvt_total{1}.gy(eyeUsed+1);
                tempCount = tempCount + 1;
            end
        end
    end
    if tempCount > 0
        tempMeanGaze_x = tempMeanGaze_x / tempCount;
        tempMeanGaze_y = tempMeanGaze_y / tempCount;
        elEvt.gx = tempMeanGaze_x;
        elEvt.gy = tempMeanGaze_y;
    end
    %fprintf('toc(t0_getMeanGaze)=%.4f, tempCount=%d!\n', toc(t0_getMeanGaze), tempCount);
        
%     for tempi=1:meanGazeNum
% %         while Eyelink('NewFloatSampleAvailable') < 1
% %         end
% %         elEvt_total{tempi} = Eyelink('NewestFloatSample');
% %         tempMeanGaze_x = tempMeanGaze_x + elEvt_total{tempi}.gx(eyeUsed+1);
% %         tempMeanGaze_y = tempMeanGaze_y + elEvt_total{tempi}.gy(eyeUsed+1);
%         if Eyelink('NewFloatSampleAvailable') > 0         
%             elEvt_total{tempi} = Eyelink('NewestFloatSample');
%             tempMeanGaze_x = tempMeanGaze_x + elEvt_total{tempi}.gx(eyeUsed+1);
%             tempMeanGaze_y = tempMeanGaze_y + elEvt_total{tempi}.gy(eyeUsed+1);
%         else
%             tempCount = tempCount + 1;
%             %fprintf('No NewFloatSample, tempCount=%d!\n', tempCount);
%         end
%     end
%     if meanGazeNum > tempCount
%         tempMeanGaze_x = tempMeanGaze_x / (meanGazeNum-tempCount);
%         tempMeanGaze_y = tempMeanGaze_y / (meanGazeNum-tempCount);
%         elEvt.gx = tempMeanGaze_x;
%         elEvt.gy = tempMeanGaze_y;
%     else
%         tempMeanGaze_x = 0;
%         tempMeanGaze_y = 0;
%     end
% %     elEvt.gx = tempMeanGaze_x;
% %     elEvt.gy = tempMeanGaze_y;

    x_avg = elEvt.gx + reCenterGazeBias_x;
    y_avg = elEvt.gy + reCenterGazeBias_y;
    

    if if_mouse_mannual_eyeResponse == 1
        [mx, my, buttons] = GetMouse(window);
        x_avg = mx;
        y_avg = my;
    end  
    
    fixOn_flag = 1;
    fixOn_flag_delay2 = 1;
    %fixOn_tolerenceFlag = 0;
    if if_freeview0_fixation1 == 1
        fixDis = norm([x_avg y_avg] - [crossCenter(1) crossCenter(2)]);
        if fixDis < fixWinSize_DisplayPC
            %Screen('DrawDots', window, [x_avg y_avg], 20, [255 0 0], [], 2);
            fixOn_flag = 1;
            fixOn_tolerenceFlag = 0;
            if isShowing == 1 || isDelay == 1
                currentTrial_fixPos.x = [currentTrial_fixPos.x x_avg];
                currentTrial_fixPos.y = [currentTrial_fixPos.y y_avg];
            end
        else
            %Screen('DrawDots', window, [x_avg y_avg], 20, [0 255 0], [], 2);
            if if_fixTolerence == 1 && (isShowing == 1 || isDelay == 1)
                if fixOn_tolerenceFlag == 0                
                    fixOn_tolerenceFlag = 1;
                    t0_fixTolerence = tic;    
                    fprintf("FixTolerence start!\n");
                else
                    if toc(t0_fixTolerence) > fixTolerenceTime/1000
                        fixOn_flag = 0;
                    end
                end
            else
                fixOn_flag = 0;
            end
        end
        
            if fixDis < fixWinSize_delay2
                fixOn_flag_delay2 = 1;
                fixOn_tolerenceFlag_delay2 = 0;
            else
                if if_fixTolerence == 1 && isTwoDelay == 1
                    if fixOn_tolerenceFlag_delay2 == 0                
                        fixOn_tolerenceFlag_delay2 = 1;
                        t0_fixTolerence_delay2 = tic;    
                        fprintf("FixTolerence start!\n");
                    else
                        if toc(t0_fixTolerence_delay2) > fixTolerenceTime_delay2/1000
                            fixOn_flag_delay2 = 0;
                        end
                    end
                else
                    fixOn_flag_delay2 = 0;
                end
            end    
    
    else
        %Screen('DrawDots', window, [x_avg y_avg], 20, [255 255 255], [], 2);
    end

    
    
    if keyCode(applicationKey)
        isManual_last = isManual;
        isManual = 1;
        fixOn_flag = 1;
    else
        isManual_last = isManual;
        isManual = 0;
    end
    
%     if keyCode(leftaltKey) == 1 && keyCode(questionKey) == 1
%         isManual_last = isManual;
%         isManual = 1;
%         fixOn_flag = 1;
%     end       

    if keyCode(leftaltKey) == 1 && keyCode(questionKey) == 1 && lastKeyCode(questionKey) == 0
        if_disableRightSelection = 1 - if_disableRightSelection;
        fprintf('if_disableRightSelection=%d.\n',if_disableRightSelection);
    end     
    
    
    if keyCode(leftaltKey) == 1 && keyCode(orKey) == 1 && lastKeyCode(orKey) == 0
        old_data = m_marker.Data;
        m_marker.Data = uint8(123); 
        WaitSecs(0.1);
        giveWater(initGiveWaterTime);
        m_marker.Data = uint8(old_data); 
    end    

    if keyCode(leftaltKey) == 1 && keyCode(minusKey) == 1 && lastKeyCode(minusKey) == 0
        smallWaterInterval = smallWaterInterval - 100;
        if smallWaterInterval < 100
            smallWaterInterval = 100;
        end
        smallWater = totalWater/7000*smallWaterInterval;
        %smallWater = 3 * smallWater^(2/3);
        smallWater = smallWater^(0.75);%0.95-->0.9-->0.85
        fprintf("smallWaterInterval = %d, smallWater = %.2f\n",smallWaterInterval, smallWater);
    end
    if keyCode(leftaltKey) == 1 && keyCode(equalKey) == 1 && lastKeyCode(equalKey) == 0
        smallWaterInterval = smallWaterInterval + 100;
        smallWater = totalWater/7000*smallWaterInterval;
        smallWater = smallWater^(0.75);
        fprintf("smallWaterInterval = %d, smallWater = %.2f\n",smallWaterInterval, smallWater);        
    end    
    
%     if keyCode(leftaltKey) == 1 && keyCode(minusKey) == 1 && lastKeyCode(minusKey) == 0
%         preFixationTime = preFixationTime - 200;
% %         if preFixationTime < 0
% %             preFixationTime = 0;
% %         end
%         fprintf("preFixationTime = %d\n",mean(preFixationTime));
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(equalKey) == 1 && lastKeyCode(equalKey) == 0
%         preFixationTime = preFixationTime + 200;
%         fprintf("preFixationTime = %d\n",mean(preFixationTime));
%     end
    
%     if keyCode(leftaltKey) == 1 && keyCode(minusKey) == 1 && lastKeyCode(minusKey) == 0
%         initGiveWaterTime = initGiveWaterTime - 10;
%         if initGiveWaterTime < 2
%             initGiveWaterTime = 2;
%         end
%         tempWaterTime = initGiveWaterTime + floor(trial_count/100)*waterStepPer100;
%         fprintf("initGiveWaterTime = %d, tempWaterTime = %d\n",initGiveWaterTime, tempWaterTime);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(equalKey) == 1 && lastKeyCode(equalKey) == 0
%         initGiveWaterTime = initGiveWaterTime + 10;
%         tempWaterTime = initGiveWaterTime + floor(trial_count/100)*waterStepPer100;
%         fprintf("initGiveWaterTime = %d, tempWaterTime = %d\n",initGiveWaterTime, tempWaterTime);
%     end
    
%     if keyCode(leftaltKey) == 1 && keyCode(nineKey) == 1 && lastKeyCode(nineKey) == 0
%         fixationTime = fixationTime - 100;
%         if fixationTime(trial_count) < 0 
%             fixationTime = zeros(1, trial_num);
%         end
%         fprintf("fixationTime = %d\n",mean(fixationTime));
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(zeroKey) == 1 && lastKeyCode(zeroKey) == 0
%         fixationTime = fixationTime + 100;
%         fprintf("fixationTime = %d\n",mean(fixationTime));
%     end
%     %InternalMemoryWaterRatio
%     if keyCode(leftaltKey) == 1 && keyCode(nineKey) == 1 && lastKeyCode(nineKey) == 0
%         InternalMemoryWaterRatio = InternalMemoryWaterRatio - 0.05;
%         if InternalMemoryWaterRatio < 0
%             InternalMemoryWaterRatio = 0;
%         end
%         fprintf("InternalMemoryWaterRatio = %.2f\n",InternalMemoryWaterRatio);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(zeroKey) == 1 && lastKeyCode(zeroKey) == 0
%         InternalMemoryWaterRatio = InternalMemoryWaterRatio + 0.05;
%         fprintf("InternalMemoryWaterRatio = %.2f\n",InternalMemoryWaterRatio);
%     end
     

%     if keyCode(leftaltKey) == 1 && keyCode(nineKey) == 1 && lastKeyCode(nineKey) == 0
%         stimulusExtendTime = stimulusExtendTime - 100;
%         if stimulusExtendTime < 0
%             stimulusExtendTime = 0;
%         end
%         fprintf("stimulusExtendTime = %d\n",stimulusExtendTime);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(zeroKey) == 1 && lastKeyCode(zeroKey) == 0
%         stimulusExtendTime = stimulusExtendTime + 100;
%         fprintf("stimulusExtendTime = %d\n",stimulusExtendTime);
%     end
    
    
%     if keyCode(leftaltKey) == 1 && keyCode(oKey) == 1 && lastKeyCode(oKey) == 0
%         InternalMemoryWaterRatio = InternalMemoryWaterRatio - 0.05;
%         if InternalMemoryWaterRatio < 0
%             InternalMemoryWaterRatio = 0;
%         end
%         fprintf("InternalMemoryWaterRatio = %.2f\n",InternalMemoryWaterRatio);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(pKey) == 1 && lastKeyCode(pKey) == 0
%         InternalMemoryWaterRatio = InternalMemoryWaterRatio + 0.05;
%         fprintf("InternalMemoryWaterRatio = %.2f\n",InternalMemoryWaterRatio);
%     end

%     if keyCode(leftaltKey) == 1 && keyCode(oKey) == 1 && lastKeyCode(oKey) == 0
%         dynamicFixThreshold = dynamicFixThreshold - 0.05;
%         if dynamicFixThreshold < 0
%             dynamicFixThreshold = 0;
%         end
%         fprintf("dynamicFixThreshold = %.2f\n",dynamicFixThreshold);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(pKey) == 1 && lastKeyCode(pKey) == 0
%         dynamicFixThreshold = dynamicFixThreshold + 0.05;
%         fprintf("dynamicFixThreshold = %.2f\n",dynamicFixThreshold);
%     end

%     if keyCode(leftaltKey) == 1 && keyCode(oKey) == 1 && lastKeyCode(oKey) == 0
%         submitHoldingTime = submitHoldingTime - 5;
%         if submitHoldingTime < 0
%             submitHoldingTime = 0;
%         end
%         fprintf("submitHoldingTime = %d\n",submitHoldingTime);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(pKey) == 1 && lastKeyCode(pKey) == 0
%         submitHoldingTime = submitHoldingTime + 5;
%         fprintf("submitHoldingTime = %d\n",submitHoldingTime);
%     end
    if keyCode(leftaltKey) == 1 && keyCode(oKey) == 1 && lastKeyCode(oKey) == 0
        eyeResponseTime = eyeResponseTime - 10;
        if eyeResponseTime < 0
            eyeResponseTime = 0;
        end
        fprintf("eyeResponseTime = %d\n",eyeResponseTime);
    end
    if keyCode(leftaltKey) == 1 && keyCode(pKey) == 1 && lastKeyCode(pKey) == 0
        eyeResponseTime = eyeResponseTime + 10;
        fprintf("eyeResponseTime = %d\n",eyeResponseTime);
    end

%     if keyCode(leftaltKey) == 1 && keyCode(oKey) == 1 && lastKeyCode(oKey) == 0
%         fixTolerenceTime = fixTolerenceTime - 10;
%         if fixTolerenceTime < 0
%             fixTolerenceTime = 0;
%         end
%         fprintf("fixTolerenceTime = %d\n",fixTolerenceTime);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(pKey) == 1 && lastKeyCode(pKey) == 0
%         fixTolerenceTime = fixTolerenceTime + 10;
%         fprintf("fixTolerenceTime = %d\n",fixTolerenceTime);
%     end


%     if keyCode(leftaltKey) == 1 && keyCode(oKey) == 1 && lastKeyCode(oKey) == 0
%         pointShowTime = pointShowTime - 10;
%         if pointShowTime < 0
%             pointShowTime = 0;
%         end
%         fprintf("pointShowTime = %d\n",pointShowTime);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(pKey) == 1 && lastKeyCode(pKey) == 0
%         pointShowTime = pointShowTime + 10;
%         fprintf("pointShowTime = %d\n",pointShowTime);
%     end

%     if keyCode(leftaltKey) == 1 && keyCode(oKey) == 1 && lastKeyCode(oKey) == 0
%         endBet_fastSlowBoundary = endBet_fastSlowBoundary - 25;
%         if endBet_fastSlowBoundary < 0
%             endBet_fastSlowBoundary = 0;
%         end
%         fprintf("endBet_fastSlowBoundary = %d\n",endBet_fastSlowBoundary);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(pKey) == 1 && lastKeyCode(pKey) == 0
%         endBet_fastSlowBoundary = endBet_fastSlowBoundary + 25;
%         fprintf("endBet_fastSlowBoundary = %d\n",endBet_fastSlowBoundary);
%     end
    
    if keyCode(leftaltKey) == 1 && keyCode(leftBracketsKey) == 1 && lastKeyCode(leftBracketsKey) == 0
        consecutiveCorrectLimit = consecutiveCorrectLimit - 1;
        if consecutiveCorrectLimit < 0
            consecutiveCorrectLimit = 0;
        end
        fprintf("consecutiveCorrectLimit = %d\n",consecutiveCorrectLimit);
    end
    if keyCode(leftaltKey) == 1 && keyCode(rightBracketsKey) == 1 && lastKeyCode(rightBracketsKey) == 0
        consecutiveCorrectLimit = consecutiveCorrectLimit + 1;
        fprintf("consecutiveCorrectLimit = %d\n",consecutiveCorrectLimit);
    end    
    
    
% %     twoFixationTime
%     if keyCode(leftaltKey) == 1 && keyCode(commaKey) == 1 && lastKeyCode(commaKey) == 0
%         twoFixationTime = twoFixationTime - 50;
%         if twoFixationTime(trial_count) < 0 
%             twoFixationTime = zeros(1, trial_num);
%         end
%         %fprintf("twoFixationTime = %d\n",twoFixationTime(trial_count));
%         fprintf("twoFixationTime = %d\n",mean(twoFixationTime));
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(fullstopKey) == 1 && lastKeyCode(fullstopKey) == 0
%         twoFixationTime = twoFixationTime + 50;
%         %fprintf("twoFixationTime = %d\n",twoFixationTime(trial_count));
%         fprintf("twoFixationTime = %d\n",mean(twoFixationTime));
%     end
    
    changeFlag_cumulativeCorrectLimit_seqSwitch = 0;
%     if keyCode(leftaltKey) == 1 && keyCode(commaKey) == 1 && lastKeyCode(commaKey) == 0
%         cumulativeCorrectLimit_seqLengthSwitch = cumulativeCorrectLimit_seqLengthSwitch - 1;
%         if cumulativeCorrectLimit_seqLengthSwitch < 0
%             cumulativeCorrectLimit_seqLengthSwitch = 0;
%         end
%         changeFlag_cumulativeCorrectLimit_seqSwitch = 1;
%         fprintf("cumulativeCorrectLimit_seqLengthSwitch = %d\n",cumulativeCorrectLimit_seqLengthSwitch);        
%     end    
%     if keyCode(leftaltKey) == 1 && keyCode(fullstopKey) == 1 && lastKeyCode(fullstopKey) == 0
%         cumulativeCorrectLimit_seqLengthSwitch = cumulativeCorrectLimit_seqLengthSwitch + 1;
%         changeFlag_cumulativeCorrectLimit_seqSwitch = 1;
%         fprintf("cumulativeCorrectLimit_seqLengthSwitch = %d\n",cumulativeCorrectLimit_seqLengthSwitch);        
%     end
    if changeFlag_cumulativeCorrectLimit_seqSwitch == 1
        if ifSeq_length_3to1 == 1
            temp_seq_length = [1 4];
            temp_seq_length2 = cell(1, seq_length_rangeTail-seq_length_rangeHead+1);
            for i=1:seq_length_rangeTail-seq_length_rangeHead+1
                temp_seq_length2{i} = repmat(temp_seq_length(i), 1, cumulativeCorrectLimit_seqLengthSwitch);
            end
        elseif ifSeq_length_3to1 == 0
            temp_seq_length = seq_length_rangeHead:seq_length_rangeTail;
            temp_seq_length2 = cell(1, seq_length_rangeTail-seq_length_rangeHead+1);
            for i=1:seq_length_rangeTail-seq_length_rangeHead+1
                temp_seq_length2{i} = repmat(temp_seq_length(i), 1, cumulativeCorrectLimit_seqLengthSwitch);
            end
        end
        
        temp_blockIndex = repmat(temp_seq_length, 1, ceil(trial_num/cumulativeCorrectLimit_seqLengthSwitch/length(temp_seq_length)));
        temp_blockIndex = temp_blockIndex(randperm(length(temp_blockIndex)));
        psudo_seq_length = zeros(1, trial_num);
        for i=1:length(temp_blockIndex)
            for j=1:cumulativeCorrectLimit_seqLengthSwitch
%                 if (i-1)*cumulativeCorrectLimit_seqLengthSwitch+j >= trial_num
%                     break
%                 end
                psudo_seq_length((i-1)*cumulativeCorrectLimit_seqLengthSwitch+j) = temp_seq_length2{find(temp_seq_length==temp_blockIndex(i))}(j);
            end
        end
        
        seq_length(trial_count: trial_num) = psudo_seq_length(trial_count: trial_num);
        isRelaSeqChange = 1;
    end
   
    
    
    if keyCode(leftaltKey) == 1 && keyCode(commaKey) == 1 && lastKeyCode(commaKey) == 0
        numFrames_rangeHead = numFrames_rangeHead - 1;
        numFrames_rangeTail = numFrames_rangeTail - 1;
        if numFrames_rangeHead < seq_length_rangeHead
            numFrames_rangeHead = seq_length_rangeHead;
        end
        if numFrames_rangeTail < seq_length_rangeHead
            numFrames_rangeTail = seq_length_rangeHead;
        end
        
        %         sort_correctCount = zeros(1, numFrames_rangeTail);
        %         sort_allCount = zeros(1, numFrames_rangeTail);
        %         sort_isCorrect = zeros(1, numFrames_rangeTail);
        
        sort_correctCount_last40 = zeros(40, numFrames_rangeTail);
        sort_allCount_last40 = zeros(40, numFrames_rangeTail);
        sort_correctCount_last40 = zeros(40, numFrames_rangeTail);        
        
        isRelaSeqChange = 1;
        
        temp_numFrames = numFrames_rangeHead:numFrames_rangeTail;
        temp_numFrames2 = repmat(temp_numFrames,1,ceil((trial_num)/(length(temp_numFrames))));
        numFrames = temp_numFrames2(randperm((trial_num)));
        
        for temp_trial_count=trial_count:trial_num
            redFillProbWeight{temp_trial_count} = ones(1, numFrames(temp_trial_count)) / numFrames(temp_trial_count);
        end
        
        fprintf("numFrames_rangeHead = %d, numFrames_rangeTail = %d\n",...
            numFrames_rangeHead, numFrames_rangeTail);
    end
    if keyCode(leftaltKey) == 1 && keyCode(fullstopKey) == 1 && lastKeyCode(fullstopKey) == 0
        numFrames_rangeHead = numFrames_rangeHead + 1;
        numFrames_rangeTail = numFrames_rangeTail + 1;
        if numFrames_rangeHead > numSquares
            numFrames_rangeHead = numSquares;
        end
        if numFrames_rangeTail > numSquares
            numFrames_rangeTail = numSquares;
        end
        
        sort_correctCount = [ sort_correctCount 0]; %#ok<AGROW>
        sort_allCount = [ sort_allCount 0]; %#ok<AGROW>
        sort_isCorrect = [ sort_isCorrect nan]; %#ok<AGROW>
        
        sort_correctCount_last40 = zeros(40, numFrames_rangeTail);
        sort_allCount_last40 = zeros(40, numFrames_rangeTail);
        sort_correctCount_last40 = zeros(40, numFrames_rangeTail);        
        
        isRelaSeqChange = 1;
        
        temp_numFrames = numFrames_rangeHead:numFrames_rangeTail;
        temp_numFrames2 = repmat(temp_numFrames,1,ceil((trial_num)/(length(temp_numFrames))));
        numFrames = temp_numFrames2(randperm((trial_num)));
        
        for temp_trial_count=trial_count:trial_num
            redFillProbWeight{temp_trial_count} = ones(1, numFrames(temp_trial_count)) / numFrames(temp_trial_count);
        end
        
        fprintf("numFrames_rangeTail = %d, numFrames_rangeTail = %d\n",...
            numFrames_rangeHead, numFrames_rangeTail);
    end
    
%     %freeTouch_RT_limit
%     if keyCode(leftaltKey) == 1 && keyCode(commaKey) == 1 && lastKeyCode(commaKey) == 0
%         freeTouch_RT_limit = freeTouch_RT_limit - 100;
%         if freeTouch_RT_limit < 0
%             freeTouch_RT_limit = 0;
%         end
%         fprintf("freeTouch_RT_limit = %d\n",freeTouch_RT_limit);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(fullstopKey) == 1 && lastKeyCode(fullstopKey) == 0
%         freeTouch_RT_limit = freeTouch_RT_limit + 100;
%         fprintf("freeTouch_RT_limit = %d\n",freeTouch_RT_limit);
%     end

%     %probFreeTouch
%     if keyCode(leftaltKey) == 1 && keyCode(commaKey) == 1 && lastKeyCode(commaKey) == 0
%         probFreeTouch = probFreeTouch - 0.05;
%         if probFreeTouch < 0
%             probFreeTouch = 0;
%         end
%         temp = ifFreeTouch_origin * ones(1, trial_num) .* (randArray4<probFreeTouch);
%         ifFreeTouch(trial_count:end) = temp(trial_count:end);        
%         fprintf("probFreeTouch = %.2f\n",probFreeTouch);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(fullstopKey) == 1 && lastKeyCode(fullstopKey) == 0
%         probFreeTouch = probFreeTouch + 0.05;
%         temp = ifFreeTouch_origin * ones(1, trial_num) .* (randArray4<probFreeTouch);
%         ifFreeTouch(trial_count:end) = temp(trial_count:end);        
%         fprintf("probFreeTouch = %.2f\n",probFreeTouch);
%     end
    

%     if keyCode(leftaltKey) == 1 && keyCode(zKey) == 1 && lastKeyCode(zKey) == 0
%         choiceCondition_boardLine1 = choiceCondition_boardLine1 - 0.03;
%         if choiceCondition_boardLine1 < 0 
%             choiceCondition_boardLine1 = 0;
%         end
%         fprintf("choiceCondition_boardLine1 = %.2f\n",choiceCondition_boardLine1);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(xKey) == 1 && lastKeyCode(xKey) == 0
%         choiceCondition_boardLine1 = choiceCondition_boardLine1 + 0.03;
%         if choiceCondition_boardLine1 > choiceCondition_boardLine2
%             choiceCondition_boardLine1 = choiceCondition_boardLine2;
%         end
%         fprintf("choiceCondition_boardLine1 = %.2f\n",choiceCondition_boardLine1);
%     end
% 
%     if keyCode(leftaltKey) == 1 && keyCode(zKey) == 1 && lastKeyCode(zKey) == 0
%         cumulativeCorrectLimit_greenRedSwitch = cumulativeCorrectLimit_greenRedSwitch - 1;
%         if cumulativeCorrectLimit_greenRedSwitch < 0
%             cumulativeCorrectLimit_greenRedSwitch = 0;
%         end
%         fprintf("cumulativeCorrectLimit_greenRedSwitch = %d\n",cumulativeCorrectLimit_greenRedSwitch);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(xKey) == 1 && lastKeyCode(xKey) == 0
%         cumulativeCorrectLimit_greenRedSwitch = cumulativeCorrectLimit_greenRedSwitch + 1;
%         fprintf("cumulativeCorrectLimit_greenRedSwitch = %d\n",cumulativeCorrectLimit_greenRedSwitch);
%     end
    
    %     if keyCode(leftaltKey) == 1 && keyCode(zKey) == 1 && lastKeyCode(zKey) == 0
    %         propExplore = propExplore - 0.1;
    %         if propExplore < 0
    %             propExplore = 0;
    %         end
    %         fprintf("propExplore = %.1f\n",propExplore);
    %     end
    %     if keyCode(leftaltKey) == 1 && keyCode(xKey) == 1 && lastKeyCode(xKey) == 0
    %         propExplore = propExplore + 0.1;
    %         if propExplore > 1
    %             propExplore = 1;
    %         end
    %         fprintf("propExplore = %.1f\n",propExplore);
    %     end
    
%     if keyCode(leftaltKey) == 1 && keyCode(qKey) == 1 && lastKeyCode(qKey) == 0
%         isProbWeightChange = 1;
%         isProbWeightEven = 0;
%         isManual_rFPW = 1;        
%         addTarget_probWeight = 1;
%         fprintf("addTarget_probWeight = %d\n",addTarget_probWeight);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(qKey) == 1 && lastKeyCode(qKey) == 0
%         isProbWeightChange = 1;
%         isProbWeightEven = 0;
%         isManual_rFPW = 1;        
%         addTarget_probWeight = numFrames_rangeTail;
%         fprintf("addTarget_probWeight = %d\n",addTarget_probWeight);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(wKey) == 1 && lastKeyCode(wKey) == 0
%         isProbWeightChange = 1;
%         isProbWeightEven = 0;
%         isManual_rFPW = 1;
%         addTarget_probWeight = 2;
%         fprintf("addTarget_probWeight = %d\n",addTarget_probWeight);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(eKey) == 1 && lastKeyCode(eKey) == 0
%         isProbWeightChange = 1;
%         isProbWeightEven = 0;
%         isManual_rFPW = 1;
%         addTarget_probWeight = 3;
%         fprintf("addTarget_probWeight = %d\n",addTarget_probWeight);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(rKey) == 1 && lastKeyCode(rKey) == 0
%         isProbWeightChange = 1;
%         isProbWeightEven = 0;
%         isManual_rFPW = 1;
%         addTarget_probWeight = 4;
%         fprintf("addTarget_probWeight = %d\n",addTarget_probWeight);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(tKey) == 1 && lastKeyCode(tKey) == 0
%         isProbWeightChange = 1;
%         isProbWeightEven = 0;
%         isManual_rFPW = 1;
%         addTarget_probWeight = 5;
%         fprintf("addTarget_probWeight = %d\n",addTarget_probWeight);
%     end
    
    if isProbWeightEven == 1
        for temp_trial_count=trial_count:trial_num
            for i=1:numFrames(temp_trial_count)
                redFillProbWeight{temp_trial_count}(i) = 1/numFrames(temp_trial_count);
            end
        end
        
        isProbWeightEven = 0;
        isProbWeightChange = 0;
        isRelaSeqChange = 1;
    end
    
    
    if isProbWeightChange == 1
        %update redFillProbWeight
        if isManual_rFPW == 1
            for temp_trial_count=trial_count:trial_num
                for i=addTarget_probWeight:numFrames(temp_trial_count)
                    if i == addTarget_probWeight
                        redFillProbWeight{temp_trial_count}(i) = redFillProbWeight{temp_trial_count}(i) + 0.1;
                        if redFillProbWeight{temp_trial_count}(i) > 1
                            redFillProbWeight{temp_trial_count} = zeros(1, numFrames(temp_trial_count));
                            redFillProbWeight{temp_trial_count}(i) = 1;
                        end
                    else
                        temp = 0.1 * redFillProbWeight{temp_trial_count}(i) / (1 - redFillProbWeight{temp_trial_count}(addTarget_probWeight));
                        if redFillProbWeight{temp_trial_count}(i) > temp
                            redFillProbWeight{temp_trial_count}(i) = redFillProbWeight{temp_trial_count}(i) - temp;
                        end
                        redFillProbWeight{temp_trial_count} = redFillProbWeight{temp_trial_count} / ...
                            sum(redFillProbWeight{temp_trial_count});
                    end
                end
                for i=1:addTarget_probWeight
                    if i < addTarget_probWeight
                        temp = 0.1 * redFillProbWeight{temp_trial_count}(i) / (1 - redFillProbWeight{temp_trial_count}(addTarget_probWeight));
                        if redFillProbWeight{temp_trial_count}(i) > temp
                            redFillProbWeight{temp_trial_count}(i) = redFillProbWeight{temp_trial_count}(i) - temp;
                        end
                        redFillProbWeight{temp_trial_count} = redFillProbWeight{temp_trial_count} / ...
                            sum(redFillProbWeight{temp_trial_count});
                    end
                end
            end
        elseif isManual_rFPW == 0
            %temp_rFPW = 1 - sort_isCorrect(1:numFrames(trial_count));
            if trial_count <= 40
                temp_rFPW = (1 - sort_isCorrect(1:numFrames(trial_count))) .^2.5;
            else
                temp_rFPW = (1 - sort_accuracy_last40(1:numFrames(trial_count))) .^2.5;
            end
            temp_rFPW(temp_rFPW<0.07) = 0.07;
            
            
            for i=1:length(temp_rFPW)
                if isnan(temp_rFPW(i))
                    temp_rFPW(i) = 1;
                end
            end
            for temp_trial_count=trial_count+1:trial_num
                if numFrames(temp_trial_count) >= numFrames(trial_count)
                    temp_rFPW2 = temp_rFPW./sum(temp_rFPW);
                    temp_rFPW3 = (temp_rFPW2 + redFillProbWeight{trial_count})/2;
                    redFillProbWeight{temp_trial_count} = temp_rFPW3;
                    %elseif numFrames(temp_trial_count) > numFrames(trial_count)
                elseif numFrames(temp_trial_count) < numFrames(trial_count)
                    temp_rFPW41 = temp_rFPW(1:numFrames(temp_trial_count));
                    temp_rFPW42 = temp_rFPW41./sum(temp_rFPW41);
                    redFillProbWeight{temp_trial_count} = temp_rFPW42;
                end
            end
        end
        
        if isManual_rFPW==1
            fprintf("Current redFillProbWeight = ");
            for i=1:numFrames(trial_count+1)
                fprintf('%.2f, ', redFillProbWeight{trial_count+1}(i));
            end
            fprintf("\n");
        end
        
        isProbWeightChange = 0;
        isRelaSeqChange = 1;
    end
    
    if max(seq_length) <= 2
        if isRelaSeqChange == 1
            %update relativeSequence
            for temp_trial_count=trial_count:trial_num
                for seq_index=1:seq_length(temp_trial_count)
                    for i=1:numFrames(temp_trial_count)
                        if i==1
                            if temprand(seq_index, temp_trial_count) < redFillProbWeight{temp_trial_count}(1)
                                relativeSequence{temp_trial_count}(seq_index) = 1;
                            end
                        else
                            if temprand(seq_index, temp_trial_count) > sum(redFillProbWeight{temp_trial_count}(1:i-1))...
                                    && temprand(seq_index, temp_trial_count) < sum(redFillProbWeight{temp_trial_count}(1:i))
                                relativeSequence{temp_trial_count}(seq_index) = i;
                            end
                        end
                    end
                    %Change same red fill so that they are unique
                    if length(relativeSequence{temp_trial_count})-length(unique(relativeSequence{temp_trial_count})) > 0
                        tempValue = relativeSequence{temp_trial_count}(end);
                        seg_temprand = temprand(seq_index, temp_trial_count);
                        if tempValue > 1
                            seg_temprand = seg_temprand - sum(redFillProbWeight{temp_trial_count}(1:tempValue-1));
                        end
                        tempValueWeight = redFillProbWeight{temp_trial_count}(tempValue);
                        for i=1:numFrames(temp_trial_count)
                            if i < tempValue
                                if i == 1
                                    if seg_temprand/tempValueWeight < ...
                                            redFillProbWeight{temp_trial_count}(1)/(1-tempValueWeight)
                                        relativeSequence{temp_trial_count}(end) = 1;
                                    end
                                else
                                    if seg_temprand/tempValueWeight > ...
                                            sum(redFillProbWeight{temp_trial_count}(1:i-1))/(1-tempValueWeight)...
                                            && seg_temprand/tempValueWeight <...
                                            sum(redFillProbWeight{temp_trial_count}(1:i))/(1-tempValueWeight)
                                        relativeSequence{temp_trial_count}(end) = i;
                                    end
                                end
                            elseif i > tempValue
                                if seg_temprand/tempValueWeight > ...
                                        (sum(redFillProbWeight{temp_trial_count}(1:i-1))-tempValueWeight)/(1-tempValueWeight)...
                                        && seg_temprand/tempValueWeight <...
                                        (sum(redFillProbWeight{temp_trial_count}(1:i))-tempValueWeight)/(1-tempValueWeight)
                                    relativeSequence{temp_trial_count}(end) = i;
                                end
                                
                            end
                        end
                    end
                end
            end
            isRelaSeqChange = 0;
        end
    elseif max(seq_length) >= 3
        %add this if condition to fix up a break bug
        if isRelaSeqChange == 1
            for temp_trial_count=1:trial_num
                relativeSequence{temp_trial_count} = randperm(numFrames(temp_trial_count), seq_length(temp_trial_count));
            end
            isRelaSeqChange = 0;
        end
    end

%     if keyCode(leftaltKey) == 1 && keyCode(cKey) == 1 && lastKeyCode(cKey) == 0
%         choiceCondition_boardLine2 = choiceCondition_boardLine2 - 0.03;
%         if choiceCondition_boardLine2 < 0 
%             choiceCondition_boardLine2 = 0;
%         end
%         fprintf("choiceCondition_boardLine2 = %.2f\n",choiceCondition_boardLine2);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(vKey) == 1 && lastKeyCode(vKey) == 0
%         choiceCondition_boardLine2 = choiceCondition_boardLine2 + 0.03;
%         if choiceCondition_boardLine2 > 1
%             choiceCondition_boardLine2 = 1;
%         end
%         fprintf("choiceCondition_boardLine2 = %.2f\n",choiceCondition_boardLine2);
%     end    

    if keyCode(leftaltKey) == 1 && keyCode(xKey) == 1 && lastKeyCode(xKey) == 0
        isManual_rFPW = 1;
        isProbWeightChange = 1;
        isProbWeightEven = 1;
        fprintf("Now adjust rFPW evenly!\n");
    end
    if keyCode(leftaltKey) == 1 && keyCode(cKey) == 1 && lastKeyCode(cKey) == 0
        isManual_rFPW = 0;
        isProbWeightChange = 1;
        isProbWeightEven = 0;
        fprintf("Now adjust rFPW automatically!\n");
    end
    
    if keyCode(leftaltKey) == 1 && keyCode(vKey) == 1 && lastKeyCode(vKey) == 0
        fprintf("isManual_rFPW = %d, Current redFillProbWeight = ", isManual_rFPW);
        for i=1:numFrames(trial_count+1)
            fprintf('%.2f, ', redFillProbWeight{trial_count+1}(i));
        end
        fprintf("\n");
    end
    
%     %sKey = KbName('s');%switch between force-choice-g, force-choice-r and free-choice-g/r
%     if keyCode(leftaltKey) == 1 && keyCode(sKey) == 1 && lastKeyCode(sKey) == 0
%         if choiceMode == 0
%             choiceMode = 1;
%             choiceMode_backup = 1;
%             fprintf("Current condition is force-choice-r! ");
%         elseif choiceMode == 1
%             choiceMode = 2;
%             choiceMode_backup = 2;
%             fprintf("Current condition is free-choice-g/r! ");
%         elseif choiceMode == 2
%             choiceMode = 0;
%             choiceMode_backup = 0;
%             fprintf("Current condition is force-choice-g! ");
%         end
%         fprintf("Now choiceMode = %d!\n", choiceMode);
%     end
%     if keyCode(leftaltKey) == 1 && keyCode(sKey) == 1 && lastKeyCode(sKey) == 0
%         if frameMode == 0
%             frameMode = 1;
%         elseif frameMode == 1
%             frameMode = 2;
%         elseif frameMode == 2
%             frameMode = 3;
%         elseif frameMode == 3
%             frameMode = 0;
%         end
%         fprintf("Now frameMode = %d!\n", frameMode);
%     end
    
%     if keyCode(leftaltKey) == 1 && keyCode(dKey) == 1 && lastKeyCode(dKey) == 0
%         if taskMode_noChoice_Green0_Red1 == 0
%             taskMode_noChoice_Green0_Red1 = 1;
%             fprintf("Now taskMode = noChoice_Red1!\n");
%         elseif taskMode_noChoice_Green0_Red1 == 1
%             taskMode_noChoice_Green0_Red1 = 0;
%             fprintf("Now taskMode = noChoice_Green0!\n");
%         end
%     end
    
    
    if keyCode(leftaltKey) == 1 && keyCode(lKey) == 1 && lastKeyCode(lKey) == 0
        if ifLinkTriangleSquare == 0
            ifLinkTriangleSquare = 1;
            fprintf("Now ifLinkTriangleSquare in two trials!\n");
            %         elseif ifLinkTriangleSquare == 1
            %             ifLinkTriangleSquare = 2;
            %             fprintf("Now ifLinkTriangleSquare in one trial!\n");
        elseif ifLinkTriangleSquare == 1
            ifLinkTriangleSquare = 0;
            fprintf("Now ifLinkTriangleSquare = 0!\n");
        end
    end
    
    if keyCode(leftaltKey) == 1 && keyCode(mKey) == 1 && lastKeyCode(mKey) == 0
        if if_mouse_mannual_eyeResponse == 0
            if_mouse_mannual_eyeResponse = 1;            
        elseif if_mouse_mannual_eyeResponse == 1
            if_mouse_mannual_eyeResponse = 0;
        end
        fprintf("Now if_mouse_mannual_eyeResponse = %d!\n",if_mouse_mannual_eyeResponse);
    end    
    

    if keyCode(escapeKey)
        exitProgram = 1;
    end
    
    lastKeyCode = keyCode;
    if trial_count == trial_num
        exitProgram = 1;
    end    
    
    %tic
    %joyState = getJoyState();%time-consuming, about 100ms!!!
    %toc
    if if_eyeResponse == 0
        if joyMode_JJBJoy0_JoyMex1 == 0
            joyState = m.Data;%very quick, using memmapfile
        elseif joyMode_JJBJoy0_JoyMex1 == 1
            tempJS = getAxisState(js);
            if tempJS == 2
                joyState = 1;
            elseif tempJS == 1
                joyState = 0;
            end
        end
    elseif if_eyeResponse == 1
        joyState = 0;
    end
    
    if isManual == 1
        joyState = isManual;
    end
    if isManual_last == 1 && isManual == 0
        joyState = 0;
    end
    newState = joyState;  
    
    % if release the lever, joyState == 0
    % if hold the lever, joyState == 1    
    if joyState == 0
        isHold = 0;
    elseif joyState == 1
        isHold = 1;       
    end
    
    if newState ~=  oldState                         
        %from release to hold
        if newState > oldState
            isReleaseToHold = 1;
        %from hold to release                                                                 
        else 
            isHoldToRelease = 1;
        end
        oldState = newState;
    else
        isReleaseToHold = 0;
        isHoldToRelease = 0;
    end  
         
    
    %% task state machine
    %isShowingInitial-->isShowing-->isDelay-->isReleaseForTouch-->isSelecting...
    %if g:-->isTwoHolding-->isTwoDelay-->isTouching-->isItiInital-->isIti
    %if r:-->isTouching(just offload)-->isTwoHolding-->isTwoDelay-->isTouching-->isItiInital-->isIti
    %---------------------------------------------------------------------------------------------        
    if isShowingInitial == 1
        %% isShowingInitial
        [mx, my, buttons] = GetMouse(window);
        rng('default');%rng('shuffle');
        
        
        %fixationBase_delay2 = [-fixWinSize_delay2 -fixWinSize_delay2 ...
        %   fixWinSize_delay2 fixWinSize_delay2];
        %fixationRect_delay2 = CenterRectOnPointd(fixationBase_delay2, crossCenter(1), crossCenter(2));
        %Screen('FrameRect', window, [1 1 1], fixationRect_delay2);        
        
        if if_mouse_mannual_eyeResponse == 0
            %Screen('FillRect', window, [0.3 0.3 0.3], centeredVerticalCrossRect);
            %Screen('FillRect', window, [0.3 0.3 0.3], centeredHorizontalCrossRect);
            Screen('FillRect', window, [0.45 0.45 0.45], centeredVerticalCrossRect);
            Screen('FillRect', window, [0.45 0.45 0.45], centeredHorizontalCrossRect);
        end
        
        if initialFlag == 0
            holdEvtMappedX = 0;
            holdEvtMappedY = 0;
            holdInside = 0;
        end
        
        if initialFlag == 0
            initialFlag = 1;                                                                          
            
            selectFlag_newSeq(trial_count) = 0;   

            flag_if_breakWithNewSeq = 0;
            if if_breakWithNewSeq == 1 && isWaitGocueBreak(trial_count) == 1
                flag_if_breakWithNewSeq = 1;
                                                    
                %if ifSpecializedTraining_G == 1 && (consecutiveCorrect+cumulativeError) > 0
                if if_cancelBreakForCCL == 1 && (consecutiveCorrect+cumulativeError) > 0
                    flag_if_breakWithNewSeq = 0;
                    fprintf('Cancel break with new seq because of training!\n');
                end
            end
            if trial_count == 1 || consecutiveCorrect >= temp_consecutiveCorrectLimit...
                    || cumulativeError >= cumulativeErrorLimit(trial_count) ...
                    || ( cumulativeError >= encourageCumuErrorLimit && consecutiveCorrect == 1 ) ...
                    || flag_if_breakWithNewSeq == 1
                    %|| ( if_breakWithNewSeq == 1 && isWaitGocueBreak(trial_count) == 1 )
                
                        
                consecutiveCorrect = 0;
                cumulativeError = 0;
                
                generateFramesOK = 0;
                while generateFramesOK == 0
                    frames = randperm(numSquares,numFrames(trial_count));
                    if frameMode == 0
                        tempFlag = min(frames)<=1;
                    elseif frameMode == 1
                        tempFlag = max(frames)>=8;
                    elseif frameMode == 2
                        tempFlag = min(frames)<=1 || max(frames)>=8;
                    elseif frameMode == 3
                        tempFlag = 1;
                    end
                    if max(frames)-min(frames)+1 == numFrames(trial_count) && tempFlag == 1
                        generateFramesOK = 1;      
                        currentFrames = sort(frames);

                        if if_breakWithNewSeq == 1 && isWaitGocueBreak(trial_count) == 1
                            % Get new seq with NEW seq_length
                            if ifFormalTask == 1 
                                randArray5(break_count(trial_count));
                                
                                if if_seqList == 1
                                    
                                    if_newLength = 1;
                                    
                                    if if_newLength == 0
                                        % Get new seq with old seq_length, 20220719
                                        temp_trialNum_oneRun;
                                        temp_trialCount_oneRun = mod(trial_count,temp_trialNum_oneRun);
                                        if temp_trialCount_oneRun == 0
                                            temp_trialCount_oneRun = temp_trialNum_oneRun;
                                        end
                                        
                                        temp_num = temp_trialNum_oneRun;
                                        if (trial_num-trial_count+1) < temp_trialNum_oneRun
                                            temp_num = mod(trial_num,temp_trialNum_oneRun);
                                        end
                                        
                                        temp_range = trial_count+(0:(temp_num-temp_trialCount_oneRun));
                                        temp_range_length = seq_length(temp_range);
                                        temp_range_choice = choiceCondition_flag_list(temp_range);
                                        
                                        temp_range_valid = temp_range((temp_range_length == seq_length(trial_count)) & ...
                                            (temp_range_choice == choiceCondition_flag_list(trial_count)));
                                        
                                        temp_range_valid;
                                        temp_swapIndex = temp_range_valid(...
                                            ceil( randArray5(break_count(trial_count))*length(temp_range_valid) )...
                                            );
                                    end
                                    
                                    if if_newLength == 1
                                        %temp_swapIndex = trial_count + ceil(randArray5(break_count(trial_count)) *...
                                        %    (temp_trialNum_oneRun-temp_trialCount_oneRun));
                                        
                                        temp_trialNum_oneRun;
                                        if if_length1_alone == 0
                                            temp_trialCount_oneRun = mod(trial_count,temp_trialNum_oneRun);
                                            if temp_trialCount_oneRun == 0
                                                temp_trialCount_oneRun = temp_trialNum_oneRun;
                                            end
                                            
                                            temp_num = temp_trialNum_oneRun;
                                            if (trial_num-trial_count+1) < temp_trialNum_oneRun
                                                temp_num = mod(trial_num,temp_trialNum_oneRun);
                                            end
                                            temp_swapIndex = trial_count + ceil(randArray5(break_count(trial_count)) *...
                                                (temp_num-temp_trialCount_oneRun));
                                            
                                        elseif if_length1_alone == 1
                                            trialNum_length1Run;
                                            if trial_count <= trialNum_length1Run
                                                temp_trialCount_oneRun = trial_count;
                                                temp_num = trialNum_length1Run;
                                            elseif trial_count > trialNum_length1Run
                                                temp_trial_count = trial_count - trialNum_length1Run;
                                                
                                                temp_trialCount_oneRun = mod(temp_trial_count,temp_trialNum_oneRun);
                                                if temp_trialCount_oneRun == 0
                                                    temp_trialCount_oneRun = temp_trialNum_oneRun;
                                                end
                                                
                                                temp_num = temp_trialNum_oneRun;
                                                if (trial_num-trial_count+1) < temp_trialNum_oneRun
                                                    temp_num = mod(trial_num-trialNum_length1Run,temp_trialNum_oneRun);
                                                end
                                            end
                                            
                                            temp_swapIndex = trial_count + ceil(randArray5(break_count(trial_count)) *...
                                                (temp_num-temp_trialCount_oneRun));                                            
                                            
                                            
                                        end
                                                                                
                                    end
                                    
                                    %choiceCondition_flag_list
                                    temp_swapFlag = choiceCondition_flag_list(trial_count);
                                    choiceCondition_flag_list(trial_count) = choiceCondition_flag_list(temp_swapIndex);
                                    choiceCondition_flag_list(temp_swapIndex) = temp_swapFlag;
                                    
                                elseif if_seqList == 0
                                    temp_swapIndex = trial_count + ceil(randArray5(break_count(trial_count)) * (length(relativeSequence)-trial_count));
                                end
                                
                            % Get new seq with OLD seq_length    
                            elseif ifFormalTask == 0 
                                randArray5(break_count(trial_count));
                                
                                temp_same_seq_length = find(seq_length == seq_length(trial_count));
                                temp_same_seq_length(temp_same_seq_length <= trial_count) = [];
                                temp_same_seq_length(temp_same_seq_length > trial_num) = [];
                                
                                if length(temp_same_seq_length) >= 20
                                    temp_next_20 = temp_same_seq_length(1:20);
                                else
                                    temp_next_20 = temp_same_seq_length;
                                end
                                
                                temp_swapIndex = temp_next_20(ceil(randArray5(break_count(trial_count)) * length(temp_next_20)));
                            end
                            
                            temp_swapIndex;
                            
                            temp_swapSeq = relativeSequence{trial_count};
                            relativeSequence{trial_count} = relativeSequence{temp_swapIndex};
                            relativeSequence{temp_swapIndex} = temp_swapSeq;
                            
                            seq_length(temp_swapIndex) = seq_length(trial_count);
                            seq_length(trial_count) = length(relativeSequence{trial_count});                           
                            
                            fprintf("Swap with No.%d trial.\n", temp_swapIndex);
                        end
                        
%                         if seq_length(trial_count) == 2
%                             temp_rand;
%                             relativeSequence{trial_count} = [temp_rand(trial_count) temp_rand(trial_count)+1];
%                         end
                        
                        currentSequence = relativeSequence{trial_count} + currentFrames(1) - 1;                        
                        currentSequence = sort(currentSequence);%so that seqs are clockwise
                                                
                        if if_length4_freetouchCCL0 == 1 && length(currentSequence) == 4
                            temp_consecutiveCorrectLimit = 0;
                            ifFreeTouch(trial_count) = 1;
                        elseif if_length1_moreCCL == 1 && length(currentSequence) == 1
                            %temp_consecutiveCorrectLimit = consecutiveCorrectLimit + 1;
                            %temp_consecutiveCorrectLimit = 2;
                            temp_consecutiveCorrectLimit = 1;
                        elseif if_length3_lessCCL == 1 && length(currentSequence) == 3
                            temp_consecutiveCorrectLimit = consecutiveCorrectLimit - 1;
                            if temp_consecutiveCorrectLimit < 0
                                temp_consecutiveCorrectLimit = 0;
                            end
                        else
                            temp_consecutiveCorrectLimit = consecutiveCorrectLimit;
                        end                        
                        
                        selectFlag_newSeq(trial_count) = 1;
                                                
                    end
                end
            else 
                %add this if condition to fix up a break bug
                if isWaitGocueBreak(trial_count) == -1
                    relativeSequence{trial_count} = relativeSequence{trial_count-1};
                    seq_length(trial_count) = seq_length(trial_count-1);
                    ifFreeTouch(trial_count) = ifFreeTouch(trial_count-1);
                end
                
                selectFlag_newSeq(trial_count) = 0;
            end  
            %currentSequence %#ok<NOPTS>   
            %home % To scrolling down in command window
            fprintf("=================currentSequence = ");
            %fprintf("currentSequence = ");
            for tempi=1:length(currentSequence)
                fprintf("%d ", currentSequence(tempi));
            end
            fprintf("\n");
            
            fprintf("Current ifFreeTouch = %d\n", ifFreeTouch(trial_count));
            
            
            if ifGreenRedAutoSwitch == 1
                if cumulativeCorrect_greenRedSwitch >= cumulativeCorrectLimit_greenRedSwitch
                    cumulativeCorrect_greenRedSwitch = 0;
                    if randArray2(trial_count) <= 0.5
                        taskMode_noChoice_Green0_Red1 = 0;
                    elseif randArray2(trial_count) > 0.5
                        taskMode_noChoice_Green0_Red1 = 1;
                    end
                    fprintf("Current randArray2 = %.2f, ", randArray2(trial_count)); 
                    fprintf("taskMode_noChoice_Green0_Red1 = %d\n", taskMode_noChoice_Green0_Red1);
                end
            end
            
            if ifLinkTriangleSquare == 1
                if consecutiveCorrect == 0
                    taskMode_noChoice_Green0_Red1 = 1;
                    tempPointShowTime = 500;
                    tempStimulusExtendTime = stimulusExtendTime;
                elseif consecutiveCorrect == 1
                    taskMode_noChoice_Green0_Red1 = 0;
                    tempPointShowTime = pointShowTime;
                    tempStimulusExtendTime = 10000;                    
                end
            elseif ifLinkTriangleSquare == 0
                tempPointShowTime = pointShowTime;
                tempStimulusExtendTime = stimulusExtendTime;
            end                            
            
            isFrame = zeros(1, numSquares);
            isFill = zeros(1, numSquares);
            isRedFill = zeros(1, numSquares);
            isFilled = zeros(1, numSquares);
            frameColorMatrices = repmat(frameColor',1,numSquares);
            frameColorMatrices2 = repmat(frameColor2',1,numSquares);            
            initialColorMatrices = repmat(normalColor',1,numSquares);
            currentColorMatrices = initialColorMatrices;            
            for i = 1:length(currentFrames)
                if taskMode_noChoice_Green0_Red1 == 0
                    if ismember(currentFrames(i), currentSequence) == 1
                        isRedFill(currentFrames(i)) = 1;
                    else
                        isFill(currentFrames(i)) = 1;
                    end
                elseif taskMode_noChoice_Green0_Red1 == 1
                    if ismember(currentFrames(i), currentSequence) == 1
                        isRedFill(currentFrames(i)) = 1;
                        isFill(currentFrames(i)) = 1;
                    end
                end
                isFrame(currentFrames(i)) = 1;
            end            
            
        end       
        


        if if_fix_noHold0_hold1 == 0                 
            isHold = 1 - isHold;
            if if_eyeResponse == 0
                if TouchEventAvail(dev) > 0
                    isHold = 0;
                end
            end
%             [mx, my, buttons] = GetMouse(window);
%             holdInside = max(IsInRect(mx, my, holdRects_dection'), holdInside); %#ok<*NASGU>
%             while TouchEventAvail(dev) > 0
%                 evt = TouchEventGet(dev, window);
% 
%                 % for touchhold
%                 if evt.Type == 2 || evt.Type == 3 || evt.Type == 4
%                     holdEvtMappedX = evt.MappedX;
%                     holdEvtMappedY = evt.MappedY;
%                 end
%                 holdInside = max( IsInRect(holdEvtMappedX, holdEvtMappedY, holdRects_dection'), holdInside );
% 
%                 if evt.Type == 4 
%                     holdEvtMappedX = 0;
%                     holdEvtMappedY = 0;
%                 end
%             end       
%             holdInside = max( IsInRect(holdEvtMappedX, holdEvtMappedY, holdRects_dection'), holdInside );
%             isHold = holdInside;
%             Screen('FrameRect', window, frameColor, holdRects, 3);
        end
        
        if if_freeview0_fixation1 == 1
            currentColorMatrices = repmat(specialColor',1,numSquares);
            %Screen('FillOval', window, currentColorMatrices(:, currentSequence), showRects(:, currentSequence));
            %Screen('FrameOval', window, frameColorMatrices(:, 1:2), showRects(:, 1:2), frameWidth, frameWidth);

            if fixOn_flag == 1
                %toc(t0_fixOn)
                
                Screen('FillRect', window, [1 1 1]*0.75, centeredVerticalCrossRect);%0.6-->0.75
                Screen('FillRect', window, [1 1 1]*0.75, centeredHorizontalCrossRect);
                Screen('FrameOval', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);

                if toc(t0_fixOn) > preFixationTime(trial_count)/1000
                    %Screen('FillRect', window, [1 1 1]*0.75, centeredVerticalCrossRect);%0.6-->0.75
                    %Screen('FillRect', window, [1 1 1]*0.75, centeredHorizontalCrossRect);
                    preFixationEndFlag = 1;
                    %Screen('FillOval', window, currentColorMatrices(:, currentSequence), showRects(:, currentSequence));

%                     if preFixationTime(trial_count) > 1
%                         initGiveWaterTime = 140;
%                         giveWater(initGiveWaterTime);
%                         isShowingInitial = 0;
%                         isIti = 1;                              
%                         tempIti = ItiTime(trial_count);
%                     end
                    
                end
                %if ifDrawDots_eye == 1
                if ifDrawDots_eye == 1 && isHold == 1
                    %Screen('DrawDots', window, [x_avg y_avg], 20, [255 0 0], [], 2);
                    Screen('DrawDots', window, [x_avg y_avg], 20, [0.5 0.5 0.5], [], 2);
                end
            elseif fixOn_flag == 0
                %if ifDrawDots_eye == 1
                if ifDrawDots_eye == 1 && isHold == 1
                    %Screen('DrawDots', window, [x_avg y_avg], 20, [0 255 0], [], 2);
                end
                %isReleaseToHold = 0;
                t0_fixOn = tic;
            end
            
            %Screen('FrameOval', window, frameColorMatrices(:, 1:numFrames_rangeTail), showRects(:, 1:numFrames_rangeTail), frameWidth, frameWidth);
            
            if last_preFixationEndFlag == 0 && preFixationEndFlag == 1
                fprintf("current preFixationTime = %d ms\n",round(preFixationTime(trial_count)));
            end
            
            
%             %if isReleaseToHold == 1 && preFixationEndFlag == 0  
%             if isHold == 1 && preFixationEndFlag == 0 
%                 isWaitGocueBreak(trial_count) = 1;
%                 isItiInitial_break = 1;
%                 isShowingInitial = 0;
%             end

%             if TouchEventAvail(dev) > 0
%                 isWaitGocueBreak(trial_count) = 1;
%                 isItiInitial_break = 1;
%                 isShowingInitial = 0;
%                 t0_isShowing = tic;
%                 continue
%             end
           
        elseif if_freeview0_fixation1 == 0
            preFixationEndFlag = 1;
        end
        vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);
        

        
        %if if_fix_noHold0_hold1 == 0
        %   isHold = 1;
        %end
        %if isReleaseToHold == 1 && isItiInitial_break == 0
        if isHold == 1 && isItiInitial_break == 0 && preFixationEndFlag == 1
            if if_eyeResponse == 0
                while TouchEventAvail(dev) > 0
                    TouchEventGet(dev, window);
                end
            end
            if if_freeview0_fixation1 == 0 || toc(t0_fixOn) > preFixationTime(trial_count)/1000
                
                preFixationEndFlag = 0;
                
                holdEvtMappedX = 0;
                holdEvtMappedY = 0;
                holdInside = 0;

                isShowingInitial = 0;
                isShowing = 1;
                isShowingReady = 0;
                showing_count = 1;
                t0_isShowing = tic;
                newTouchPoint = 0;
                timeFlag1 = 0;
                
                isTriangletrialMiddle = 0;          

                m_marker.Data = uint8(1);% The trial begin marker 
                trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms
                %Eyelink('Message', 'TRIALID %d', trial_count);
                message = sprintf('TRIALID %d', trial_count);
                Eyelink('Message', message);                
                trial_para.marker.count = trial_para.marker.count + 1;
                m_markerCount.Data = int32(trial_para.marker.count);
                trial_para.marker.content = [trial_para.marker.content, {message}];                
                flag_targetOn_marker = 0;
            end
        end     
        
        last_preFixationEndFlag = preFixationEndFlag;        
    %---------------------------------------------------------------------------------------------   
    elseif isShowing == 1
        %% isShowing
        if if_fix_noHold0_hold1 == 0                 
            isHold = 1 - isHold;
            if if_eyeResponse == 0
                if TouchEventAvail(dev) > 0
                    isHold = 0;
                end
            end
            
%             [mx, my, buttons] = GetMouse(window);
%             holdInside = max(IsInRect(mx, my, holdRects_dection'), holdInside); %#ok<*NASGU>
%             while TouchEventAvail(dev) > 0
%                 evt = TouchEventGet(dev, window);
% 
%                 % for touchhold
%                 if evt.Type == 2 || evt.Type == 3 || evt.Type == 4
%                     holdEvtMappedX = evt.MappedX;
%                     holdEvtMappedY = evt.MappedY;
%                 end
%                 holdInside = max( IsInRect(holdEvtMappedX, holdEvtMappedY, holdRects_dection'), holdInside );
% 
%                 if evt.Type == 4 
%                     holdEvtMappedX = 0;
%                     holdEvtMappedY = 0;
%                 end
%             end            
%             holdInside = max( IsInRect(holdEvtMappedX, holdEvtMappedY, holdRects_dection'), holdInside );            
%             isHold = holdInside;
%             Screen('FrameRect', window, frameColor, holdRects, 3);
        end
        
        if if_freeview0_fixation1 == 1 && fixOn_flag == 0         
            isHold = 0;
        end         
        if isHold == 0            
            %isCorrect(trial_count) = -1;
            isWaitGocueBreak(trial_count) = 1;
            isItiInitial_break = 1;
            isShowing = 0;
            continue
        end
        
        if toc(t0_isShowing) >= 0 && isShowingReady == 0

            if seq_length_rangeTail > 0
%                 if toc < 0.3
                if toc(t0_isShowing) < 0
                    if ifCircle == 1
                        Screen('FrameOval', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
                        
                        Screen('FillRect', window, [1 1 1]*0.75, centeredVerticalCrossRect);
                        Screen('FillRect', window, [1 1 1]*0.75, centeredHorizontalCrossRect);
                        
%                         if if_freeview0_fixation1 == 1
%                             if fixOn_flag == 1
%                                 Screen('DrawDots', window, [x_avg y_avg], 20, [255 0 0], [], 2);
%                             elseif fixOn_flag == 0
%                                 Screen('DrawDots', window, [x_avg y_avg], 20, [0 255 0], [], 2);
%                             end
%                         end
                        
                        vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);
                        
                    else
                        
                        if taskMode_noChoice_Green0_Red1 == 0
                            if ifStep0 == 0
                                Screen('FrameRect', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), 3);
                            end
                        elseif taskMode_noChoice_Green0_Red1 == 1
                            for i=1:length(currentFrames)
                                Screen('FramePoly', window, frameColorMatrices(:, currentFrames(i)), showTriangles{currentFrames(i)}, 1);
                            end
                        end
                        
                        Screen('FillRect', window, [1 1 1]*0.75, centeredVerticalCrossRect);
                        Screen('FillRect', window, [1 1 1]*0.75, centeredHorizontalCrossRect);                        
                        
                        vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);
                    end
                else
                    isShowingReady = 1;
                    t0_showReady = tic;
                    t0_smallWater = tic;
                end
            elseif  seq_length_rangeTail == 0
                isShowingReady = 1;
                t0_showReady = tic;
                t0_smallWater = tic;
            end

        end
        
        currentColorMatrices = repmat(specialColor',1,numSquares);
        if isShowingReady == 1  
            t1_showReady = toc(t0_showReady);
            t1_smallWater = toc(t0_smallWater);                 
            if t1_smallWater > smallWaterInterval/1000 && ifSmallWater == 1
                giveWater(smallWater);
                t0_smallWater = tic;
            end
            if ifCircle == 0
                if taskMode_noChoice_Green0_Red1 == 0
                    if isempty(currentSequence) == 0
                        if ifStep0 == 0
                            if tempStimulusExtendTime == 0 && ...
                                    t1_showReady >= tempPointShowTime/1000 * (showing_count-1+pointShowPWM)
                                Screen('FillRect', window, currentColorMatrices(:, currentSequence(showing_count)), showRects(:, currentSequence(showing_count)));
                            else                                
                                %Screen('FillRect', window, currentColorMatrices(:, currentSequence(showing_count)), showRects(:, currentSequence(showing_count)));
                            end
                        end
                    end
                    if ifStep0 == 0
                        Screen('FrameRect', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), 3);
                    end
                elseif taskMode_noChoice_Green0_Red1 == 1
                    if tempStimulusExtendTime == 0 && ...
                            t1_showReady >= tempPointShowTime/1000 * (showing_count-1+pointShowPWM)
                        Screen('FillPoly', window, currentColorMatrices(:, currentSequence(showing_count))', showTriangles{currentSequence(showing_count)}, 1);                        
                    else
                        %Screen('FillPoly', window, currentColorMatrices(:, currentSequence(showing_count))', showTriangles{currentSequence(showing_count)}, 1);
                    end
                    for i=1:length(currentFrames)
                        Screen('FramePoly', window, frameColorMatrices(:, currentFrames(i)), showTriangles{currentFrames(i)}, 1);
                    end
                end
            elseif ifCircle == 1
                if tempStimulusExtendTime == 0 && ...
                        t1_showReady >= tempPointShowTime/1000 * (showing_count-1+pointShowPWM)
                    Screen('FillOval', window, currentColorMatrices(:, currentSequence(showing_count)), showRects(:, currentSequence(showing_count)));                    
                    if flag_targetOn_marker == 0
                       flag_targetOn_marker = 1; 
                       m_marker.Data = uint8(2);% Inner trial markers
                       trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms
                       %Eyelink('Message', 'TARGET %d ON', showing_count);
                       %message = sprintf('TARGET %d ON', showing_count);
                       message = sprintf('TARGET %d ITEM %d ON', showing_count, currentSequence(showing_count));
                       Eyelink('Message', message);
                       trial_para.marker.count = trial_para.marker.count + 1;
                       m_markerCount.Data = int32(trial_para.marker.count);
                       trial_para.marker.content = [trial_para.marker.content, {message}];
                       
                    end                    
                else
                    flag_targetOn_marker = 0;
                    %Screen('FillOval', window, currentColorMatrices(:, currentSequence(showing_count)), showRects(:, currentSequence(showing_count)));
                end
                                
                Screen('FrameOval', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
            end
            Screen('FillRect', window, [1 1 1]*0.75, centeredVerticalCrossRect);
            Screen('FillRect', window, [1 1 1]*0.75, centeredHorizontalCrossRect);
            
            if ifDrawDots_eye == 1
                if if_freeview0_fixation1 == 1
                    if fixOn_flag == 1
                        %Screen('DrawDots', window, [x_avg y_avg], 20, [255 0 0], [], 2);
                        Screen('DrawDots', window, [x_avg y_avg], 20, [0.5 0.5 0.5], [], 2);
                    elseif fixOn_flag == 0
                        %Screen('DrawDots', window, [x_avg y_avg], 20, [0 255 0], [], 2);
                    end
                end
            end
            
            if tempStimulusExtendTime == 0
                vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);
            elseif tempStimulusExtendTime > 0
                vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi, 1);
            end
  
            
            
            if t1_showReady >= tempPointShowTime/1000 * showing_count || seq_length(trial_count) == 0
                if showing_count >= seq_length(trial_count)   
                    holdEvtMappedX = 0;
                    holdEvtMappedY = 0;
                    holdInside = 0;

                    isShowing = 0;
                    isDelay = 1;
                    choiceCondition_flag(trial_count) = 2;
                    t0_delay = tic;
                    %t0_smallWater = tic;
                    m_marker.Data = uint8(2);% Inner trial markers
                    trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms
                    %Eyelink('Message', 'DELAY 1 ON');
                    message = sprintf('DELAY 1 ON');
                    Eyelink('Message', message);
                    trial_para.marker.count = trial_para.marker.count + 1;
                    m_markerCount.Data = int32(trial_para.marker.count);
                    trial_para.marker.content = [trial_para.marker.content, {message}];
                    flag_LongWater1 = 0;
                    flag_LongWater2 = 0;
                    flag_LongWater3 = 0;
                    flag_LongWater4 = 0;                    
                end
                showing_count = showing_count + 1;
            end
            
        end
    %---------------------------------------------------------------------------------------------    
    elseif isDelay == 1 
        %% isDelay
        if if_fix_noHold0_hold1 == 0                 
            isHold = 1 - isHold;
            if if_eyeResponse == 0
                if TouchEventAvail(dev) > 0
                    isHold = 0;
                end
            end
            
%             [mx, my, buttons] = GetMouse(window);
%             holdInside = max(IsInRect(mx, my, holdRects_dection'), holdInside); %#ok<*NASGU>
%             while TouchEventAvail(dev) > 0
%                 evt = TouchEventGet(dev, window);
% 
%                 % for touchhold
%                 if evt.Type == 2 || evt.Type == 3 || evt.Type == 4
%                     holdEvtMappedX = evt.MappedX;
%                     holdEvtMappedY = evt.MappedY;
%                 end
%                 holdInside = max( IsInRect(holdEvtMappedX, holdEvtMappedY, holdRects_dection'), holdInside );
% 
%                 if evt.Type == 4 
%                     holdEvtMappedX = 0;
%                     holdEvtMappedY = 0;
%                 end
%             end            
%             holdInside = max( IsInRect(holdEvtMappedX, holdEvtMappedY, holdRects_dection'), holdInside );            
%             isHold = holdInside;
%             Screen('FrameRect', window, frameColor, holdRects, 3);
        end
        
        if if_freeview0_fixation1 == 1 && fixOn_flag == 0         
            isHold = 0;
        end          
        if isHold == 0
            %isCorrect(trial_count) = -1;
            isWaitGocueBreak(trial_count) = 1;            
            isItiInitial_break = 1;
            isDelay = 0;
            continue
        end
        
        t1_smallWater = toc(t0_smallWater);
        if t1_smallWater > smallWaterInterval/1000 && ifSmallWater == 1
            tempSmallWater = smallWater;
            if toc(t0_isShowing) > 8000/1000
                tempSmallWater = 2 * smallWater;
            end
            %giveWater(tempSmallWater);
            t0_smallWater = tic;  
            temp_flag = 0;
            if ifLongWater == 1
                if toc(t0_isShowing) > 7000/1000 && flag_LongWater1 == 0
                    flag_LongWater1 = 1;
                    fprintf('\nGood fixation!\n');
                    giveWater(delay1Water);     
                    temp_flag = 1;
                end
                if toc(t0_isShowing) > 11000/1000 && flag_LongWater2 == 0
                    flag_LongWater2 = 1;                    
                    fprintf('\nVery good fixation!\n');
                    giveWater(delay1Water*1.3);
                    temp_flag = 1;
                end     
                if toc(t0_isShowing) > 15000/1000 && flag_LongWater3 == 0
                    flag_LongWater3 = 1;                    
                    fprintf('\nExcellent fixation!\n');
                    giveWater(delay1Water*1.8);
                    temp_flag = 1;
                end        
                if toc(t0_isShowing) > 19000/1000 && flag_LongWater4 == 0
                    flag_LongWater4 = 1;                    
                    fprintf('\nBest fixation!\n');
                    giveWater(delay1Water*2.5);
                    temp_flag = 1;
                end                   
            end
            if temp_flag == 0
                giveWater(tempSmallWater);                
            end
        end
        
        if cumulativeError >= encourageCumuErrorLimit
            tempStimET = encourageStimET;
        else
            tempStimET = tempStimulusExtendTime;
        end
        
%         if toc(t0_delay) < tempStimET/1000   
%             if ifCircle == 0
%                 if taskMode_noChoice_Green0_Red1 == 0
%                     if isempty(currentSequence) == 0
%                         Screen('FillRect', window, currentColorMatrices(:, currentSequence), showRects(:, currentSequence));
%                     end
%                     if ifStep0 == 0
%                         Screen('FrameRect', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), 3);
%                     end
%                 elseif taskMode_noChoice_Green0_Red1 == 1
%                     for i=1:length(currentSequence)
%                         Screen('FillPoly', window, currentColorMatrices(:, currentSequence(i)), showTriangles{currentSequence(i)}, 1);
%                     end
%                     for i=1:length(currentFrames)
%                         Screen('FramePoly', window, frameColorMatrices(:, currentFrames(i)), showTriangles{currentFrames(i)}, 1);
%                     end
%                 end
%             elseif ifCircle == 1
%                 Screen('FillOval', window, currentColorMatrices(:, currentSequence), showRects(:, currentSequence));
%                 Screen('FrameOval', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
%             end
% 
%         end
        
        currentFixationTime = fixationTime(trial_count);
        if ifSelection == 0
            if toc(t0_delay) < 0.5*currentFixationTime/1000
                Screen('FillRect', window, [1 1 1]*0.75, centeredVerticalCrossRect);
                Screen('FillRect', window, [1 1 1]*0.75, centeredHorizontalCrossRect);
            elseif toc(t0_delay) > 0.5*currentFixationTime/1000 && toc(t0_delay) < currentFixationTime/1000
                if taskMode_noChoice_Green0_Red1 == 0
                    Screen('FillRect', window, [0 0.5 0 ], centeredFixSquare);
                elseif taskMode_noChoice_Green0_Red1 == 1
                    Screen('FillPoly', window, [0.5 0 0 ], [fixTri_xPosVector; fixTri_yPosVector]', 1);
                end
            end
        elseif ifSelection == 1
            if toc(t0_delay) < currentFixationTime/1000
                Screen('FillRect', window, [1 1 1]*0.75, centeredVerticalCrossRect);
                Screen('FillRect', window, [1 1 1]*0.75, centeredHorizontalCrossRect);
                
                if ifDrawDots_eye == 1
                    if if_freeview0_fixation1 == 1
                        if fixOn_flag == 1
                            %Screen('DrawDots', window, [x_avg y_avg], 20, [255 0 0], [], 2);
                            Screen('DrawDots', window, [x_avg y_avg], 20, [0.5 0.5 0.5], [], 2);
                        elseif fixOn_flag == 0
                            %Screen('DrawDots', window, [x_avg y_avg], 20, [0 255 0], [], 2);
                        end
                    end
                end
            end
        end

        
        if ifCircle == 0
            if taskMode_noChoice_Green0_Red1 == 0
                if ifStep0 == 0
                    Screen('FrameRect', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), 3);
                end
            elseif taskMode_noChoice_Green0_Red1 == 1
                for i=1:length(currentFrames)
                    Screen('FramePoly', window, frameColorMatrices(:, currentFrames(i)), showTriangles{currentFrames(i)}, 1);
                end
            end
        elseif ifCircle == 1
%             if consecutiveCorrect >= 1
%                 randArray7(trial_count) = 0;
%             else
%                 randArray7(trial_count) = 1;
%             end
            %randArray7(trial_count) = 1;
            randArray7(trial_count) = 0;
            if toc(t0_delay) < pointShowTime*pointShowPWM/1000 || randArray7(trial_count) > 0.5 || if_justFix0_touch1 == 0
                %Screen('FrameOval', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
                %Screen('FrameOval', window, frameColorMatrices2(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);                
                
                frameColor;
                frameColor2;
                %temp_holdOnTime = 0/1000;
                %temp_holdOnTime = 200/1000;%200 ms
                temp_holdOnTime = (pointShowTime*pointShowPWM-frameDisappearTime)/1000;
                if toc(t0_delay) < temp_holdOnTime
                    temp_W = 0;
                else
                    temp_W = (toc(t0_delay)-temp_holdOnTime) / (pointShowTime*pointShowPWM/1000-temp_holdOnTime);
                end
                %temp_W = toc(t0_delay) / (pointShowTime*pointShowPWM/1000);
                frameColorX = (1-temp_W) .* frameColor + temp_W .* frameColor2;
                frameColorMatricesX = repmat(frameColorX',1,numSquares);
                Screen('FrameOval', window, frameColorMatricesX(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
            end
            if toc(t0_delay) >= pointShowTime*pointShowPWM/1000
                Screen('FrameOval', window, frameColorMatrices2(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
            end
        end
        
        if toc(t0_delay) > currentFixationTime/1000
            Screen('FillRect', window, [1 0 1], centeredVerticalCrossRect);
            Screen('FillRect', window, [1 0 1], centeredHorizontalCrossRect);            
            
            if if_eyeResponse == 0
                while TouchEventAvail(dev) > 0
                    TouchEventGet(dev, window);
                end
            end
            fprintf("\ncurrentFixationTime = %d ms\n",round(currentFixationTime));
%             fprintf("tempVar * currentFixationTime = %d\n",(tempVar * currentFixationTime));
            if taskMode_noChoice_Green0_Red1 == 0
                %Screen('FillRect', window, [0 0.5 0 ], centeredFixSquare);
            elseif taskMode_noChoice_Green0_Red1 == 1
                %Screen('FillPoly', window, [0.5 0 0 ], [fixTri_xPosVector; fixTri_yPosVector]', 1);
                
                %%%%%OffloadBeginEnd%%%%%%
                if taskMode_noChoice_Green0_Red1 == 1
                    sortedCurrentSeq = sort(currentSequence);%let begin=left, end=right in random stimulus condition
                    if ifBeginEnd == 1
                        currentMiddleSeq = sortedCurrentSeq(2:end-1);
                        currentBeginEndSeq = [ sortedCurrentSeq(1) sortedCurrentSeq(end) ];
                    elseif ifBeginEnd == 0
                        currentMiddleSeq = [];
                        currentBeginEndSeq = sortedCurrentSeq;
                    end
                    if ~isempty(currentMiddleSeq)
                        if ifCircle == 0
                            for i=1:length(currentMiddleSeq)                                
                                Screen('FillPoly', window, currentColorMatrices(:, currentMiddleSeq(i)), showTriangles{currentMiddleSeq(i)}, 1);                                
                            end
                        elseif ifCircle == 1
                            if ifSelection == 0
                                Screen('FillOval', window, currentColorMatrices(:, currentMiddleSeq), showRects(:, currentMiddleSeq));
                            end
                        end
                    end
                    temp_isFill_OffloadBeginEnd = zeros(1, length(isFill));
                    temp_isFill_OffloadBeginEnd(currentBeginEndSeq) = 1;
                    %temp_isFill = temp_isFill_OffloadBeginEnd;
                    isFill = temp_isFill_OffloadBeginEnd;
                end                
                
            end
            if ifSelection == 1
                if selectFlag_newSeq(trial_count) == 0
                    isGreenRectLeft(trial_count) = isGreenRectLeft(trial_count - 1);
                    choiceCondition_flag(trial_count) = choiceCondition_flag(trial_count - 1);                                       
                    if ifSelectOffloading == 0
                        choiceMode = 0;
                    elseif ifSelectOffloading == 1
                        choiceMode = 1;
                    end
                end
                if selectFlag_newSeq(trial_count) == 1
                    if choiceMode_backup ~= 2
                        choiceMode = choiceMode_backup;
                    elseif choiceMode_backup == 2
                        if if_seqList == 1 && ifFormalTask == 1
                            choiceMode = choiceCondition_flag_list(trial_count);
                            choiceCondition_flag(trial_count) = choiceCondition_flag_list(trial_count);
                            
                            if if_choiceToR == 1
                                if choiceCondition_flag_list(trial_count) == 2
                                    choiceMode = 1;
                                    choiceCondition_flag(trial_count) = 1;
                                end
                            end
                            
                            if choiceMode == 0
                                fprintf("No choice condition!\n");
                            elseif choiceMode == 2
                                fprintf("Choice condition!\n");
                            end
                        else
                            %randArray2
                            fprintf("Current randArray2 = %.2f, ", randArray2(trial_count));
                            if randArray2(trial_count) <= choiceCondition_boardLine1
                                choiceMode = 0;
                                choiceCondition_flag(trial_count) = 0;
                                fprintf("g/r random --> g\n");
                            elseif randArray2(trial_count) <= choiceCondition_boardLine2
                                choiceMode = 1;
                                choiceCondition_flag(trial_count) = 1;
                                fprintf("g/r random --> r\n");
                            else
                                choiceMode = 2;
                                choiceCondition_flag(trial_count) = 2;
                                fprintf("g/r random --> g/r\n");
                            end
                        end
                    end
                end
                
                %selectionTriangles
                if isGreenRectLeft(trial_count) == 1
                    if choiceMode == 0
                        Screen('FillRect', window, normalColor, selectionRects(:,1));
                    elseif choiceMode == 1
                        Screen('FillPoly', window, specialColor_decision, selectionTriangles{2}, 1);
                    elseif choiceMode == 2
                        Screen('FillRect', window, normalColor, selectionRects(:,1));
                        Screen('FillPoly', window, specialColor_decision, selectionTriangles{2}, 1);
                    end
                    
                else
                    if choiceMode == 0
                        Screen('FillRect', window, normalColor, selectionRects(:,2));
                    elseif choiceMode == 1
                        Screen('FillPoly', window, specialColor_decision, selectionTriangles{1}, 1);
                    elseif choiceMode == 2
                        Screen('FillPoly', window, specialColor_decision, selectionTriangles{1}, 1);
                        Screen('FillRect', window, normalColor, selectionRects(:,2));
                    end
                    
                end                
                
            end
            
            holdEvtMappedX = 0;
            holdEvtMappedY = 0;
            holdInside = 0;

            isDelay = 0;
            %isWaitGocueBreak(trial_count) = -1;
            isReleaseForTouch = 1;
            tempFlag2 = 0;
            
            %giveWater(90);
            
%             if preFixationTime(trial_count) < 0

            if if_justFix0_touch1 == 0
                isReleaseForTouch = 0;
                isItiInitial = 1;
                tempIti = ItiTime(trial_count);
                isCorrect(trial_count) = 1;
                ifSelectOffloading = 1;
                endingHold_RT = 0;                
            else
%                 if ifDelay1OffWater == 1
%                     giveWater(delay1Water);
%                 end
            end
            if ifDelay1OffWater == 1
                giveWater(delay1Water);
            end
            
                
%             else
%                 giveWater(90);
%             end
            
            
            t0_isReleaseForTouch = tic;
            
            %Screen('FrameOval', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
            Screen('FrameOval', window, frameColorMatrices2(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);            
        end
        


        vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);        
    %---------------------------------------------------------------------------------------------    
    elseif isReleaseForTouch == 1
        %% isReleaseForTouch
        if if_fix_noHold0_hold1 == 0
            isHold = 0;
            
%             [mx, my, buttons] = GetMouse(window);
%             while TouchEventAvail(dev) > 0
%                 evt = TouchEventGet(dev, window);
% 
%                 % for touchhold
%                 if evt.Type == 2 || evt.Type == 3
%                     holdEvtMappedX = evt.MappedX;
%                     holdEvtMappedY = evt.MappedY;
%                 end
% 
%                 holdInside = IsInRect(mx, my, holdRects_dection'); %#ok<*NASGU>
%                 holdInside = max( IsInRect(holdEvtMappedX, holdEvtMappedY, holdRects_dection'), holdInside );
% 
%                 holdEvtMappedX = 0;
%                 holdEvtMappedY = 0;
%             end
% 
%             isHold = holdInside;
%             Screen('FrameRect', window, frameColor, holdRects, 3);            
        end
        
        if isHold == 1 && tempFlag2 == 0
            fprintf('isReleaseForTouch=%d, isHold=%d\n',isReleaseForTouch,isHold);
            tempFlag2 = 1;
        end
        if isHold == 0
            fprintf('isReleaseForTouch=%d, isHold=%d\n',isReleaseForTouch,isHold);
        end
        
        if isHold == 0
            isReleaseForTouch = 0;
            
            m_marker.Data = uint8(2);% Inner trial markers
            trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms
            %Eyelink('Message', 'RELEASE 1 ON');
            %message = sprintf('RELEASE 1 ON');
            message = sprintf('DELAY 1 OFF');
            Eyelink('Message', message);
            trial_para.marker.count = trial_para.marker.count + 1;
            m_markerCount.Data = int32(trial_para.marker.count);
            trial_para.marker.content = [trial_para.marker.content, {message}];
            
            if ifSelection == 0
                
                newTouchEnableFlag = 1;
                tempKeycode = 0;       

                isTouching = 1;
                lastTouchKeycode = 0;
                isTouchingDone = 0;
                isTouching_tempFlag1 = 0;
                isTouching_tempFlag2 = 0;
                last_freeTouch_toc = 0;
                tempFlag_toc = 0;
                


            elseif ifSelection == 1
                isSelecting = 1;
                raw_selectInside = zeros(1, 2);
                last_selectInside = zeros(1, 2);              
                
                newTouchEnableFlag = 1;
                tempKeycode = 0;                                
                currentKeycode = 1;
                %tic
            end
            newTouchPoint = [];

            %tic
        end

        if ifTurnOff == 1 && ifSelection == 0
            if taskMode_noChoice_Green0_Red1 == 1
                tempStimET = encourageStimET;
            elseif taskMode_noChoice_Green0_Red1 == 0
                tempStimET = tempStimulusExtendTime;
            end
        end
        
        if if_fix_noHold0_hold1 == 1
            if toc(t0_isReleaseForTouch) < (tempStimET - fixationTime(trial_count))/1000
                if ifCircle == 0
                    if taskMode_noChoice_Green0_Red1 == 0
                        Screen('FillRect', window, currentColorMatrices(:, currentSequence), showRects(:, currentSequence));
                    elseif taskMode_noChoice_Green0_Red1 == 1
                        for i=1:length(currentSequence)
                            Screen('FillPoly', window, currentColorMatrices(:, currentSequence(i)), showTriangles{currentSequence(i)}, 1);
                        end
                    end
                elseif ifCircle == 1
                    Screen('FillOval', window, currentColorMatrices(:, currentSequence), showRects(:, currentSequence));
                end     

                if taskMode_noChoice_Green0_Red1 == 0
                    %Screen('FillRect', window, [0 0.5 0 ], centeredFixSquare);
                elseif taskMode_noChoice_Green0_Red1 == 1
                    %Screen('FillPoly', window, [0.5 0 0 ], [fixTri_xPosVector; fixTri_yPosVector]', 1);

                    %%%%%OffloadBeginEnd%%%%%%
                    if taskMode_noChoice_Green0_Red1 == 1
                        sortedCurrentSeq = sort(currentSequence);%let begin=left, end=right in random stimulus condition
                        if ifBeginEnd == 1
                            currentMiddleSeq = sortedCurrentSeq(2:end-1);
                            currentBeginEndSeq = [ sortedCurrentSeq(1) sortedCurrentSeq(end) ];
                        elseif ifBeginEnd == 0
                            currentMiddleSeq = [];
                            currentBeginEndSeq = sortedCurrentSeq;
                        end
                        if ~isempty(currentMiddleSeq)
                            if ifCircle == 0
                                for i=1:length(currentMiddleSeq)                                
                                    Screen('FillPoly', window, currentColorMatrices(:, currentMiddleSeq(i)), showTriangles{currentMiddleSeq(i)}, 1);                                
                                end
                            elseif ifCircle == 1
                                if ifSelection == 0
                                    Screen('FillOval', window, currentColorMatrices(:, currentMiddleSeq), showRects(:, currentMiddleSeq));
                                end
                            end
                        end
                        temp_isFill_OffloadBeginEnd = zeros(1, length(isFill));
                        temp_isFill_OffloadBeginEnd(currentBeginEndSeq) = 1;
                        %temp_isFill = temp_isFill_OffloadBeginEnd;
                        isFill = temp_isFill_OffloadBeginEnd;
                    end


                end

                if ifCircle == 0
                    if taskMode_noChoice_Green0_Red1 == 0
                        if ifStep0 == 0
                            Screen('FrameRect', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), 3);
                        end
                    elseif taskMode_noChoice_Green0_Red1 == 1
                        for i=1:length(currentFrames)
                            Screen('FramePoly', window, frameColorMatrices(:, currentFrames(i)), showTriangles{currentFrames(i)}, 1);
                        end
                    end
                elseif ifCircle == 1
                    Screen('FrameOval', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
                end

                if ifSelection == 1
                    %selectionTriangles
                    if isGreenRectLeft(trial_count) == 1
                        if choiceMode == 0
                            Screen('FillRect', window, normalColor, selectionRects(:,1));
                        elseif choiceMode == 1
                            Screen('FillPoly', window, specialColor_decision, selectionTriangles{2}, 1);
                        elseif choiceMode == 2
                            Screen('FillRect', window, normalColor, selectionRects(:,1));
                            Screen('FillPoly', window, specialColor_decision, selectionTriangles{2}, 1);
                        end

                    else
                        if choiceMode == 0
                            Screen('FillRect', window, normalColor, selectionRects(:,2));
                        elseif choiceMode == 1
                            Screen('FillPoly', window, specialColor_decision, selectionTriangles{1}, 1);
                        elseif choiceMode == 2
                            Screen('FillPoly', window, specialColor_decision, selectionTriangles{1}, 1);
                            Screen('FillRect', window, normalColor, selectionRects(:,2));
                        end

                    end                
                end

                vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);  
            end 
        end
    %---------------------------------------------------------------------------------------------           
    elseif isSelecting == 1       
        %% isSelecting
        if toc(t0_isReleaseForTouch) > selecting_RT_limit/1000
            isWaitGocueBreak(trial_count) = 1;
            isItiInitial_break = 1;
            isSelecting = 0;
            continue
        end
        
        Screen('FrameOval', window, frameColorMatrices2(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
        
        % Get the current position of the mouse
        [mx, my, buttons] = GetMouse(window);          
                
        tempEvtMappedX = 0;
        tempEvtMappedY = 0;
        if if_eyeResponse == 0
            while TouchEventAvail(dev) > 0
                evt = TouchEventGet(dev, window);
                %evt.Type == 2 means new touch point
                %evt.Type == 3 means moving touch point
                %evt.Type == 4 means touch finish
                if evt.Type == 4
                    if evt.Keycode == tempKeycode
                        newTouchEnableFlag = 1;
                    end
                end
                if evt.Type == 4
                    %if evt.Type == 2
                    %if evt.Type == 2 || evt.Type == 3
                    %if evt.Type == 2 || evt.Type == 3
                    % Draw a white dot where the touch point is
                    %Screen('DrawDots', window, [evt.MappedX evt.MappedY], 50, white, [], 2);
                    tempEvtMappedX = evt.MappedX;
                    tempEvtMappedY = evt.MappedY;
                    
                    if isempty(mytouchQueueBuffer_selection{1}) == 1 || touchBufferSize == 1
                        mytouchQueueBuffer_selection{1} = [evt.MappedX evt.MappedY];
                    elseif isempty(mytouchQueueBuffer_selection{1}) == 0
                        tempDis = norm([evt.MappedX evt.MappedY] - mytouchQueueBuffer_selection{1});
                        if tempDis > selectbaseRect_detection(4)*2
                            mytouchQueueBuffer_selection = cell(1, touchBufferSize);
                            mytouchQueueBuffer_selection{1} = [evt.MappedX evt.MappedY];
                        else
                            if isempty(mytouchQueueBuffer_selection{touchBufferSize}) == 1
                                for i=2:touchBufferSize
                                    if isempty(mytouchQueueBuffer_selection{i}) == 1
                                        for j=i:-1:2
                                            mytouchQueueBuffer_selection{j} = mytouchQueueBuffer_selection{j-1};
                                        end
                                        mytouchQueueBuffer_selection{1} = [evt.MappedX evt.MappedY];
                                        break
                                    end
                                end
                            elseif isempty(mytouchQueueBuffer_selection{touchBufferSize}) == 0
                                for j=touchBufferSize:-1:2
                                    mytouchQueueBuffer_selection{j} = mytouchQueueBuffer_selection{j-1};
                                end
                                mytouchQueueBuffer_selection{1} = [evt.MappedX evt.MappedY];
                            end
                        end
                    end
                    tempTouchBufferLen = 0;
                    for i=1:touchBufferSize
                        if isempty(mytouchQueueBuffer_selection{i}) == 0
                            tempTouchBufferLen = tempTouchBufferLen + 1;
                        end
                    end
                    
                    if newTouchEnableFlag == 1
                        % Draw a white dot where the touch point is
                        % Screen('DrawDots', window, [evt.MappedX evt.MappedY], 50, white, [], 2);
                        %newTouchEnableFlag = 0;
                        tempKeycode = evt.Keycode;
                    end
                    currentKeycode = evt.Keycode;
                end
            end
        end
            
        %selectionTriangles
        if isGreenRectLeft(trial_count) == 1
            if choiceMode == 0
                Screen('FillRect', window, normalColor, selectionRects(:,1));
            elseif choiceMode == 1
                Screen('FillPoly', window, specialColor_decision, selectionTriangles{2}, 1);
            elseif choiceMode == 2
                Screen('FillRect', window, normalColor, selectionRects(:,1));
                Screen('FillPoly', window, specialColor_decision, selectionTriangles{2}, 1);
            end
            
        else
            if choiceMode == 0
                Screen('FillRect', window, normalColor, selectionRects(:,2));
            elseif choiceMode == 1
                Screen('FillPoly', window, specialColor_decision, selectionTriangles{1}, 1);
            elseif choiceMode == 2
                Screen('FillPoly', window, specialColor_decision, selectionTriangles{1}, 1);
                Screen('FillRect', window, normalColor, selectionRects(:,2));
            end
            
        end
        %end
        
        %Screen('FrameRect', window, normalColor, selectionRects_detection);
        
        %Screen('FrameOval', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
        %Screen('FrameOval', window, frameColorMatrices2(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
        
        %if ifDrawDots_eye == 1 && if_eyeResponse == 1
        if if_eyeResponse == 1
            %Screen('DrawDots', window, [x_avg y_avg], 20, [0.5 0.5 0.5], [], 2);
        end        
        vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);
        
        tempBufferMeanXpos_selection = 0;
        tempBufferMeanYpos_selection = 0;
        if isempty(mytouchQueueBuffer_selection{1}) == 0
            for tempj=1:tempTouchBufferLen
                tempBufferMeanXpos_selection = tempBufferMeanXpos_selection + mytouchQueueBuffer_selection{tempj}(1);
                tempBufferMeanYpos_selection = tempBufferMeanYpos_selection + mytouchQueueBuffer_selection{tempj}(2);
            end
            tempBufferMeanXpos_selection = tempBufferMeanXpos_selection/tempTouchBufferLen;
            tempBufferMeanYpos_selection = tempBufferMeanYpos_selection/tempTouchBufferLen;
        end
        
        if touchBuffer_enable == 0
            tempBufferMeanXpos_selection = 0;
            tempBufferMeanYpos_selection = 0;
        end
        
        %tempBufferMeanXpos_selection tempBufferMeanYpos_selection
        selectInside = [0 0];
        for i = 1:2
            if choiceMode == 0
                if isGreenRectLeft(trial_count) == 1
                    if i == 1
                        if if_mouse_mannual_eyeResponse == 0
                            selectInside(i) = IsInRect(mx, my, selectionRects_detection(:, i)');
                        end
                        if currentKeycode == tempKeycode
                            selectInside(i) = max( IsInRect(tempEvtMappedX, tempEvtMappedY, selectionRects_detection(:, i)'), selectInside(i) );
                        end
                        selectInside(i) = max( IsInRect(tempBufferMeanXpos_selection, tempBufferMeanYpos_selection,...
                            selectionRects_detection(:, i)'), selectInside(i) );
                        if if_eyeResponse == 1
                            %selectInside(i) = max( IsInRect(x_avg, y_avg,...
                            %    selectionRects_detection(:, i)'), selectInside(i) );   
                            
                            raw_selectInside(i) = IsInRect(x_avg, y_avg, selectionRects_detection(:, i)');         
                            if raw_selectInside(i) == 1 && last_selectInside(i) == 1

                            else
                                t0_selectHolding = tic;
                            end              
                            if raw_selectInside(i) == 1         
                                if toc(t0_selectHolding) > eyeResponseTime/1000
                                    selectInside(i) = 1;
                                end
                            end
                            last_selectInside(i) = raw_selectInside(i);                            
                        end
                    end
                elseif isGreenRectLeft(trial_count) == 0
                    if i == 2
                        if if_mouse_mannual_eyeResponse == 0
                            selectInside(i) = IsInRect(mx, my, selectionRects_detection(:, i)');
                        end
                        if currentKeycode == tempKeycode
                            selectInside(i) = max( IsInRect(tempEvtMappedX, tempEvtMappedY, selectionRects_detection(:, i)'), selectInside(i) );
                        end
                        selectInside(i) = max( IsInRect(tempBufferMeanXpos_selection, tempBufferMeanYpos_selection,...
                           selectionRects_detection(:, i)'), selectInside(i) );  
                        if if_eyeResponse == 1
                            %selectInside(i) = max( IsInRect(x_avg, y_avg,...
                            %    selectionRects_detection(:, i)'), selectInside(i) ); 
                            raw_selectInside(i) = IsInRect(x_avg, y_avg, selectionRects_detection(:, i)');         
                            if raw_selectInside(i) == 1 && last_selectInside(i) == 1

                            else
                                t0_selectHolding = tic;
                            end              
                            if raw_selectInside(i) == 1         
                                if toc(t0_selectHolding) > eyeResponseTime/1000
                                    selectInside(i) = 1;
                                end
                            end
                            last_selectInside(i) = raw_selectInside(i);                               
                        end                       
                    end
                end
            elseif choiceMode == 1
                if isGreenRectLeft(trial_count) == 1
                    if i == 2
                        if if_mouse_mannual_eyeResponse == 0
                            selectInside(i) = IsInRect(mx, my, selectionRects_detection(:, i)');
                        end
                        if currentKeycode == tempKeycode
                            selectInside(i) = max( IsInRect(tempEvtMappedX, tempEvtMappedY, selectionRects_detection(:, i)'), selectInside(i) );
                        end
                        selectInside(i) = max( IsInRect(tempBufferMeanXpos_selection, tempBufferMeanYpos_selection,...
                            selectionRects_detection(:, i)'), selectInside(i) );  
                        if if_eyeResponse == 1
                            %selectInside(i) = max( IsInRect(x_avg, y_avg,...
                            %    selectionRects_detection(:, i)'), selectInside(i) ); 
                            raw_selectInside(i) = IsInRect(x_avg, y_avg, selectionRects_detection(:, i)');         
                            if raw_selectInside(i) == 1 && last_selectInside(i) == 1

                            else
                                t0_selectHolding = tic;
                            end              
                            if raw_selectInside(i) == 1         
                                if toc(t0_selectHolding) > eyeResponseTime/1000
                                    selectInside(i) = 1;
                                end
                            end
                            last_selectInside(i) = raw_selectInside(i);                               
                        end                        
                    end
                elseif isGreenRectLeft(trial_count) == 0
                    if i == 1
                        if if_mouse_mannual_eyeResponse == 0
                            selectInside(i) = IsInRect(mx, my, selectionRects_detection(:, i)');
                        end
                        if currentKeycode == tempKeycode
                            selectInside(i) = max( IsInRect(tempEvtMappedX, tempEvtMappedY, selectionRects_detection(:, i)'), selectInside(i) );
                        end
                        selectInside(i) = max( IsInRect(tempBufferMeanXpos_selection, tempBufferMeanYpos_selection,...
                            selectionRects_detection(:, i)'), selectInside(i) );    
                        if if_eyeResponse == 1
                            %selectInside(i) = max( IsInRect(x_avg, y_avg,...
                            %    selectionRects_detection(:, i)'), selectInside(i) ); 
                            raw_selectInside(i) = IsInRect(x_avg, y_avg, selectionRects_detection(:, i)');         
                            if raw_selectInside(i) == 1 && last_selectInside(i) == 1

                            else
                                t0_selectHolding = tic;
                            end              
                            if raw_selectInside(i) == 1         
                                if toc(t0_selectHolding) > eyeResponseTime/1000
                                    selectInside(i) = 1;
                                end
                            end
                            last_selectInside(i) = raw_selectInside(i);                               
                        end                        
                    end
                end
            elseif choiceMode == 2
                if if_mouse_mannual_eyeResponse == 0
                    selectInside(i) = IsInRect(mx, my, selectionRects_detection(:, i)');
                end
                if currentKeycode == tempKeycode
                    selectInside(i) = max( IsInRect(tempEvtMappedX, tempEvtMappedY, selectionRects_detection(:, i)'), selectInside(i) );
                end
                selectInside(i) = max( IsInRect(tempBufferMeanXpos_selection, tempBufferMeanYpos_selection,...
                    selectionRects_detection(:, i)'), selectInside(i) );
                if if_eyeResponse == 1
                    %selectInside(i) = max( IsInRect(x_avg, y_avg,...
                    %    selectionRects_detection(:, i)'), selectInside(i) );     
                    raw_selectInside(i) = IsInRect(x_avg, y_avg, selectionRects_detection(:, i)');         
                    if raw_selectInside(i) == 1 && last_selectInside(i) == 1

                    else
                        if i == 1
                            t0_selectHolding_choice_i1= tic;
                        elseif i == 2
                            t0_selectHolding_choice_i2= tic;
                        end
                    end              
                    if raw_selectInside(i) == 1   
                        if i == 1
                            if toc(t0_selectHolding_choice_i1) > eyeResponseTime/1000
                                selectInside(i) = 1;
                            end
                        elseif i == 2
                            %if toc(t0_selectHolding_choice_i2) > eyeResponseTime/1000
                            if toc(t0_selectHolding_choice_i2) > eyeResponseTime/1000 && if_disableRightSelection == 0                            
                                selectInside(i) = 1;
                            end                            
                        end
                        %if toc(t0_selectHolding) > eyeResponseTime/1000
                        %    selectInside(i) = 1;
                        %end
                    end
                    last_selectInside(i) = raw_selectInside(i);                       
                end                
            end                                                                                               
                            
            
            if selectInside(i) == 1    
                %Clear mytouchQueueBuffer_selection once selectInside(i) == 1
                mytouchQueueBuffer_selection = cell(1, touchBufferSize);

                trial_para.selecting_RT(trial_count) = toc(t0_isReleaseForTouch);
                trial_para.selecting_point(trial_count) = i;
                if isGreenRectLeft(trial_count) == 1
                    selecting_point = i;
                else
                    selecting_point = 3-i;
                end                
                isSelecting = 0;
                
                newTouchEnableFlag = 1;
                tempKeycode = 0;
                
                t0_isTwoHolding = tic;

                if selecting_point == 1
                    ifSelectOffloading = 0;
                    trial_para.ifSelectOffloading(trial_count) = -1;
                    taskMode_noChoice_Green0_Red1 = 0;
                    
                    if if_R_mustErrorStop == 1 && isWaitGocueBreak(trial_count) == 1
                        ifFreeTouch(trial_count) = ifFreeTouch_raw(trial_count);
                    end  
                    
                    %isTouching = 1;
                    isTwoHolding = 1;
                    isStartTwoHolding = 0;
                    
                    
                    %to correct isFill as taskMode_noChoice_Green0_Red1
                    %changed from 1 to 0, while selecting_point == 2 needn't do it,
                    %because "%%%%%OffloadBeginEnd%%%%%%" in isTouching
                    %will solve this problem
                    if last_ifSelectOffloading == 1
                        isFill = ~isFill;
                        isFill(length(currentFrames)+1:numSquares) = 0;
                        isFill(currentMiddleSeq) = 0;
                    end
                    
                elseif selecting_point == 2
                    ifSelectOffloading = 1;
                    trial_para.ifSelectOffloading(trial_count) = 1;
                    taskMode_noChoice_Green0_Red1 = 1;
                    
                    if if_R_mustErrorStop == 1
                        ifFreeTouch(trial_count) = 0;
                    end
                    
                    %isTouching = 1;
                    isTwoHolding = 1;
                    isStartTwoHolding = 0;
                    
                    %isTwoDelay = 0;
%                     isTouching = 1;
                    lastTouchKeycode = 0;
                    ifTwoHoldFlag = 1;
                    isTouchingDone = 0;
                    %isHoldTouchingDone = 0;
                    isTouching_tempFlag1 = 0;
                    isTouching_tempFlag2 = 0;
                    correctLock_flag = 0;
                    last_freeTouch_toc = 0;
                    tempFlag_toc = 0;
                    isWaitGocueBreak(trial_count) = -1;
                    
                    pause(100/1000)
                    if if_eyeResponse == 0
                        while TouchEventAvail(dev) > 0
                            TouchEventGet(dev, window);
                        end
                    end
                    
                    %tic
                    
                    
                    
                    %stimulusExtendTime                   
                    %tempStimulusExtendTime
                    if ifTurnOff == 1
                        tempStimET = encourageStimET;
                    end
                    %SetMouse(0, 0, window);
                end
                
                m_marker.Data = uint8(2);% Inner trial markers
                trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms
                % Eyelink('Message', 'SELECTING ON');
                %message = sprintf('SELECTING ON');
                %message = sprintf('DECISIONMAKING');
                message = sprintf('SELECTING AND DELAY 2 ON');
                Eyelink('Message', message);
                trial_para.marker.count = trial_para.marker.count + 1;
                m_markerCount.Data = int32(trial_para.marker.count);
                trial_para.marker.content = [trial_para.marker.content, {message}];
                
                
                home % To scrolling down in command window
                fprintf("------------------------------------------------------\n");
                fprintf("currentSequence = ");
                for tempi=1:length(currentSequence)
                    fprintf("%d ", currentSequence(tempi));
                end
                fprintf("\n");
            
                fprintf("ifSelectOffloading = %d!\n", ifSelectOffloading);
                last_ifSelectOffloading = ifSelectOffloading;
                last_selectInside = selectInside;
            end
            
        end
        
        
    %---------------------------------------------------------------------------------------------        
    elseif isTwoHolding == 1
        %% isTwoHolding
        if if_eyeResponse == 1                 
            isReleaseToHold = 1;
        end
        
        if isReleaseToHold == 1
            %isStartTwoHolding = 1;
            isTwoHolding = 0;
            isTwoDelay = 1;
            t0_isTwoDelay = tic;
            
            %m_marker.Data = uint8(2);% Inner trial markers
            % % Eyelink('Message', 'DELAY 2 ON');
            %message = sprintf('DELAY 2 ON');
            %Eyelink('Message', message);
            %trial_para.marker.count = trial_para.marker.count + 1;
            %m_markerCount.Data = int32(trial_para.marker.count);
            %trial_para.marker.content = [trial_para.marker.content, {message}];
            
            continue
        end
                
        temp_t = 0;%100
        if toc(t0_isTwoHolding) < temp_t/1000
            temp_W = toc(t0_isTwoHolding) / (temp_t/1000);
            normalColor_end = normalColor;
            normalColor_end(4) = 0;
            specialColor_decision_end = specialColor_decision;
            specialColor_decision_end(4) = 0;            
            normalColorX = (1-temp_W) .* normalColor + temp_W .* normalColor_end;
            specialColor_decisionX = (1-temp_W) .* specialColor_decision + temp_W .* specialColor_decision_end;                
            %selectionTriangles
            if isGreenRectLeft(trial_count) == 1
                if choiceMode == 0
                    Screen('FillRect', window, normalColorX, selectionRects(:,1));
                elseif choiceMode == 1
                    Screen('FillPoly', window, specialColor_decisionX, selectionTriangles{2}, 1);
                elseif choiceMode == 2
                    Screen('FillRect', window, normalColorX, selectionRects(:,1));
                    Screen('FillPoly', window, specialColor_decisionX, selectionTriangles{2}, 1);
                end                
            else
                if choiceMode == 0
                    Screen('FillRect', window, normalColorX, selectionRects(:,2));
                elseif choiceMode == 1
                    Screen('FillPoly', window, specialColor_decisionX, selectionTriangles{1}, 1);
                elseif choiceMode == 2
                    Screen('FillPoly', window, specialColor_decisionX, selectionTriangles{1}, 1);
                    Screen('FillRect', window, normalColorX, selectionRects(:,2));
                end                
            end        
        end
        
        Screen('FrameOval', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
        %Screen('FrameOval', window, frameColorMatrices2(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
        
        Screen('FillRect', window, [0.3 0.3 0.3], centeredVerticalCrossRect);
        Screen('FillRect', window, [0.3 0.3 0.3], centeredHorizontalCrossRect);
        vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);
        
               
        
    %---------------------------------------------------------------------------------------------    
    elseif isTwoDelay == 1
        %% isTwoDelay
        if if_eyeResponse == 1                 
            isHold = 1;            
            %if fixDis > fixWinSize_delay2
            if fixOn_flag_delay2 == 0
                isHold = 0;
                fprintf('Delay2 break!\n');
            end      
        end              
        if isHold == 0
            %isCorrect(trial_count) = -1;
            isWaitGocueBreak(trial_count) = 1;
            isItiInitial_break = 1;
            isTwoDelay = 0;
            continue
        end
        
        temp_t = 0;%100
        if toc(t0_isTwoHolding) < temp_t/1000
            temp_W = toc(t0_isTwoHolding) / (temp_t/1000);
            normalColor_end = normalColor;
            normalColor_end(4) = 0;
            specialColor_decision_end = specialColor_decision;
            specialColor_decision_end(4) = 0;            
            normalColorX = (1-temp_W) .* normalColor + temp_W .* normalColor_end;
            specialColor_decisionX = (1-temp_W) .* specialColor_decision + temp_W .* specialColor_decision_end;                
            %selectionTriangles
            if isGreenRectLeft(trial_count) == 1
                if choiceMode == 0
                    Screen('FillRect', window, normalColorX, selectionRects(:,1));
                elseif choiceMode == 1
                    Screen('FillPoly', window, specialColor_decisionX, selectionTriangles{2}, 1);
                elseif choiceMode == 2
                    Screen('FillRect', window, normalColorX, selectionRects(:,1));
                    Screen('FillPoly', window, specialColor_decisionX, selectionTriangles{2}, 1);
                end                
            else
                if choiceMode == 0
                    Screen('FillRect', window, normalColorX, selectionRects(:,2));
                elseif choiceMode == 1
                    Screen('FillPoly', window, specialColor_decisionX, selectionTriangles{1}, 1);
                elseif choiceMode == 2
                    Screen('FillPoly', window, specialColor_decisionX, selectionTriangles{1}, 1);
                    Screen('FillRect', window, normalColorX, selectionRects(:,2));
                end                
            end        
        end
        
        %Screen('FrameOval', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
        Screen('FrameOval', window, frameColorMatrices2(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
        
        Screen('FillRect', window, [1 1 1]*0.75, centeredVerticalCrossRect);
        Screen('FillRect', window, [1 1 1]*0.75, centeredHorizontalCrossRect);
        vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);
 
        if toc(t0_isTwoDelay) > twoFixationTime(trial_count)/1000
            if if_eyeResponse == 0
                while TouchEventAvail(dev) > 0
                    TouchEventGet(dev, window);
                end
            end
 
            fprintf("twoFixationTime = %d ms\n", round(twoFixationTime(trial_count)));
            isTwoDelay = 0;
            isTouching = 1;            
            
            isHoldToRelease_first = 0;            
            isHoldToRelease_notFirst = 0;
            HoldToRelease_notFirst_enable = 0;
            
            lastTouchKeycode = 0;
            ifTwoHoldFlag = 1;
            isTouchingDone = 0;
            %isHoldTouchingDone = 0;
            isTouching_tempFlag1 = 0;
            isTouching_tempFlag2 = 0;
            correctLock_flag = 0;
            last_freeTouch_toc = 0;
            tempFlag_toc = 0;
            isWaitGocueBreak(trial_count) = -1;
            
            newTouchEnableFlag = 1;
            newTouchEnableFlag_sumbit = 1;
            tempKeycode = 0;
            submitKeycode = 0;
            submitEvtMappedX = 0;
            submitEvtMappedY = 0;
            submitInside = 0;
            last_submitInside = 0;
            raw_submitInside = 0;
            flag_first_submitInside = 0;
            
            
            
            m_marker.Data = uint8(2);% Inner trial markers
            trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms
            %Eyelink('Message', 'RESPONSE ON');
            message = sprintf('RESPONSE ON');
            Eyelink('Message', message);
            trial_para.marker.count = trial_para.marker.count + 1;
            m_markerCount.Data = int32(trial_para.marker.count);
            trial_para.marker.content = [trial_para.marker.content, {message}];
            
            t0_isTouching = tic;
            
            last_raw_eye_inside = zeros(1, numSquares);            
            raw_eye_inside = zeros(1, numSquares);
            t0_eyeTouching = uint64(zeros(1, numSquares));            
            
            % dequeues all samples and events from the link. AKA. Clear Buffer!
            %[samplesIn, eventsIn, drained] = Eyelink('GetQueuedData');
            
            isTouching_firstFlag = 0;  
            completedFilled = false(1, numSquares);
            completedFilled_raw = completedFilled;
            lock_flag = 0;
            first_lock_flag = 0;
        end
        
        
        
    %---------------------------------------------------------------------------------------------    
    elseif isTouching == 1
        %% isTouching        
        if isHoldToRelease == 1            
            isHoldToRelease_first = 1;            
        end        
        
        if isHoldToRelease_first == 1
            if isHold == 1
                HoldToRelease_notFirst_enable = 1;                
            end
        end
        
        if HoldToRelease_notFirst_enable == 1 && isHoldToRelease == 1
            isHoldToRelease_notFirst = 1;
        end
        
        if isTouching_firstFlag == 0
           isTouching_firstFlag = 1;
           % dequeues all samples and events from the link. AKA. Clear Buffer!
           %[samplesIn, eventsIn, drained] = Eyelink('GetQueuedData');            
        end
        
        
        
        if toc(t0_isTouching) < (tempStimET - fixationTime(trial_count))/1000 %&& ifTwoHoldFlag == 0
            if ifCircle == 0
                if taskMode_noChoice_Green0_Red1 == 0
                    Screen('FillRect', window, currentColorMatrices(:, currentSequence), showRects(:, currentSequence));
                elseif taskMode_noChoice_Green0_Red1 == 1
                    for i=1:length(currentSequence)
                        Screen('FillPoly', window, currentColorMatrices(:, currentSequence(i)), showTriangles{currentSequence(i)}, 1);
                    end
                end
            elseif ifCircle == 1
                Screen('FillOval', window, currentColorMatrices(:, currentSequence), showRects(:, currentSequence));
            end     
        end
        if ifSelection == 1 && toc(t0_isTouching) < 150/1000 && ifTwoHoldFlag == 0
            %selectionTriangles
            if isGreenRectLeft(trial_count) == 1
                if choiceMode == 0
                    Screen('FillRect', window, normalColor, selectionRects(:,1));
                    Screen('FrameRect', window, [1 1 1], selectionRects(:,1), 4);
                    
                elseif choiceMode == 1
                    Screen('FillPoly', window, specialColor_decision, selectionTriangles{2}, 1);
                    Screen('FramePoly', window, [1 1 1], selectionTriangles{2}, 4);                                        
                elseif choiceMode == 2
                    if ifSelectOffloading == 0
                        Screen('FillRect', window, normalColor, selectionRects(:,1));
                        Screen('FrameRect', window, [1 1 1], selectionRects(:,1), 4);
                    elseif ifSelectOffloading == 1
                        Screen('FillPoly', window, specialColor_decision, selectionTriangles{2}, 1);
                        Screen('FramePoly', window, [1 1 1], selectionTriangles{2}, 4);                        
                    end                                       
                end                
            else
                if choiceMode == 0
                    Screen('FillRect', window, normalColor, selectionRects(:,2));
                    Screen('FrameRect', window, [1 1 1], selectionRects(:,2), 4);
                elseif choiceMode == 1
                    Screen('FillPoly', window, specialColor_decision, selectionTriangles{1}, 1);
                    Screen('FramePoly', window, [1 1 1], selectionTriangles{1}, 4); 
                elseif choiceMode == 2
%                     Screen('FillPoly', window, specialColor_decision, selectionTriangles{1}, 1);
%                     Screen('FillRect', window, normalColor, selectionRects(:,2));                    
                    if ifSelectOffloading == 0
                        Screen('FillRect', window, normalColor, selectionRects(:,2));
                        Screen('FrameRect', window, [1 1 1], selectionRects(:,2), 4);
                    elseif ifSelectOffloading == 1
                        Screen('FillPoly', window, specialColor_decision, selectionTriangles{1}, 1);
                        Screen('FramePoly', window, [1 1 1], selectionTriangles{1}, 4);
                    end                    
                end                
            end
        end
        
        tempColorMatrices = repmat(normalColor_offload',1,numSquares);
        
        
        if isTriangletrialMiddle == 1 || ((isTouching_tempFlag1 == 1 || isTouching_tempFlag2 == 1) && taskMode_noChoice_Green0_Red1 == 1)
            if ifCircle == 0
                Screen('FillRect', window, tempColorMatrices(:, currentSequence), showRects(:, currentSequence));
            elseif ifCircle == 1
                Screen('FillOval', window, tempColorMatrices(:, currentSequence), showRects(:, currentSequence));
            end
        end
        
        %%%%%OffloadBeginEnd%%%%%%
        if taskMode_noChoice_Green0_Red1 == 1 && isTouching_tempFlag1 == 0 && isTouching_tempFlag2 == 0
            sortedCurrentSeq = sort(currentSequence);%let begin=left, end=right in random stimulus condition
            if ifBeginEnd == 1
                currentMiddleSeq = sortedCurrentSeq(2:end-1);
                currentBeginEndSeq = [ sortedCurrentSeq(1) sortedCurrentSeq(end) ];
                if ~isempty(currentMiddleSeq)
                    if ifCircle == 0
                        for i=1:length(currentMiddleSeq)
                            Screen('FillPoly', window, currentColorMatrices(:, currentMiddleSeq(i)), showTriangles{currentMiddleSeq(i)}, 1);
                        end
                    elseif ifCircle == 1
                        Screen('FillOval', window, currentColorMatrices(:, currentMiddleSeq), showRects(:, currentMiddleSeq));
                    end
                end
            elseif ifBeginEnd == 0
                currentMiddleSeq = [];
                currentBeginEndSeq = sortedCurrentSeq;
                if ifTurnOff == 1
                    if ifCircle == 0
                        for i=1:length(currentMiddleSeq)
                            Screen('FillPoly', window, currentColorMatrices(:, currentBeginEndSeq(i)), showTriangles{currentBeginEndSeq(i)}, 1);
                        end
                    elseif ifCircle == 1
                        Screen('FillOval', window, currentColorMatrices(:, currentBeginEndSeq), showRects(:, currentBeginEndSeq));
                    end
                end
            end
            
            
            temp_isFill_OffloadBeginEnd = zeros(1, length(isFill));
            temp_isFill_OffloadBeginEnd(currentBeginEndSeq) = 1;
            %temp_isFill = temp_isFill_OffloadBeginEnd;
            isFill = temp_isFill_OffloadBeginEnd;
        end
        
        
        
        [mx, my, buttons] = GetMouse(window);
        inside = zeros(1, numSquares);
        if if_mouse_mannual_eyeResponse == 0
            for i = 1:numSquares
                if ifCircle == 0
                    inside(i) = IsInRect(mx, my, detectRects(:, i)');
                elseif ifCircle == 1
                    tempDis = sqrt((squareXpos(i)-mx)^2 + (squareYpos(i)-my)^2);
                    if tempDis <= detectBaseRect(4)/2
                        inside(i) = 1;
                        inside(i) = 1;
                    else
                        inside(i) = 0;
                    end
                end
                if inside(i) == 1 && isFilled(i) == 0 && isFrame(i) == 1  

                    last_freeTouch_toc = toc(t0_isTouching);

                    if isTriangletrialMiddle == 0
                        isFilled(i) = 1;                    
                        if taskMode_noChoice_Green0_Red1 == 0                        
                            if ismember(i,newTouchPoint) == 0
                                newTouchPoint = [newTouchPoint i]; %#ok<AGROW>
                            end
                        elseif taskMode_noChoice_Green0_Red1 == 1
                            if ismember(i,currentMiddleSeq) == 0
                                if ismember(i,newTouchPoint) == 0
                                    newTouchPoint = [newTouchPoint i]; %#ok<AGROW>
                                end                        
                            elseif ismember(i,currentMiddleSeq) == 1
                                if ifAllowRedMiddle == 1
                                    isFilled(i) = 0;
                                elseif ifAllowRedMiddle == 0
                                    if ismember(i,newTouchPoint) == 0
                                        newTouchPoint = [newTouchPoint i]; %#ok<AGROW>
                                    end
                                end
                            end
                        end
                    elseif isTriangletrialMiddle == 1
                        if isFill(i) == 1
                            isFilled(i) = 1;
                            if ismember(i,newTouchPoint) == 0
                                newTouchPoint = [newTouchPoint i]; %#ok<AGROW>
                            end
                        end
                    end

                end
            end
        end
        
        if if_mouse_mannual_eyeResponse == 0
            raw_submitInside = IsInRect(mx, my, submitRects_dection');
        end
        
        if if_eyeResponse == 1
            eye_inside = zeros(1, numSquares);
            % for response            home % To scrolling down in command window
            for i = 1:numSquares
                if completedFilled(i) == true
                    continue
                end
                if completedFilled_raw(i) == true
                    continue
                end                
                if ifCircle == 0
                    raw_eye_inside(i) = IsInRect(x_avg, y_avg, detectRects(:, i)');
                elseif ifCircle == 1
                    tempDis = sqrt((squareXpos(i)-x_avg)^2 + (squareYpos(i)-y_avg)^2);
                    if tempDis <= detectBaseRect(4)/2
                        raw_eye_inside(i) = 1;
                    else
                        raw_eye_inside(i) = 0;
                    end
                end
                
                if raw_eye_inside(i) == 1 && last_raw_eye_inside(i) == 1
                    a = 1;
                else
                    t0_eyeTouching(i) = tic;
                end     
                
if raw_eye_inside(i) == 0 && lock_flag == i    
    eyeTouchingLevel = floor(fixupdateTimeFoldNum/2);
    Screen('FillOval', window, frameColor, showRects_fold(eyeTouchingLevel,:, i)');
end
                
                if raw_eye_inside(i) == 1
                    %if toc(t0_eyeTouching(i)) > eyeResponseTime/1000 || ...
                    %        (toc(t0_eyeTouching(i)) > fixupdateTimeBound_first/1000 && sum(completedFilled) == 0)
                    
                    temp_TimeBound = 0;
                    if sum(completedFilled) == 0
                        if taskMode_noChoice_Green0_Red1 == 0
                            temp_TimeBound = fixupdateTimeBound_first/1000;
                        elseif taskMode_noChoice_Green0_Red1 == 1
                            temp_TimeBound = fixupdateTimeBound/1000;
                        end
                    elseif sum(completedFilled) > 0
                        temp_TimeBound = fixupdateTimeBound/1000;
                    end
                    
                    %if (toc(t0_eyeTouching(i)) > fixupdateTimeBound_first/1000 && sum(completedFilled) == 0) || ...
                    %        (toc(t0_eyeTouching(i)) > fixupdateTimeBound/1000 && sum(completedFilled) > 0)
                    %if toc(t0_eyeTouching(i)) > temp_TimeBound
                    if toc(t0_eyeTouching(i)) > temp_TimeBound && lock_flag == i
                        if isFilled(i) == 0 && eye_inside(i)  == 0
                            eye_inside(i) = 1;
                            fprintf('Fixation time (%d)=%.1f ms.',i,toc(t0_eyeTouching(i))*1000); 
                            
                            Screen('FrameOval', window, frameColor3, showRects_fold(eyeTouchingLevel,:, i)', frameWidth*2, frameWidth*2);
                            
                            lock_flag = 0;
                        end
                    else
                        
                        if toc(t0_eyeTouching(i)) > fixupdateTimeBound_lock/1000 && lock_flag == 0
                            lock_flag = i;
                            fprintf('Locked at item %d.',i);
                            
                            first_lock_flag = 1;
                        end
                        
                        if lock_flag == 0 || lock_flag == i
                            eyeTouchingLevel = find(toc(t0_eyeTouching(i))>fixupdateTimeBound_fold/1000,1,'last');
                            eyeTouchingLevel = min(eyeTouchingLevel, fixupdateTimeFoldNum);

                            if isempty(eyeTouchingLevel) == 0
                                Screen('FillOval', window, frameColor, showRects_fold(eyeTouchingLevel,:, i)');
                            end                                                        
                        end

                    end  

%                     if if_mouse_mannual_eyeResponse == 0 && sum(completedFilled) > 0
%                         if isFilled(i) == 0 && eye_inside(i)  == 0
%                             
%                             valid_tempi = 0;
%                             while valid_tempi == 0
%                                 [samplesIn, eventsIn, drained] = Eyelink('GetQueuedData');
%                                 if isempty(eventsIn)
%                                     continue
%                                 end
%                                 queuedEvents_type = eventsIn(2,:);
%                                 
%                                 valid_tempi = 0;
%                                 temp_validBoolIndex = (...
%                                     queuedEvents_type == el.FIXUPDATE ...
%                                     | queuedEvents_type == el.ENDFIX ...
%                                     | queuedEvents_type == el.ENDSACC ...
%                                     );
%                                 for tempi=length(temp_validBoolIndex):-1:1
%                                     if temp_validBoolIndex(tempi) == true
%                                         valid_tempi = tempi;
%                                         valid_type = queuedEvents_type(tempi);
%                                         break
%                                     end
%                                 end
%                                 
%                             end
%                             valid_tempi;
%                             valid_type;
%                             
%                             %fixupdateTimeBound = 205; %150-->175-->205
%                             if toc(t0_eyeTouching(i))*1000 > fixupdateTimeBound
%                                 temp_fixupdateFlag = 1;
%                             else
%                                 temp_fixupdateFlag = 0;
%                             end
% 
%                             
%                             if valid_type == el.FIXUPDATE || valid_type == el.ENDFIX
%                                 temp_x_avg = eventsIn(19,valid_tempi) + reCenterGazeBias_x;%display gaze position average x
%                                 temp_y_avg = eventsIn(20,valid_tempi) + reCenterGazeBias_y;%display gaze position average y
%                             elseif valid_type == el.ENDSACC
%                                 temp_x_avg = eventsIn(14,valid_tempi) + reCenterGazeBias_x;%display gaze position ending point x
%                                 temp_y_avg = eventsIn(15,valid_tempi) + reCenterGazeBias_y;%display gaze position ending point y
%                             end
%                             
%                             if temp_fixupdateFlag == 1
%                                 temp_Dis = sqrt((squareXpos(i)-temp_x_avg)^2 + (squareYpos(i)-temp_y_avg)^2);
%                                 if temp_Dis <= detectBaseRect(4)/2
%                                     eye_inside(i) = 1;
%                                     fprintf('Fixation time (%d)=%.1f ms, ',i,toc(t0_eyeTouching(i))*1000);
%                                     fprintf('eye_evtype=%d.',valid_type);
%                                 end
%                             end
%                             
%                             
%                         end
%                     end
                    
                end
                last_raw_eye_inside(i) = raw_eye_inside(i);                
                
                if eye_inside(i) == 1 && isFilled(i) == 0 && isFrame(i) == 1  

                    last_freeTouch_toc = toc(t0_isTouching);

                    if sum(completedFilled) == 0
                        %home % To scrolling down in command window
                    end
                    
if eye_inside(i) == 1 && sum(isFilled(i+1:end)) > 0 && isFrame(i) == 1 && ifNoBackTouch == 1    
    isCorrect(trial_count) = -1;
    isTouchFrameCorrect(trial_count) = -1;
    isBackTouch(trial_count) = 1;
    isTouchingDone = 1;
    fprintf("BackTouch error!\n");
end   

                    
                    if isTriangletrialMiddle == 0
                        isFilled(i) = 1;
                        if taskMode_noChoice_Green0_Red1 == 0
                            if ismember(i,newTouchPoint) == 0
                                
old_data = m_marker.Data;
m_marker.Data = uint8(123);
% WaitSecs(0.1);
% giveWater(20);
giveWater(saccadeWater,waterDropNum,0);
m_marker.Data = uint8(old_data);  

% WaitSecs(saccadeWater/1000);

                                newTouchPoint = [newTouchPoint i]; %#ok<AGROW>
                                fprintf("Touch greed point %d !\n", i);
                                completedFilled(i) = true;

                                m_marker.Data = uint8(2);% Inner trial markers
                                trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms
                                %Eyelink('Message', 'GREEN RESPONSE');
                                %message = sprintf('GREEN RESPONSE');
                                message = sprintf('GREEN RESPONSE %d', i);
                                Eyelink('Message', message);
                                trial_para.marker.count = trial_para.marker.count + 1;
                                m_markerCount.Data = int32(trial_para.marker.count);
                                trial_para.marker.content = [trial_para.marker.content, {message}];                                    
                            end
                        elseif taskMode_noChoice_Green0_Red1 == 1
                            if ismember(i,currentMiddleSeq) == 0
                                temp_touch_flag = 1;
                                if if_R_mustCorrect == 1
                                    if ismember(i,newTouchPoint) == 0
                                        if currentSequence(length(newTouchPoint)+1) ~= i
                                            isFilled(i) = 0;
                                            temp_touch_flag = 0;                                            
                                            fprintf('R_mustCorrect!\n');
                                            if if_eyeResponse == 1
                                                t0_eyeTouching(i) = tic;
                                            end
                                        end
                                    end
                                end
                                
                                %if ismember(i,newTouchPoint) == 0
                                if ismember(i,newTouchPoint) == 0 && temp_touch_flag == 1
                                    
                                
% old_data = m_marker.Data;
% m_marker.Data = uint8(123);
% %WaitSecs(0.1);
% % giveWater(20);
% giveWater(saccadeWater,waterDropNum,0);
% m_marker.Data = uint8(old_data); 

WaitSecs(saccadeWater/1000);
                                    
                                    newTouchPoint = [newTouchPoint i]; %#ok<AGROW>
                                    fprintf("Touch red point %d !\n", i);
                                    completedFilled(i) = true;

                                    m_marker.Data = uint8(2);% Inner trial markers
                                    trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms
                                    %Eyelink('Message', 'RED RESPONSE');
                                    %message = sprintf('RED RESPONSE');
                                    message = sprintf('RED RESPONSE %d', i);
                                    Eyelink('Message', message);
                                    trial_para.marker.count = trial_para.marker.count + 1;
                                    m_markerCount.Data = int32(trial_para.marker.count);
                                    trial_para.marker.content = [trial_para.marker.content, {message}];                                                                            
                                end
                            elseif ismember(i,currentMiddleSeq) == 1
                                if ifAllowRedMiddle == 1
                                    isFilled(i) = 0;
                                elseif ifAllowRedMiddle == 0
                                    if ismember(i,newTouchPoint) == 0
                                        newTouchPoint = [newTouchPoint i]; %#ok<AGROW>
                                    end
                                end
                            end
                        end
                    elseif isTriangletrialMiddle == 1
                        if isFill(i) == 1
                            isFilled(i) = 1;
                            if ismember(i,newTouchPoint) == 0
                                
old_data = m_marker.Data;
m_marker.Data = uint8(123);
%WaitSecs(0.1);
% giveWater(20);
giveWater(saccadeWater,waterDropNum,0);
m_marker.Data = uint8(old_data); 

% WaitSecs(saccadeWater/1000);

                                newTouchPoint = [newTouchPoint i]; %#ok<AGROW>
                                fprintf("Touch greed point %d !\n", i);
                                completedFilled(i) = true;

                                m_marker.Data = uint8(2);% Inner trial markers
                                trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms
                                %Eyelink('Message', 'GREEN RESPONSE');
                                %message = sprintf('GREEN RESPONSE');
                                message = sprintf('GREEN RESPONSE %d', i);
                                Eyelink('Message', message);
                                trial_para.marker.count = trial_para.marker.count + 1;
                                m_markerCount.Data = int32(trial_para.marker.count);
                                trial_para.marker.content = [trial_para.marker.content, {message}];                                                                        
                            end
                        end
                    end
                end
            end
            % response end

            % for submit                        
            raw_submitInside = max( IsInRect(x_avg, y_avg, submitRects_dection'), raw_submitInside );            
            
            if first_lock_flag == 0
                raw_submitInside = 0;
            end
            
            if lock_flag ~= 0
                raw_submitInside = 0;
            end
            
            if raw_submitInside == 1 && last_submitInside == 1

            else
               t0_submitHolding = tic;
            %    %fprintf("submitHolding start! \n");
            end              
            %if raw_submitInside == 1 && last_submitInside == 0
            %    t0_submitHolding = tic;
            %    fprintf("submitHolding start! \n");
            %end
            if raw_submitInside == 1         
                if toc(t0_submitHolding) > submitHoldingTime/1000
                    submitInside = 1;
                end
            end
            last_submitInside = raw_submitInside;
            raw_submitInside = 0;
            % submit end            
        end
            
        %newTouchPoint = [];
        %tempi = 0;
        if if_eyeResponse == 0
            while TouchEventAvail(dev) > 0
                evt = TouchEventGet(dev, window);
                %evt.Type == 2 means new touch point
                %evt.Type == 3 means moving touch point
                %evt.Type == 4 means touch finish
                
                % for submit
                %if evt.Type == 2 || evt.Type == 3 || evt.Type == 4
                %if evt.Type == 2 || (flag_first_submitInside == 1 && evt.Type == 3 && submitKeycode == evt.Keycode)
                if evt.Type == 2 || (flag_first_submitInside == 1 && evt.Type == 3)
                    submitEvtMappedX = evt.MappedX;
                    submitEvtMappedY = evt.MappedY;
                end
                
                %raw_submitInside = IsInRect(mx, my, submitRects_dection');
                raw_submitInside = max( IsInRect(submitEvtMappedX, submitEvtMappedY, submitRects_dection'), raw_submitInside );
                if raw_submitInside == 1 && last_submitInside == 0
                    t0_submitHolding = tic;
                    fprintf("submitHolding start! \n");
                end
                if raw_submitInside == 1
                    if evt.Type == 2
                        submitKeycode = evt.Keycode;
                    end
                    flag_first_submitInside = 1;
                    if toc(t0_submitHolding) > submitHoldingTime/1000
                        submitInside = 1;
                    end
                end
                last_submitInside = raw_submitInside;
                
                submitEvtMappedX = 0;
                submitEvtMappedY = 0;
                % submit end
                
                if evt.Type == 4
                    if evt.Keycode == tempKeycode
                        newTouchEnableFlag = 1;
                    end
                end
                %             if evt.Type == 2
                if evt.Type == 2 || evt.Type == 3
                    if evt.Type == 2
                        tempFlag_toc = 1;
                    end
                    if evt.Type == 3
                        if evt.Keycode == tempKeycode
                            newTouchEnableFlag = 1;
                        end
                    end
                    
                    if flag_first_submitInside == 1
                        %newTouchEnableFlag = 0;
                    end
                    
                    if isempty(mytouchQueueBuffer{1}) == 1 || touchBufferSize == 1
                        mytouchQueueBuffer{1} = [evt.MappedX evt.MappedY];
                    elseif isempty(mytouchQueueBuffer{1}) == 0
                        tempDis = norm([evt.MappedX evt.MappedY] - mytouchQueueBuffer{1});
                        if tempDis > detectBaseRect(4)*2
                            mytouchQueueBuffer = cell(1, touchBufferSize);
                            mytouchQueueBuffer{1} = [evt.MappedX evt.MappedY];
                        else
                            if isempty(mytouchQueueBuffer{touchBufferSize}) == 1
                                for i=2:touchBufferSize
                                    if isempty(mytouchQueueBuffer{i}) == 1
                                        for j=i:-1:2
                                            mytouchQueueBuffer{j} = mytouchQueueBuffer{j-1};
                                        end
                                        mytouchQueueBuffer{1} = [evt.MappedX evt.MappedY];
                                        break
                                    end
                                end
                            elseif isempty(mytouchQueueBuffer{touchBufferSize}) == 0
                                for j=touchBufferSize:-1:2
                                    mytouchQueueBuffer{j} = mytouchQueueBuffer{j-1};
                                end
                                mytouchQueueBuffer{1} = [evt.MappedX evt.MappedY];
                            end
                        end
                    end
                    tempTouchBufferLen = 0;
                    for i=1:touchBufferSize
                        if isempty(mytouchQueueBuffer{i}) == 0
                            tempTouchBufferLen = tempTouchBufferLen + 1;
                        end
                    end
                    
                    %allow radiant detect
                    if enableRadiantDetect == 1
                        distanceToOrigin = norm([evt.MappedX evt.MappedY] - [originalPointXpos originalPointYpos]);
                        tempCondition1 = distanceToOrigin>=radiantDetectR1 && distanceToOrigin<=radius-max(detectBaseRect)/2;
                        tempCondition2 = distanceToOrigin>=radius+max(detectBaseRect)/2 && distanceToOrigin<=radiantDetectR2;
                        tempCondition = tempCondition1 == 1 || tempCondition2 == 1;
                        if tempCondition == 1
                            tempDis = 0;
                            minDis = 2000;
                            minFlag = 1;
                            for i = 1:numSquares
                                tempDis = norm([evt.MappedX evt.MappedY] - [squareXpos(i) squareYpos(i)]);
                                if tempDis < minDis
                                    minDis = tempDis;
                                    minFlag = i;
                                end
                            end
                            inside(minFlag) = 1;
                        end
                    else
                        tempCondition = 0;
                    end
                    
                    for i = 1:numSquares
                        if tempCondition == 0
                            if ifCircle == 0  && newTouchEnableFlag == 1
                                %inside(i) = IsInRect(mx, my, detectRects(:, i)');
                                inside(i) = max( IsInRect(evt.MappedX, evt.MappedY, detectRects(:, i)'), inside(i) );
                            elseif ifCircle == 1
                                %tempDis = sqrt((squareXpos(i)-evt.MappedX)^2 + (squareYpos(i)-evt.MappedY)^2);
                                tempDis = norm([evt.MappedX evt.MappedY] - [squareXpos(i) squareYpos(i)]);
                                if tempDis <= detectBaseRect(4)/2 && newTouchEnableFlag == 1
                                    inside(i) = 1;
                                end
                                if isempty(mytouchQueueBuffer{1}) == 0
                                    tempBufferMeanXpos = 0;
                                    tempBufferMeanYpos = 0;
                                    for tempj=1:tempTouchBufferLen
                                        tempBufferMeanXpos = tempBufferMeanXpos + mytouchQueueBuffer{tempj}(1);
                                        tempBufferMeanYpos = tempBufferMeanYpos + mytouchQueueBuffer{tempj}(2);
                                    end
                                    tempBufferMeanXpos = tempBufferMeanXpos/tempTouchBufferLen;
                                    tempBufferMeanYpos = tempBufferMeanYpos/tempTouchBufferLen;
                                    
                                    if touchBuffer_enable == 0
                                        tempBufferMeanXpos = 0;
                                        tempBufferMeanYpos = 0;
                                    end
                                    
                                    tempDis2 = norm([tempBufferMeanXpos tempBufferMeanYpos] - [squareXpos(i) squareYpos(i)]);
                                    if tempDis2 <= detectBaseRect(4)/2
                                        inside(i) = 1;
                                    end
                                end
                                
                            end
                        end
                        
                        
                        if inside(i) == 1 && sum(isFilled(i+1:end)) > 0 && isFrame(i) == 1 && ifNoBackTouch == 1
                            
                            isCorrect(trial_count) = -1;
                            isTouchFrameCorrect(trial_count) = -1;
                            isBackTouch(trial_count) = 1;
                            isTouchingDone = 1;
                            fprintf("BackTouch error!\n");
                        end
                        
                        if inside(i) == 1 && isFilled(i) == 1 && isFrame(i) == 1
                            tempFlag_toc = 0;
                        end
                        
                        if inside(i) == 1 && isFilled(i) == 0 && isFrame(i) == 1 && evt.Keycode ~= lastTouchKeycode
                            
                            %submitEvtMappedX = 0;
                            %submitEvtMappedY = 0;
                            
                            lastTouchKeycode = evt.Keycode;
                            
                            last_freeTouch_toc = toc(t0_isTouching);
                            
                            %Clear mytouchQueueBuffer once inside(i) == 1
                            mytouchQueueBuffer = cell(1, touchBufferSize);
                            
                            if isTriangletrialMiddle == 0
                                isFilled(i) = 1;
                                if taskMode_noChoice_Green0_Red1 == 0
                                    if ismember(i,newTouchPoint) == 0
                                        newTouchPoint = [newTouchPoint i]; %#ok<AGROW>
                                        fprintf("Touch greed point %d !\n", i);
                                        
                                        m_marker.Data = uint8(2);% Inner trial markers
                                        trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms
                                        %Eyelink('Message', 'GREEN RESPONSE');
                                        %message = sprintf('GREEN RESPONSE');
                                        message = sprintf('GREEN RESPONSE %d', i);
                                        Eyelink('Message', message);
                                        trial_para.marker.count = trial_para.marker.count + 1;
                                        m_markerCount.Data = int32(trial_para.marker.count);
                                        trial_para.marker.content = [trial_para.marker.content, {message}];
                                    end
                                elseif taskMode_noChoice_Green0_Red1 == 1
                                    if ismember(i,currentMiddleSeq) == 0
                                        if ismember(i,newTouchPoint) == 0
                                            newTouchPoint = [newTouchPoint i]; %#ok<AGROW>
                                            fprintf("Touch red point %d !\n", i);
                                            
                                            m_marker.Data = uint8(2);% Inner trial markers
                                            trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms
                                            %Eyelink('Message', 'RED RESPONSE');
                                            %message = sprintf('RED RESPONSE');
                                            message = sprintf('RED RESPONSE %d', i);
                                            Eyelink('Message', message);
                                            trial_para.marker.count = trial_para.marker.count + 1;
                                            m_markerCount.Data = int32(trial_para.marker.count);
                                            trial_para.marker.content = [trial_para.marker.content, {message}];
                                        end
                                    elseif ismember(i,currentMiddleSeq) == 1
                                        if ifAllowRedMiddle == 1
                                            isFilled(i) = 0;
                                        elseif ifAllowRedMiddle == 0
                                            if ismember(i,newTouchPoint) == 0
                                                newTouchPoint = [newTouchPoint i]; %#ok<AGROW>
                                            end
                                        end
                                    end
                                end
                            elseif isTriangletrialMiddle == 1
                                if isFill(i) == 1
                                    isFilled(i) = 1;
                                    if ismember(i,newTouchPoint) == 0
                                        newTouchPoint = [newTouchPoint i]; %#ok<AGROW>
                                        fprintf("Touch greed point %d !\n", i);
                                        
                                        m_marker.Data = uint8(2);% Inner trial markers
                                        trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms
                                        %Eyelink('Message', 'GREEN RESPONSE');
                                        message = sprintf('GREEN RESPONSE');
                                        Eyelink('Message', message);
                                        trial_para.marker.count = trial_para.marker.count + 1;
                                        m_markerCount.Data = int32(trial_para.marker.count);
                                        trial_para.marker.content = [trial_para.marker.content, {message}];
                                    end
                                end
                            end
                            
                        end
                    end
                    if newTouchEnableFlag == 1
                        % Draw a white dot where the touch point is
                        % Screen('DrawDots', window, [evt.MappedX evt.MappedY], 50, white, [], 2);
                        newTouchEnableFlag = 0;
                        tempKeycode = evt.Keycode;
                    end
                end
            end
        end
        
%         if if_R_mustCorrect == 1
%             if taskMode_noChoice_Green0_Red1 == 1
%                 for i = 1:numSquares
%                     if isFill(i) == 0
%                         isFilled(i) = 0;
%                         if ismember(i,newTouchPoint) == 1
%                             newTouchPoint(find(newTouchPoint==i)) = []; %#ok<SAGROW>
%                         end
%                     end
%                 end
%                 for tempi=1:sum(isFill)
%                     tempj = find(isFill==1);
%                     if tempi >= 2
%                         if isFilled(tempj(tempi)) == 1 && isFilled(tempj(tempi-1)) == 0
%                             isFilled(tempj(tempi)) = 0;
%                             newTouchPoint(find(newTouchPoint==tempj(tempi))) = []; %#ok<SAGROW>
%                         end
%                     end
%                 end
%             end
%         end

        
        if  ifNoBackTouch == 0
            if length(newTouchPoint) > 1
                if newTouchPoint(end) < newTouchPoint(end-1)
                    isFilled(newTouchPoint(end)) = 0;
                    newTouchPoint = newTouchPoint(1:end-1);
                end
            end
        end
        
        if ifCorrectLock == 1
            if sum(isFill == isFilled) == length(isFill) && taskMode_noChoice_Green0_Red1 == 0
                correctLock_flag = 1;                
            end
            if correctLock_flag == 1
                
                if ifCorrectLock_last_toErrorStop == 1 %&& numSquares == 8
                    isFilled = [isFill(1:numSquares-1) isFilled(numSquares)];
                else
                    isFilled = isFill;
                end
                
                if isempty(newTouchPoint) == 0
                    newTouchPoint = newTouchPoint(1:sum(isFilled));
                end
            end
        end
        
        for i = 1:numSquares
            if isFilled(i) == 1 %&& isTouching_tempFlag1 == 0
                if ismember(i,newTouchPoint) == 1
                    %Error report
                    if sum(isFilled(1:i)==isFill(1:i)) ~= i
                        if ifExplore == 0
                            if ifCircle == 0
                                if taskMode_noChoice_Green0_Red1 == 0
                                    Screen('FillRect', window, errorColor, showRects(:, i)');
                                elseif taskMode_noChoice_Green0_Red1 == 1
                                    Screen('FillPoly', window, errorColor, showTriangles{i}, 1);
                                end
                            elseif ifCircle == 1
                                if ifFreeTouch(trial_count) == 0
                                    %isTouchingDone = 1;
                                    Screen('FillOval', window, errorColor, showRects(:, i)');
                                elseif ifFreeTouch(trial_count) == 1
                                    if taskMode_noChoice_Green0_Red1 == 0
                                        Screen('FillOval', window, normalColor, showRects(:, i)');
                                    elseif taskMode_noChoice_Green0_Red1 == 1
                                        Screen('FillOval', window, errorColor, showRects(:, i)');
                                        %isTouchingDone = 1;
                                    end
                                end
                            end
                            isCorrect(trial_count) = -1;
                            isTouchFrameCorrect(trial_count) = -1;
                            %isTouching = 0;
                            %isItiInitial = 1;
                            if ifFreeTouch(trial_count) == 0
                                isTouchingDone = 1;
                                %Screen('FillOval', window, errorColor, showRects(:, i)');
                            elseif ifFreeTouch(trial_count) == 1
                                if taskMode_noChoice_Green0_Red1 == 0
                                    %Screen('FillOval', window, normalColor, showRects(:, i)');
                                    
                                    sortedCurrentSeq = sort(currentSequence);                                    
                                    if toc(t0_isTouching) < (tempStimET - fixationTime(trial_count))/1000
                                        if ismember(i, sortedCurrentSeq)
                                            isTouchingDone = 1;
                                            Screen('FillOval', window, errorColor, showRects(:, i)');
                                        end
                                    end
                                    if ifCorrectLock_last_toErrorStop == 1 %&& numSquares == 8
                                        if ismember(i, sortedCurrentSeq) == 1 && i == numFrames_rangeTail
                                            isTouchingDone = 1;
                                            Screen('FillOval', window, errorColor, showRects(:, i)');
                                        end
                                    end
                                    
                                elseif taskMode_noChoice_Green0_Red1 == 1
                                    %Screen('FillOval', window, errorColor, showRects(:, i)');
                                    isTouchingDone = 1;
                                end                                
                            end
                           
                        elseif ifExplore == 1
                            fprintf('randExplore(trial_count)=%.2f, propExplore=%.1f\n',...
                                randExplore(trial_count),propExplore);
                            if randExplore(trial_count) <= propExplore
                                isFilled(i) = 0;
                            else
                                if ifCircle == 0
                                    if taskMode_noChoice_Green0_Red1 == 0
                                        Screen('FillRect', window, errorColor, showRects(:, i)');
                                    elseif taskMode_noChoice_Green0_Red1 == 1
                                        Screen('FillPoly', window, errorColor, showTriangles{i}, 1);
                                    end
                                elseif ifCircle == 1
                                    Screen('FillOval', window, errorColor, showRects(:, i)');
                                end
                                isCorrect(trial_count) = -1;
                                isTouchFrameCorrect(trial_count) = -1;
                                isTouching = 0;
                                isItiInitial = 1;
                            end
                        end
                    else
                        if ifCircle == 0
                            if taskMode_noChoice_Green0_Red1 == 0
                                Screen('FillRect', window, normalColor, showRects(:, i)');
                            elseif taskMode_noChoice_Green0_Red1 == 1                                
                                Screen('FillPoly', window, normalColor_offload, showTriangles{i}, 1);
                            end
                        elseif ifCircle == 1
                            if taskMode_noChoice_Green0_Red1 == 0 || isTouching_tempFlag1 == 1 || isTouching_tempFlag2 == 1
                                Screen('FillOval', window, normalColor, showRects(:, i)');
                            elseif taskMode_noChoice_Green0_Red1 == 1
                                Screen('FillOval', window, normalColor_offload, showRects(:, i)');
                            end
                        end
                    end
                else
                    if ifCircle == 0
                        if taskMode_noChoice_Green0_Red1 == 0
                            Screen('FillRect', window, normalColor, showRects(:, i)');
                        elseif taskMode_noChoice_Green0_Red1 == 1
                            Screen('FillPoly', window, normalColor_offload, showTriangles{i}, 1);
                        end
                    elseif ifCircle == 1
                        if taskMode_noChoice_Green0_Red1 == 0
                            Screen('FillOval', window, normalColor, showRects(:, i)');
                        elseif taskMode_noChoice_Green0_Red1 == 1
                            Screen('FillOval', window, normalColor_offload, showRects(:, i)');
                        end
                    end
                end
            end
        end
                
        if ifCircle == 0
            if taskMode_noChoice_Green0_Red1 == 0
                if ifStep0 == 0
                    Screen('FrameRect', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), 3);
                end
            elseif taskMode_noChoice_Green0_Red1 == 1
                for i=1:length(currentFrames)
                    Screen('FramePoly', window, frameColorMatrices(:, currentFrames(i)), showTriangles{currentFrames(i)}, 1);
                end
            end
        elseif ifCircle == 1
            Screen('FrameOval', window, frameColorMatrices(:, currentFrames), showRects(:, currentFrames), frameWidth, frameWidth);
        end
        
%         raw_submitInside = IsInRect(mx, my, submitRects_dection');
%         raw_submitInside = max( IsInRect(submitEvtMappedX, submitEvtMappedY, submitRects_dection'), raw_submitInside );
%         
%         if raw_submitInside == 1 && last_submitInside == 0
%             t0_submitHolding = tic;   
%             fprintf("submitHolding start! \n");           
%         end
%         if raw_submitInside == 1
%             flag_first_submitInside = 1;
%             if toc(t0_submitHolding) > submitHoldingTime/1000
%                 submitInside = 1;
%             end
%         end
%         last_submitInside = raw_submitInside;
%         
%         submitEvtMappedX = 0;
%         submitEvtMappedY = 0;
        
        if first_lock_flag == 1
            Screen('FrameRect', window, frameColor, submitRects, 3);
        end
        %Screen('FrameRect', window, frameColor, submitRects_dection, 3);
        
        %if ifDrawDots_eye == 1 && if_eyeResponse == 1
        if if_eyeResponse == 1
            %Screen('DrawDots', window, [x_avg y_avg], 20, [0.5 0.5 0.5], [], 2);
        end        
        % Flip to the screen
        vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);
        if sum(isFilled) == sum(isFill) && isTouching_tempFlag1 == 0 && isTouching_tempFlag2 == 0
            %Correct report
            if sum(isFill == isFilled) == length(isFill)  
                if isCorrect(trial_count) == 0
                    isCorrect(trial_count) = 1;
                    isTouchFrameCorrect(trial_count) = 1;
                end
                
                if ifFreeTouch(trial_count) == 0 && taskMode_noChoice_Green0_Red1 == 0
                %if taskMode_noChoice_Green0_Red1 == 0
                    isTouchingDone = 1;
                end
            end
            if ifLinkTriangleSquare_inOneTrial == 1
                if taskMode_noChoice_Green0_Red1 == 0
                    
                    %pause(freeTouch_RT_limit/1000);                    
                    if ifFreeTouch(trial_count) == 0
                        isTouchingDone = 1;
                        isTouching_tempFlag2 = 1;
                    elseif ifFreeTouch(trial_count) == 1
                        isTouching_tempFlag1 = 1;
                    end
                    %isTouching = 0;
                    %isItiInitial = 1;                    
                    if isTriangletrialMiddle == 1
                        isTriangletrialMiddle = 0;
                        taskMode_noChoice_Green0_Red1 = 1;
                    end
                elseif taskMode_noChoice_Green0_Red1 == 1
                    isFill = ~isFill;
                    isFill(length(currentFrames)+1:numSquares) = 0;
                    isFill(currentMiddleSeq) = 0;
                    isFilled = zeros(1, numSquares);
                    newTouchPoint = [];
                    isTriangletrialMiddle = 1;
                    
                    completedFilled_raw = completedFilled;
                    completedFilled = false(1, numSquares);
                    
%                     isTouching = 0;
                    
                    
                    
                    if isCorrect(trial_count) ~= -1
                        isCorrect(trial_count) = 0;
                        isTouchFrameCorrect(trial_count) = 0;
                    end
                    
                    if isCorrect(trial_count) == -1
                        isTouchingDone = 1;
                    else
%                         isTwoHolding = 1;
%                         pause(200/1000);
                    end                    
                    
                    %SetMouse(0, 0, window);
                    %pause(0.5);
%                     while TouchEventAvail(dev) > 0
%                         TouchEventGet(dev, window);
%                     end
                    taskMode_noChoice_Green0_Red1 = 0;                                      
                    
                end
            elseif ifLinkTriangleSquare_inOneTrial == 0
                isTouching = 0;
                isItiInitial = 1;
            end
        end
        
        
        %if (toc(t0_isTouching)-last_freeTouch_toc) > freeTouch_RT_limit/1000 && ifFreeTouch(trial_count) == 1 && sum(isFilled) >= 1
        if (toc(t0_isTouching)-last_freeTouch_toc) > freeTouch_RT_limit/1000 && ifFreeTouch(trial_count) == 1
            isTouchingDone = 1;                                         
        end       
        
        
        if submitInside == 1 %isReleaseToHold, isHoldToRelease_notFirst
            isTouchingDone = 1;
        end
            
        if submitInside == 0 && isCorrect(trial_count) == 1
            isTouchingDone = 0;
        end
        
        if if_validHoldUntilCorrect == 1            
            if isTouchingDone == 1 && isCorrect(trial_count) == 0  
                submitInside = 0;
                isTouchingDone = 0;                
            end        
        end
        
%         if if_validHoldUntilCorrect == 0
%             %if isReleaseToHold == 1 && ifFreeTouch(trial_count) == 1
%             %if isHoldToRelease_notFirst == 1 && ifFreeTouch(trial_count) == 1
%             if submitInside == 1 && ifFreeTouch(trial_count) == 1
%                 isTouchingDone = 1;
%                 %isHoldTouchingDone = 1;
%                 %isReleaseToHold = 0;
%                 if isCorrect(trial_count) == 0
%                     isTouchingTimeOut(trial_count) = 1;
%                     isCorrect(trial_count) = -1;
%                     isTouchFrameCorrect(trial_count) = -1;
%                 end
%             end
%         end
%         
%         if if_validHoldUntilCorrect == 0
%             %             if isReleaseToHold == 1 && ifFreeTouch(trial_count) == 0 && isCorrect(trial_count) == 1% && isTouchingDone == 1
%             %if isReleaseToHold == 1 && ifFreeTouch(trial_count) == 0
%             %if isHoldToRelease_notFirst == 1 && ifFreeTouch(trial_count) == 0
%             if submitInside == 1 && ifFreeTouch(trial_count) == 0
%                 isTouchingDone = 1;
%                 if isCorrect(trial_count) == 0
%                     isTouchingTimeOut(trial_count) = 1;
%                     isCorrect(trial_count) = -1;
%                     isTouchFrameCorrect(trial_count) = -1;
%                 end
%             %elseif isReleaseToHold == 0 && ifFreeTouch(trial_count) == 0 && isCorrect(trial_count) == 1 && isTouchingDone == 1
%             %elseif isHoldToRelease_notFirst == 0 && ifFreeTouch(trial_count) == 0 && isCorrect(trial_count) == 1 && isTouchingDone == 1
%             elseif submitInside == 0 && ifFreeTouch(trial_count) == 0 && isCorrect(trial_count) == 1 && isTouchingDone == 1
%                 isTouchingDone = 0;
%             end
%             
%         elseif if_validHoldUntilCorrect == 1            
%             %if isReleaseToHold == 1 && ifFreeTouch(trial_count) == 0 && isCorrect(trial_count) == 1% && isTouchingDone == 1
%             %if isHoldToRelease_notFirst == 1 && ifFreeTouch(trial_count) == 0 && isCorrect(trial_count) == 1% && isTouchingDone == 1
%             %if submitInside == 1 && ifFreeTouch(trial_count) == 0 && isCorrect(trial_count) == 1% && isTouchingDone == 1
%             if submitInside == 1 && isCorrect(trial_count) == 1% && isTouchingDone == 1
%                 
%                 %         if isReleaseToHold == 1 && ifFreeTouch(trial_count) == 0
%                 isTouchingDone = 1;
%                 %if isCorrect(trial_count) == 0
%                 %   isCorrect(trial_count) = -1;
%                 %end
%             %elseif isReleaseToHold == 0 && ifFreeTouch(trial_count) == 0 && isCorrect(trial_count) == 1 && isTouchingDone == 1
%             %elseif isHoldToRelease_notFirst == 0 && ifFreeTouch(trial_count) == 0 && isCorrect(trial_count) == 1 && isTouchingDone == 1
%             elseif submitInside == 0 && ifFreeTouch(trial_count) == 0 && isCorrect(trial_count) == 1 && isTouchingDone == 1
%                 isTouchingDone = 0;
%             end
%         end

        
        if isTouchingDone == 1              
            if isCorrect(trial_count) == 0
                isTouchingTimeOut(trial_count) = 1;
                isCorrect(trial_count) = -1;
                isTouchFrameCorrect(trial_count) = -1;
            end
            
            %if isReleaseToHold == 1    
            %if isHoldToRelease_notFirst == 1 
            if submitInside == 1
                trial_para.endingHold_RT(trial_count) = toc(t0_isTouching)-last_freeTouch_toc;
                endingHold_RT = (toc(t0_isTouching)-last_freeTouch_toc) * 1000;
                fprintf("endingHold_RT = %.1f!\n", endingHold_RT);
                %a = toc-last_freeTouch_toc
            else
                trial_para.endingHold_RT(trial_count) = 0;
                endingHold_RT = 0;
            end
            
            % if error stop and this trial is error, then wait for some time so that submit marker can work
            if ifFreeTouch(trial_count) == 0 && isCorrect(trial_count) == -1
                WaitSecs(0.2);
                %fprintf('a\n');
            end
            
            m_marker.Data = uint8(2);% Inner trial markers
            trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms            
            %Eyelink('Message', 'Touching done');
            %message = sprintf('Touching done');
            message = sprintf('SUBMIT');
            Eyelink('Message', message);
            trial_para.marker.count = trial_para.marker.count + 1;
            m_markerCount.Data = int32(trial_para.marker.count);
            trial_para.marker.content = [trial_para.marker.content, {message}];            
            
            isTouching = 0;
            isItiInitial = 1;            
        end
    
    %---------------------------------------------------------------------------------------------    
    elseif isItiInitial == 1
        %% isItiInitial
        
        vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);
        vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);        
        
        tempWaterTime = initGiveWaterTime + floor(trial_count/100)*waterStepPer100;
        if ifSelection == 1
            if ifSelectOffloading == 0
                tempWaterTime = tempWaterTime * InternalMemoryWaterRatio;
            end
        end
        if isCorrect(trial_count) == 1                    
            %endBet_correctFast
            if endingHold_RT < endBet_fastSlowBoundary
                tempWaterTime = tempWaterTime * endBetRatio_correctFast;
                fprintf("endBet_correctFast!\n");
            %endBet_correctSlow
            else
                tempWaterTime = tempWaterTime * endBetRatio_correctSlow;
                fprintf("endBet_correctSlow!\n");
            end
        else
            %endBet_errorFast
            if endingHold_RT < endBet_fastSlowBoundary
                tempWaterTime = tempWaterTime * endBetRatio_errorFast;
                fprintf("endBet_errorFast!\n");
            %endBet_errorSlow
            else
                tempWaterTime = tempWaterTime * endBetRatio_errorSlow;
                fprintf("endBet_errorSlow!\n");
            end
        end
        if tempWaterTime > 20
            %home % To scrolling down in command window
            %WaitSecs(0.1);
            fprintf("==============================");
            %giveWater(tempWaterTime);
            
            
        old_data = m_marker.Data;
        m_marker.Data = uint8(123); 
        WaitSecs(0.1);
        %giveWater(tempWaterTime);
        giveWater(tempWaterTime,waterDropNum, waterDropInterval);
        m_marker.Data = uint8(old_data);            
            
            
            fprintf("\n");
        end
        
        if isCorrect(trial_count) == 1
%             tempWaterTime = initGiveWaterTime + floor(trial_count/100)*waterStepPer100;
%             if ifSelection == 0
%                 giveWater(tempWaterTime);
%             elseif ifSelection == 1
%                 if ifSelectOffloading == 1
%                     giveWater(tempWaterTime);
%                 elseif ifSelectOffloading == 0
%                     giveWater(tempWaterTime * InternalMemoryWaterRatio);
%                 end
%             end
            
            %home % To scrolling down in command window
            fprintf("A correct trial!\n");
            consecutiveCorrect = consecutiveCorrect + 1;
            cumulativeCorrect_greenRedSwitch = cumulativeCorrect_greenRedSwitch + 1;
            %fprintf("cumulativeCorrect_greenRedSwitch = %d, ", cumulativeCorrect_greenRedSwitch);
            fprintf("consecutiveCorrect = %d, ", consecutiveCorrect);
            fprintf("cumulativeError = %d\n", cumulativeError);    
            
            if cumulativeError == 0
                isIppatsuCorrect(trial_count) = 1;
                uniqueSeq_trial_count = uniqueSeq_trial_count + 1;
            else
                isIppatsuCorrect(trial_count) = -1;
            end
        elseif isCorrect(trial_count) == -1 
            if ifLinkTriangleSquare_inOneTrial == 1
                if isTriangletrialMiddle == 1
                    isTriangletrialMiddle = 0;
                    taskMode_noChoice_Green0_Red1 = 1;
                end
            end

            %home % To scrolling down in command window
            fprintf("An error trial!\n");
            
            if isTouchingTimeOut(trial_count) == 1
                fprintf("Touching is time out!\n");
            end
            cumulativeError = cumulativeError + 1;
            consecutiveCorrect = 0;
            %fprintf("cumulativeCorrect_greenRedSwitch = %d, ", cumulativeCorrect_greenRedSwitch);
            fprintf("consecutiveCorrect = %d, ", consecutiveCorrect);
            fprintf("cumulativeError = %d\n", cumulativeError);
            
            isIppatsuCorrect(trial_count) = -1;
            if cumulativeError == 1                
                uniqueSeq_trial_count = uniqueSeq_trial_count + 1;
            end
        end
        temp = currentSequence;
        %temp = relativeSequence{trial_count};
        if isCorrect(trial_count)==1
            sort_correctCount(temp) = sort_correctCount(temp) + 1;
        end
        sort_allCount(temp) = sort_allCount(temp) + 1;
        sort_isCorrect = sort_correctCount./sort_allCount;
        fprintf("sort_accuracy = ");
        for i=1:numFrames(trial_count)
            fprintf("%.2f ",sort_isCorrect(i));
        end
        fprintf("\n");
        
        % update sort_isCorrect_last40
        sort_correctCount_last40(2:40,:) = sort_correctCount_last40(1:39,:);
        sort_allCount_last40(2:40,:) = sort_allCount_last40(1:39,:);
        sort_correctCount_last40(1,:) = zeros(1,numFrames_rangeTail);
        sort_allCount_last40(1,:) = zeros(1,numFrames_rangeTail);        
        if isCorrect(trial_count)==1
            sort_correctCount_last40(1,currentSequence) = 1;
        end
        sort_allCount_last40(1,currentSequence) = 1;
        
        if size(sort_allCount_last40,2) ~= numFrames_rangeTail
            sort_correctCount_last40 = zeros(40, numFrames_rangeTail); %#ok<*PREALL>
            sort_allCount_last40 = zeros(40, numFrames_rangeTail);
            sort_correctCount_last40 = zeros(40, numFrames_rangeTail);
        end
        
        sort_accuracy_last40 = sum(sort_correctCount_last40,1)./sum(sort_allCount_last40,1);
        %fprintf("sort_accuracy_last40 = ");
        fprintf("sort_acc_40   = ");
        for i=1:numFrames(trial_count)
            fprintf("%.2f ",sort_accuracy_last40(i));
        end
        fprintf("\n");        
        
        
        if isManual_rFPW == 0
            isProbWeightChange = 1;
        end
        totalAccuracy = sum(isCorrect==1)/trial_count;
        ippatsuAccuracy = sum(isIppatsuCorrect==1)/uniqueSeq_trial_count;
        
        tempLen = seq_length_rangeTail-seq_length_rangeHead+1+1;
        if taskMode_noChoice_Green0_Red1 == 0
            internal_trial_count(tempLen) = internal_trial_count(tempLen) + 1;
            for tempi=1:tempLen-1
                if length(currentSequence) == tempi-1+seq_length_rangeHead
                    internal_trial_count(tempi) = internal_trial_count(tempi) + 1;                    
                end
            end                
            if isCorrect(trial_count)==1
                internal_correct_count(tempLen) = internal_correct_count(tempLen) + 1;
                for tempi=1:tempLen-1
                    if length(currentSequence) == tempi-1+seq_length_rangeHead
                        internal_correct_count(tempi) = internal_correct_count(tempi) + 1;
                    end
                end
            end
            
            if ifSeq_length_3to1 == 1
                internal_trial_count(1) = internal_trial_count(3) - internal_trial_count(2);
                internal_correct_count(1) = internal_correct_count(3) - internal_correct_count(2);
            end
            
            internalAccuracy = internal_correct_count./internal_trial_count;
            
            if selectFlag_newSeq(trial_count) == 1 && choiceCondition_flag(trial_count) == 2
                internal_trial_count_newSeq(tempLen) = internal_trial_count_newSeq(tempLen) + 1;
                for tempi=1:tempLen-1
                    if length(currentSequence) == tempi-1+seq_length_rangeHead
                        internal_trial_count_newSeq(tempi) = internal_trial_count_newSeq(tempi) + 1;
                    end
                end
                
                if ifSeq_length_3to1 == 1
                    internal_trial_count_newSeq(1) = internal_trial_count_newSeq(3) - internal_trial_count_newSeq(2);
                end
                                
                internal_trial_count_last20NewSeq(2:20) = internal_trial_count_last20NewSeq(1:19);
                offloading_trial_count_last20NewSeq(2:20) = offloading_trial_count_last20NewSeq(1:19);
                internal_trial_count_last20NewSeq(1) = 1;
                offloading_trial_count_last20NewSeq(1) = 0;
                
            end      
            
            % get internalNoChoiceAccuracy
            if choiceCondition_flag(trial_count) == 0
                internal_noChoice_trial_count(tempLen) = internal_noChoice_trial_count(tempLen) + 1;
                for tempi=1:tempLen-1
                    if length(currentSequence) == tempi-1+seq_length_rangeHead
                        internal_noChoice_trial_count(tempi) = internal_noChoice_trial_count(tempi) + 1;
                    end
                end
                if isCorrect(trial_count)==1
                    internal_noChoice_correct_count(tempLen) = internal_noChoice_correct_count(tempLen) + 1;
                    for tempi=1:tempLen-1
                        if length(currentSequence) == tempi-1+seq_length_rangeHead
                            internal_noChoice_correct_count(tempi) = internal_noChoice_correct_count(tempi) + 1;
                        end
                    end
                end
                
                if ifSeq_length_3to1 == 1
                    internal_noChoice_trial_count(1) = internal_noChoice_trial_count(3) - internal_noChoice_trial_count(2);
                    internal_noChoice_correct_count(1) = internal_noChoice_correct_count(3) - internal_noChoice_correct_count(2);
                end
                
                internalNoChoiceAccuracy = internal_noChoice_correct_count./internal_noChoice_trial_count;
            end            
            
        elseif taskMode_noChoice_Green0_Red1 == 1
            offloading_trial_count(tempLen) = offloading_trial_count(tempLen) + 1;
            for tempi=1:tempLen-1
                if length(currentSequence) == tempi-1+seq_length_rangeHead
                    offloading_trial_count(tempi) = offloading_trial_count(tempi) + 1;                    
                end
            end              
            if isCorrect(trial_count)==1
                offloading_correct_count(tempLen) = offloading_correct_count(tempLen) + 1;
                for tempi=1:tempLen-1
                    if length(currentSequence) == tempi-1+seq_length_rangeHead
                        offloading_correct_count(tempi) = offloading_correct_count(tempi) + 1;
                    end
                end                
            end
            
            if ifSeq_length_3to1 == 1
                offloading_trial_count(1) = offloading_trial_count(3) - offloading_trial_count(2);
                offloading_correct_count(1) = offloading_correct_count(3) - offloading_correct_count(2);
            end
            
            offloadingAccuracy = offloading_correct_count./offloading_trial_count;  
            
            if selectFlag_newSeq(trial_count) == 1 && choiceCondition_flag(trial_count) == 2
                offloading_trial_count_newSeq(tempLen) = offloading_trial_count_newSeq(tempLen) + 1;
                for tempi=1:tempLen-1
                    if length(currentSequence) == tempi-1+seq_length_rangeHead
                        offloading_trial_count_newSeq(tempi) = offloading_trial_count_newSeq(tempi) + 1;
                    end
                end
                
                if ifSeq_length_3to1 == 1
                    offloading_trial_count_newSeq(1) = offloading_trial_count_newSeq(3) - offloading_trial_count_newSeq(2);
                end                
                
                internal_trial_count_last20NewSeq(2:20) = internal_trial_count_last20NewSeq(1:19);
                offloading_trial_count_last20NewSeq(2:20) = offloading_trial_count_last20NewSeq(1:19);
                internal_trial_count_last20NewSeq(1) = 0;   
                offloading_trial_count_last20NewSeq(1) = 1;                
                
            end            
        end
        
        offloadingProb = offloading_trial_count_newSeq./(offloading_trial_count_newSeq+internal_trial_count_newSeq);
       
        internal_trial_count_last20NewSeq(21) = sum(internal_trial_count_last20NewSeq(1:20));
        offloading_trial_count_last20NewSeq(21) = sum(offloading_trial_count_last20NewSeq(1:20));
        last20NewSeq_rProb = offloading_trial_count_last20NewSeq(21)./(offloading_trial_count_last20NewSeq(21)+internal_trial_count_last20NewSeq(21));
        
        
        if trial_count >= 100
            last100Accuracy = sum(isCorrect(trial_count-99:trial_count)==1)/100;
        else
            last100Accuracy = -1;
        end    
        
        if trial_count >= 20
            last20Accuracy = sum(isCorrect(trial_count-19:trial_count)==1)/20;
        else
            last20Accuracy = -1;
        end
        
%         if trial_count >= 10
%             last10Accuracy = sum(isCorrect(trial_count-9:trial_count)==1)/10;
%         else
%             last10Accuracy = -1;
%         end
        
        fprintf("totalAccuracy = %.3f, ippatsuAccuracy = %.3f\n", totalAccuracy, ippatsuAccuracy);
        trial_para.totalAccuracy = totalAccuracy;
        trial_para.ippatsuAccuracy = ippatsuAccuracy;
        %fprintf("last100Accuracy = %.3f, last20Accuracy = %.3f, last10Accuracy = %.3f\n", last100Accuracy, last20Accuracy, last10Accuracy);
        
        tempLen = seq_length_rangeTail-seq_length_rangeHead+1+1;
        fprintf("gAccuracy = ");
        if tempLen-1 ~= 1
            fprintf("%.2f(total length), ", internalAccuracy(tempLen));
        end
        for tempi=1:tempLen-1         
            fprintf("%.2f, ", internalAccuracy(tempi));
        end
        if tempLen-1 ~= 1
        fprintf("\n");
        end
        trial_para.internalAccuracy = internalAccuracy;
        
        %noChoice gAccuracy
        fprintf("noChoice  = ");
        if tempLen-1 ~= 1
            fprintf("%.2f(total length), ", internalNoChoiceAccuracy(tempLen));
        end
        for tempi=1:tempLen-1         
            fprintf("%.2f, ", internalNoChoiceAccuracy(tempi));
        end
        if tempLen-1 ~= 1
        fprintf("\n");
        end
        trial_para.internalNoChoiceAccuracy = internalNoChoiceAccuracy;        
        
        
        fprintf("rAccuracy = ");
        if tempLen-1 ~= 1
            fprintf("%.2f(total length), ", offloadingAccuracy(tempLen));
        end
        for tempi=1:tempLen-1
            fprintf("%.2f, ", offloadingAccuracy(tempi));
        end
        fprintf("\n");
        trial_para.offloadingAccuracy = offloadingAccuracy;
        
        fprintf("rProb = ");
        if tempLen-1 ~= 1
            fprintf("%.2f(total length), ", offloadingProb(tempLen));
        end
        for tempi=1:tempLen-1
            fprintf("%.2f, ", offloadingProb(tempi));
        end
        fprintf("\n");
        trial_para.offloadingProb = offloadingProb;
         
        fprintf("last20NewSeq_rProb = %.2f\n", last20NewSeq_rProb);
        trial_para.last20NewSeq_rProb = last20NewSeq_rProb;
        
        %fprintf("gAccuracy = %.3f, rAccuracy = %.3f\n", internalAccuracy, offloadingAccuracy);

        fprintf("last100Accuracy = %.3f, last20Accuracy = %.3f\n", last100Accuracy, last20Accuracy);
        trial_para.last100Accuracy = last100Accuracy;
        trial_para.last20Accuracy = last20Accuracy;
        
        tempLength = find(seq_length_rangeHead:seq_length_rangeTail == length(currentSequence));
        seqLength_trial_count(tempLength) = seqLength_trial_count(tempLength) + 1;
        seqLength_trial_count(end) = seqLength_trial_count(end) + 1;
        breakRate = break_sumCount./(break_sumCount+seqLength_trial_count);
        
        allbreak = (break_count(trial_count))...
                ./(break_count(trial_count) + trial_count);
            
        if trial_count >= 100
            last100break = (break_count(trial_count) - break_count(trial_count-100+1))...
                ./( (break_count(trial_count) - break_count(trial_count-100+1)) + 100);
        else
            last100break = -1;
        end    
        
        if trial_count >= 40
            last40break = (break_count(trial_count) - break_count(trial_count-40+1))...
                ./( (break_count(trial_count) - break_count(trial_count-40+1)) + 40);
        else
            last40break = -1;
        end    
        
        if floor(trial_count/dynamicFix)*dynamicFix == trial_count
            lastXbreak = (break_count(trial_count) - break_count(trial_count-dynamicFix+1))...
                ./( (break_count(trial_count) - break_count(trial_count-dynamicFix+1)) + dynamicFix);
        else
            lastXbreak = -1;            
        end           
        if trial_count < dynamicFix
            tempLastXbreak = -1;
        end
        fprintf("allbreak = %.3f, last100break = %.3f, last40break = %.3f, lastXbreak = %.3f\n",...
            allbreak, last100break, last40break, tempLastXbreak);
        
        fprintf("breakRate = ");
        if tempLen-1 ~= 1
            fprintf("%.2f(total length), ", breakRate(tempLen));
        end
        for tempi=1:tempLen-1
            fprintf("%.2f, ", breakRate(tempi));
        end
        fprintf("\n");

        
        if lastXbreak > dynamicFixThreshold
            tempLastXbreak = lastXbreak;
            lastXbreak = -1;
            
            fixationTime = fixationTime - 24;
            if fixationTime(trial_count) < 0
                fixationTime = zeros(1, trial_num);
            end
            fprintf("fixationTime = %d\n",mean(fixationTime));
        elseif lastXbreak > 0 && lastXbreak <= dynamicFixThreshold
            tempLastXbreak = lastXbreak;
            lastXbreak = -1;
            
            fixationTime = fixationTime + 60;
            fprintf("fixationTime = %d\n",mean(fixationTime));
        end
        
        fprintf("trial_count = %d\n",trial_count);
        if if_seqList == 1    
            if if_length1_alone == 0
                temp_trialCount_oneRun = mod(trial_count,temp_trialNum_oneRun);
                if temp_trialCount_oneRun == 0
                    temp_trialCount_oneRun = temp_trialNum_oneRun;
                end
                temp_runCount = ceil(trial_count/temp_trialNum_oneRun);
                fprintf("%d/%d in run %d.\n",temp_trialCount_oneRun, temp_trialNum_oneRun, temp_runCount);
            elseif if_length1_alone == 1
                trialNum_length1Run;
                if trial_count <= trialNum_length1Run
                    temp_trialCount_oneRun = trial_count;
                    temp_runCount = 1;
                    fprintf("%d/%d in run %d.\n",temp_trialCount_oneRun, trialNum_length1Run, temp_runCount);
                    
                elseif trial_count > trialNum_length1Run
                    temp_trial_count = trial_count - trialNum_length1Run;
                    
                    temp_trialCount_oneRun = mod(temp_trial_count,temp_trialNum_oneRun);
                    if temp_trialCount_oneRun == 0
                        temp_trialCount_oneRun = temp_trialNum_oneRun;
                    end
                    temp_runCount = ceil(temp_trial_count/temp_trialNum_oneRun)+1;
                    fprintf("%d/%d in run %d.\n",temp_trialCount_oneRun, temp_trialNum_oneRun, temp_runCount);                    
                end
                
            end
        end
        
        isFilled;
        %currentResponse = find(isFilled==0);
        currentResponse = find(isFilled(1:numFrames_rangeTail)==0);
        fprintf("=================currentResponse = ");
        for tempi=1:length(currentResponse)
            fprintf("%d ", currentResponse(tempi));
        end
        fprintf("\n");
        
%         if isWaitGocueBreak(trial_count) == 1
%             fprintf("It is a break trial!\n");
%         end
        fprintf("----------------------------A trial end---------------------------\n\n");
        
        %vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);
        %vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);

        trial_para.trial_count = trial_count;
        trial_para.uniqueSeq_trial_count = uniqueSeq_trial_count;
        trial_para.isCorrect(trial_count) = isCorrect(trial_count);        
        trial_para.isIppatsuCorrect(trial_count) = isIppatsuCorrect(trial_count);
        %trial_para.isWaitGocueBreak(trial_count) = isWaitGocueBreak(trial_count);        
        trial_para.isTouchFrameCorrect(trial_count) = isTouchFrameCorrect(trial_count);
        trial_para.isTouchingTimeOut(trial_count) = isTouchingTimeOut(trial_count);
        trial_para.isBackTouch(trial_count) = isBackTouch(trial_count);
        
        trial_para.cumulativeCorrect_greenRedSwitch(trial_count) = cumulativeCorrect_greenRedSwitch;
        trial_para.consecutiveCorrect(trial_count) = consecutiveCorrect;
        trial_para.cumulativeError(trial_count) = cumulativeError;
        trial_para.currentFrames(trial_count) = {currentFrames};
        trial_para.currentSequence(trial_count) = {currentSequence};
        trial_para.seq_length(trial_count) = seq_length(trial_count);
        trial_para.isFilled(trial_count) = {isFilled(1:numFrames_rangeTail)};
        
        trial_para.taskMode_noChoice_Green0_Red1(trial_count) = taskMode_noChoice_Green0_Red1;
        trial_para.selectFlag_newSeq(trial_count) = selectFlag_newSeq(trial_count);
        trial_para.choiceMode(trial_count) = choiceMode;
        trial_para.choiceCondition_flag(trial_count) = choiceCondition_flag(trial_count);
        
        trial_para.reCenterGazeBias_x(trial_count) = reCenterGazeBias_x;
        trial_para.reCenterGazeBias_y(trial_count) = reCenterGazeBias_y;


        if trial_count >= 1
            save(filename,'basic_para');
            save(filename,'trial_para','-append');
        end
        
        t0_isIti = tic;
        isItiInitial = 0;
        isIti = 1;
          
        isPushingScreen = 0;
        
%         if isCorrect(trial_count) == 1 %#ok<*IFBDUP>
%             tempIti = ItiTime(trial_count);
%         else
%             %tempIti = ItiTime(trial_count)*2;
%             %tempIti = ItiTime(trial_count)*1.5;
%             tempIti = ItiTime(trial_count)/3;
%         end        
        %fprintf("ItiTime = %d ms\n",round(tempIti));
        if if_eyeResponse == 0
            while TouchEventAvail(dev) > 0
                TouchEventGet(dev, window);
            end
        end
    elseif isItiInitial_break == 1
        %% isItiInitial_break
        break_count(trial_count) = break_count(trial_count) + 1;
        tempLength = find(seq_length_rangeHead:seq_length_rangeTail == length(currentSequence));
        break_sumCount(tempLength) = break_sumCount(tempLength) + 1;
        break_sumCount(end) = break_sumCount(end) + 1;
        
        trial_para.isWaitGocueBreak(trial_count) = isWaitGocueBreak(trial_count); 
        
        fprintf("\n================================\n");
        fprintf("There is a break !\n");
        fprintf("break_count = %d !\n", break_count(trial_count));
        
        vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);
        vbl  = Screen('Flip', window, vbl + (shortWaitFrames - 0.5) * ifi);
        
        %breakTime(break_count(trial_count)) = toc(t0_showReady) * 1000;%#ok<SAGROW> %ms
        breakTime(break_count(trial_count)) = toc(t0_isShowing) * 1000;%#ok<SAGROW> %ms        
        breakSeqLength(break_count(trial_count)) = seq_length(trial_count); %#ok<SAGROW>
        trial_para.breakTime(break_count(trial_count)) = breakTime(break_count(trial_count));
        trial_para.breakSeqLength(break_count(trial_count)) = breakSeqLength(break_count(trial_count));
        fprintf("breakTime = %.1fms, mean(breakTime) = %.1fms!\n",...
            breakTime(break_count(trial_count)), mean(breakTime));
       
        
        if trial_count >= 1
            save(filename,'basic_para');
            save(filename,'trial_para','-append');
        end        
        
        t0_isIti = tic;
        isItiInitial_break = 0;
        isIti = 1;
                
        isPushingScreen = 0;
        
        %tempIti = ItiTime(trial_count)*4;
        %tempIti = ItiTime(trial_count);
%         tempIti = ItiTime(trial_count)/6;
        %tempIti = ItiTime(trial_count)*0.5;
        
        %fprintf("ItiTime = %d ms\n",round(tempIti));
        if if_eyeResponse == 0
            while TouchEventAvail(dev) > 0
                TouchEventGet(dev, window);
            end
        end
    %---------------------------------------------------------------------------------------------    
    elseif isIti == 1
        %% isIti                
%         if isCorrect(trial_count) == 1
%             %tempIti = ItiTime(trial_count)/2.4;
%             %tempIti = ItiTime(trial_count)/4;
%             tempIti = ItiTime(trial_count)/2;
%             %tempIti = ItiTime(trial_count);
%         elseif isCorrect(trial_count) == -1
%             %tempIti = ItiTime(trial_count)/6;
%             tempIti = ItiTime(trial_count)/3;
%         elseif isCorrect(trial_count) == 0
%             %tempIti = ItiTime(trial_count)/6;
%             tempIti = ItiTime(trial_count)/3;
%         end
        tempIti = ItiTime(trial_count);
        
        if isReleaseToHold == 1
            %fprintf('Iti pushiment begins! Iti = %d\n', tempIti);   
            fprintf('Iti begins! Iti = %.1f\n', tempIti);
        end
        if isHoldToRelease == 1
            fprintf('Iti pushiment ends!\n');
        end            
                
        isPushingScreen = 0;
        if if_eyeResponse == 0
            if TouchEventAvail(dev) > 0
                isPushingScreen = 1;
                %beep and sound could cause matlab crash
                %in Ubuntu 20.04, so hide these code.
                %beep
                %sound(sin(2*pi*25*(1:4000)/100));
                %sound(sin(2*pi*25*(100:700)/1000));
                pause(0.05)
                while TouchEventAvail(dev) > 0
                    TouchEventGet(dev, window);
                end
            end
        end
                
        if isHold == 1 || isPushingScreen == 1
            t0_isIti = tic;            
            isPushingScreen = 0;
        end        
        %tempIti = ItiTime(trial_count) + ItiTime(trial_count)*(isCorrect(trial_count)*-1+1)/2;
        
%         if isCorrect(trial_count) == 1 || isWaitGocueBreak(trial_count) == 1
%             tempIti = ItiTime(trial_count);
%         else
%             tempIti = ItiTime(trial_count)*2;
%         end

        if toc(t0_isIti) > tempIti/1000
            fprintf('Iti ends! Iti = %.1f\n', toc(t0_isIti)*1000);
            isIti = 0;
            isShowingInitial = 1;
            initialFlag = 0;            
            if if_freeview0_fixation1 == 1
                t0_fixOn = tic;                
            end                 
            
            if isWaitGocueBreak(trial_count) == -1
                trial_count = trial_count + 1;
                break_count(trial_count) = break_count(trial_count-1);                
            end
            
            temp_validSample = ceil(length(currentTrial_fixPos.x)/3):ceil(length(currentTrial_fixPos.x)*2/3);
            if sum(temp_validSample) == 0
                reCenterGazeBias_x_offset = 0;
                reCenterGazeBias_y_offset = 0;
            else
                reCenterGazeBias_x_offset = 0.2*(crossCenter(1) - mean(currentTrial_fixPos.x(temp_validSample)));
                reCenterGazeBias_y_offset = 0.2*(crossCenter(2) - mean(currentTrial_fixPos.y(temp_validSample)));
            end
            if isnan(reCenterGazeBias_x_offset) || isnan(reCenterGazeBias_y_offset)
                reCenterGazeBias_x_offset = 0;
                reCenterGazeBias_y_offset = 0;
            end            
            reCenterGazeBias_x = reCenterGazeBias_x + reCenterGazeBias_x_offset;
            reCenterGazeBias_y = reCenterGazeBias_y + reCenterGazeBias_y_offset;
            currentTrial_fixPos.x = [];
            currentTrial_fixPos.y = [];
            gazeBias_offset_count = gazeBias_offset_count + 1;
            if gazeBias_offset_count >= 10
                gazeBias_offset_count = 0;
                ELHost_changeFLag = 1;
            end            
            fprintf('reCenterGazeBias_x_offset=%.2f,reCenterGazeBias_y_offset=%.2f\n',...
                reCenterGazeBias_x_offset,reCenterGazeBias_y_offset);            
            
            if if_eyeResponse == 0
                TouchQueueStop(dev);
                TouchQueueRelease(dev);
                
                TouchQueueCreate(window, dev);
                TouchQueueStart(dev);
            end
        end              
    end
    %---------------------------------------------------------------------------------------------
    % dequeues all samples and events from the link. AKA. Clear Buffer!    
    [samplesIn, eventsIn, drained] = Eyelink('GetQueuedData');
end % End of while loop
%% End
Priority(0);

vbl  = Screen('Flip', window, vbl + 0.5 * ifi);

WaitSecs(0.5);

% Stop marker signal
m_marker.Data = uint8(255);% End signal
trial_para.marker.time = [trial_para.marker.time, GetSecs*1000];% save time in ms
message = sprintf('END');
Eyelink('Message', message);
trial_para.marker.count = trial_para.marker.count + 1;
m_markerCount.Data = int32(trial_para.marker.count);
trial_para.marker.content = [trial_para.marker.content, {message}];

m_markerCount.Data = int32(0);
WaitSecs(3);
m_marker.Data = uint8(128);% End signal 2



if trial_count >= 1
    save(filename,'basic_para');
    save(filename,'trial_para','-append');
end

% Save reCenterGazeBias
temp_fileName = [path,'/',searchName_reCenterGazeBias];
save(temp_fileName, 'reCenterGazeBias_x');
save(temp_fileName, 'reCenterGazeBias_y', '-append');


% Stop recording eye movements at the end of each trial
WaitSecs(0.1); % Add 100 msec of data to catch final events before stopping
Eyelink('StopRecording'); % Stop tracker recording

%% Eyelink STEP 6: CLOSE EDF FILE. TRANSFER EDF COPY TO DISPLAY PC. CLOSE EYELINK CONNECTION. FINISH UP
Eyelink('SetOfflineMode'); % Put tracker in idle/offline mode
Eyelink('Command', 'clear_screen 0'); % Clear Host PC backdrop graphics at the end of the experiment
WaitSecs(0.5); % Allow some time before closing and transferring file
Eyelink('CloseFile'); % Close EDF file on Host PC
% Transfer a copy of the EDF file to Display PC
% transferFile; % See transferFile function below
%%%%%%%%%%%%%%%%%%%%%%%% transferFile %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dummymode ==0 % If connected to EyeLink
    if trial_count > 1
        % Show 'Receiving data file...' text until file transfer is complete
        Screen('FillRect', window, el.backgroundcolour); % Prepare background on backbuffer
        Screen('DrawText', window, 'Receiving data file...', 5, screenYpixels-35, 0); % Prepare text
        Screen('Flip', window); % Present text
        fprintf('Receiving data file ''%s.edf''\n', edfFile); % Print some text in Matlab's Command Window
        
        % Transfer EDF file to Host PC
        % [status =] Eyelink('ReceiveFile',['src'], ['dest'], ['dest_is_path'])
        %path = [pwd, '/MATDATA/'];
        %status = Eyelink('ReceiveFile', [], path, 1);
        dest = ['edf',datestr(now, 'yyyy-mm-dd'), mfilename, '-', num2str(filenum), '.edf'];
        cd './EDF'
        status = Eyelink('ReceiveFile', [], dest);
        cd ..
        
        % Check if EDF file has been transferred successfully and print file size in Matlab's Command Window
        if status > 0
            fprintf('EDF file size: %.1f KB\n', status/1024); % Divide file size by 1024 to convert bytes to KB
        end
        % Print transferred EDF file path in Matlab's Command Window
        fprintf('Data file ''%s.edf'' can be found in ''%s''\n', edfFile, pwd);
    else
        fprintf('No EDF file saved because trial_count <= 1!\n');
    end
else
    fprintf('No EDF file saved in Dummy mode\n');
end

% Cleanup function used throughout the script above
try %#ok<TRYNC>
    Screen('CloseAll'); % Close window if it is open
end
Eyelink('Shutdown'); % Close EyeLink connection
% ListenChar(0); % Restore keyboard output to Matlab



% Clear the screen
sca;
if trial_count > 1
    fprintf('Current filename is %s\n', shortfilename);
else
    fprintf('No file saved because trial_count <= 1!\n');
end
%update log: edit in 2020-04-09, the first day back to lab after COVID-19
%01-07 update: different ITITime based on whether the trial is correct or
%not
%01-08 update: Auto rFPW 
%01-10 update: (1)frameMode;(2)bias_x = 60+150 --> 60+150-75
%01-12 update: change StimulusExtendTime 
%01-13 update: moving touch now is available
%01-14 update: detectBaseRect = 145 --> 120
%04-13 update: waterStepPer100 = 10 --> 8
%04-13 update: make change of numFrames more compatible with the code
%04-19 update: screenYpixels_bias = 50 --> 0
%04-20 update: detection method: rectangle --> rectangle & radiant
%04-20 update: multi-point touch input is available while the lever is
%under hold state
%04-21 update: radius = 486 --> 550, hope to reduce mis-touch of monkey
%04-21 update: forbid move touch, to modify touch gesture
%04-22 update: allow move touch again
%04-23 update: screenXpixels_bias = 60+150-75 --> 0
%04-24 update: waterStepPer100 = 8 --> 10 --> 12
%04-24 update: initGiveWaterTime = 160 --> 180 --> 190 --> 180
%04-24 update: add encourageStimET = 10000 and encourageCumuErrorLimit = 5,
%to encourage monkey with simplify the current task while cumulative errors
%04-28 update: add ippatsuAccuracy to calculate one shot accuracy
%05-01 update: InitAngle_arc = 15 --> -70 --> -45
%05-01 update: radius = 486 -->465
%05-01 update: screenYpixels_bias = -140 --> -125
%05-05 update: make seqs clockwise(red fill)
%05-10 update: showBaseRect = 100 --> 135
%05-10 update: detectBaseRect = 175 --> 135
%05-13 update: add white dot to show where touch point is, and fix
%a potential detection bug
%05-19 update: add firstDetectRectModifed to make the first detect rect
%larger than others so that monkey could touch more correctly.(hope so)
%06-12 update: add ifLinkTriangleSquare mechnism corresponding with
%step7)linkTriangleSquare
%06-14 update: add ifLinkTriangleSquare_inOneTrial mechnism corresponding with
%step8)linkTriangleSquare_inOneTrial
%06-15 update: let triangle and square frame&fill transfer to circle, so
%that unify the shape(step 9)
%06-21 update: modify the calculation of ippatsuAccuracy by uniqueSeq_trial_count
%06-21 update: add break_count so that break trial will not be calculated
%for kinds of accuracy
%06-29 update: ignore touching currentMiddleSeq when offloading, hope to
%improve accuracy of monkey, ifAllowRedMiddle
%07-02 update: add pointShowPWM so that there is a short break between to
%red point during stimulus period
%07-12 update: use memmapfile to quickly get joystick state
%07-22 update: add two if condition to fix up a break bug
%07-23 update: make detection field cirle
%07-23 update: set mytouchQueueBuffer to calculate mean position of touch
%points, so that monkey could hit the target more easily (hope to solve motor problem)
%08-02 update: add choiceMode, to represent force-choice-g, force-choice-r
%and free-choice-g/r with 0,1,2
%08-02 update: add choiceMode_backup so that when current seq don't update
%yet, choiceMode would temporarily become 0 or 1 from 2
%08-02 update: set mytouchQueueBuffer_selection to calculate mean position of touch
%points, just same as update in 07-23
%08-05 update: change the calculation of rProb, only involve trial under
%choiceMode == 2 condition
%08-07 update: add another TouchQueueCreate-TouchQueueRelease-pair so that
%there are two T-T-pair, one for isSelecting, another for isTouching. Hope
%to solve recent touch screen bug
%08-14 update: add last20NewSeq_rProb
%08-21 update: only permit the first touch point to touch, as well as preserve 
%mytouchQueueBuffer to calculate mean position of touch. hope to improve
%motor performance
%08-23 update: add choiceCondition_boardLine1 choiceCondition_boardLine2 to
%add some random cases into g/r condition, so that there are some
%proportion of g and r. In a word, g:r:g/r == boardLine1:boardLine2-boardLine1:1-boardLine2
%08-31 update: add cumulativeCorrectLimit_seqLengthSwitch to flexibly control
%sequence switch
%10-28 update: add ifSeqLength14_sparse to reduce frequency of seq_length ==
%1 trial
%11-09 update: add ifSpecializedTraining_R to conduct specialized training
%of R (offloading) so that relative parameters will change
%Skip 2021
%2022-06-13 update: add more eyelink marker, and use the same time params
%as previous behavioral recording experiment of Zelku.
%2022-08-11 update: add 'continue' command when detect break. This bug could occur when break happens in the
%last frame of current stage, which leads to the next stage and even the end of the tiral with a break flag!
