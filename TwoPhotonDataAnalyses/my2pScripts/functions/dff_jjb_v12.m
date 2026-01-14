function [F_dff,F0] = dff_jjb_v12(F,df_window,path)
% Modified from detrend_df_f_auto function of CaImAn-MATLAB

%% Initialization
estimate_percentile_levelName = 'estimate_percentile_level_jjb';
estimate_percentile_level_v = autoGetFunName_myScripts(estimate_percentile_levelName, [path,'\','functions']);
fun_estimate_percentile_level = str2func(estimate_percentile_level_v);
% fprintf('Choose to run %s.\n',estimate_percentile_level_v);

prctfiltName = 'prctfilt_jjb';
prctfilt_v = autoGetFunName_myScripts(prctfiltName, [path,'\','functions']);
fun_prctfilt = str2func(prctfilt_v);
% fprintf('Choose to run %s.\n',prctfilt_v);


[roiNum,T] = size(F);
F0 = zeros(roiNum,T);
% roiNum = 300;

pause(3);% wait for 2 seconds, to avoid cpu overheat

%% Main loop
parfor tempi = 1:roiNum
% parfor tempi = 251:276
% for tempi = 1:roiNum
    if isempty(df_window) || (df_window > T)
        [~,cdf] = fun_estimate_percentile_level(F(tempi,:),T,T); %#ok<*PFBNS>
        cdf_level = median(cdf);
        temp_F0 = prctile(F(tempi,:),cdf_level,2);
        F0(tempi,:) = repmat(temp_F0,1,T);
    else
        %if sum(F(tempi,:)) == 0
        %    a = 1;
        %end
        a = 1; %#ok<*NASGU>
        
        temp_F_raw = F(tempi,:);
        windowSize = 10;
        temp_F = smoothdata(temp_F_raw,2,'rloess',windowSize);
        %temp_F = smoothdata(temp_F_raw,2,'gaussian',windowSize);
        %temp_F = smoothdata(temp_F_raw,2,'sgolay' ,windowSize);
        
        
        
        [~,cdf] = fun_estimate_percentile_level(temp_F,df_window,round(df_window/2));
        cdf_level = median(cdf);
        %F0(tempi,:) = fun_prctfilt(F(tempi,:),cdf_level,df_window);
        temp_F0 = fun_prctfilt(temp_F,cdf_level,df_window);        
        
        if isnan(cdf_level) == 1
            cdf_level = 50;
        end
        
        temp_prctile = cdf_level;%50        
        while true
            if sum(temp_F0==1) == 0
                break
            end
            
            if temp_prctile > 100
                break
            end   
            
            cdf_level = temp_prctile;
            %temp_F0 = fun_prctfilt(F(tempi,:),cdf_level,df_window);
            temp_F0 = fun_prctfilt(temp_F,cdf_level,df_window);

            temp_prctile = temp_prctile + 25;         
        end
        
%         temp_prctile = cdf_level;%50
%         warningFlag = 0;
%         breakFlag = 0;
%         %warningCount = 0;
%         while breakFlag == 0
%             a = 2;
%             if sum(temp_F0==1) == 0
%                 breakFlag = 1;
%             end
%             if warningFlag == 1
%                 %fprintf('Dff warning!\n');
%                 %warningCount = warningCount + 1;
%                 %break
%                 breakFlag = 1;
%             end
%             %temp_prctile = temp_prctile + 1;
%             temp_prctile = temp_prctile + 10;
%             if temp_prctile > 100
%                 %temp_prctile = 100;
%                 warningFlag = 1; 
%             else
%                 %cdf_level = prctile(cdf,temp_prctile);
%                 cdf_level = temp_prctile;
%                 temp_F0 = fun_prctfilt(F(tempi,:),cdf_level,df_window);
%             end
%         end
%         
%         %if warningFlag == 1
%         %    fprintf('Dff warning number = %d!\n',warningCount);
%         %end
        
        
        
        F0(tempi,:) = temp_F0;
        
    end
end
F_dff = (F-F0)./F0;
F0 = single(F0);



%% END
