function categorySomaDendrite_v5(path)

if(~exist('path','var'))
    path = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230123A_Ding_Site09B_sameFOV0122\Result20230205T170425\plane0';    
end
fprintf('Now is runing categorySomaDendrite.\n');

fileName_Fall = 'Fall.mat';
fileName_iscell = 'iscell.npy';

fullFileName_Fall = [path,'\',fileName_Fall];
fullFileName_iscell = [path,'\',fileName_iscell];

load(fullFileName_Fall,'stat');
iscell = readNPY(fullFileName_iscell);

a = 1;
isdendrite = (iscell(:,2) == 1.5);
issoma = (iscell(:,2) == 2);

padSize = 15;

for tempi=1:length(stat)
    if isdendrite(tempi) == true || issoma(tempi) == true
        temp_stat = stat{tempi};
        temp_index = double([temp_stat.xpix;temp_stat.ypix]);
        I = false(512,512);
        for tempj=1:length(temp_index)
            I(temp_index(2,tempj),temp_index(1,tempj)) = true; % suite2p order is [y,x]
        end
        
        [x_min,x_max] = bounds(double(temp_stat.xpix));
        [y_min,y_max] = bounds(double(temp_stat.ypix));
        
        x_min = max(1,x_min-padSize);
        x_max = min(512,x_max+padSize);
        y_min = max(1,y_min-padSize);
        y_max = min(512,y_max+padSize);
        
        I0 = I(y_min:y_max,x_min:x_max); % suite2p order is [y,x]

        %BW0 = bwmorph(I0,'close');
        BW0 = I0;
        BW1 = imfill(BW0,'holes');

        target_npix = 50;%50-->40
        scale = max(1,target_npix/double(temp_stat.npix));% only for zoom in
        BW2 = imresize(BW1,scale,'Method','nearest');

        BW3 = bwmorph(BW2,'remove');
        BW4 = bwmorph(BW2,'skel',Inf);
        BW5 = (imdilate(BW4,ones(3,3))) & BW2;

        npix_full = sum(BW2,'all');
        npix_contour = sum(BW3,'all');
        npix_skel = sum(BW4,'all');        
        npix_skel_dilate = sum(BW5,'all');
        
        
        npix_PWM = (npix_skel+npix_contour)/npix_full;
        npix_PWM2 = (npix_skel_dilate)/npix_full;
        
        props = regionprops('table',BW2,'MajorAxisLength','MinorAxisLength');
        axisRatio = props.MajorAxisLength./props.MinorAxisLength;
        [~,temp_axisRatioIndex] = max(props.MajorAxisLength);
        axisRatio = axisRatio(temp_axisRatioIndex);

        
        if tempi == 294+1
           a = 1; 
        end
        
        % move from dendrite to soma: npix_PWM<0.8
        % move from soma to dendrite: npix_PWM>1
        if isdendrite(tempi) == true
            if double(temp_stat.npix) > target_npix
            %if double(temp_stat.npix) > 100
                %if npix_PWM < 0.85 %0.8
                if npix_PWM2 < 0.85 %0.8
                    iscell(tempi,2) = 2;
                end
            else
                if axisRatio < 2.5
                    iscell(tempi,2) = 2;
                end
            end
            
        elseif issoma(tempi) == true
            if double(temp_stat.npix) > target_npix
            %if double(temp_stat.npix) > 100
                %if npix_PWM > 0.85 %0.8
                if npix_PWM2 > 0.85 %0.8
                    iscell(tempi,2) = 1.5;
                end
            else
                if axisRatio > 2.5
                    iscell(tempi,2) = 1.5;
                end
            end
            
        end
        
    else
        continue
    end
end
writeNPY(iscell, fullFileName_iscell);

%% END
