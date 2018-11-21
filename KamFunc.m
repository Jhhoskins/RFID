function [OutDist_Phase_mat, OutDist_Dop_mat, MinPhaseOutLength, MinDopOutLength] = KamFunc(filename, ant_num)
% Clearing
% clc
% clear all
% close all

% Hard code loading the file and saving experiment number
tag = load(filename);
save(['Exp' num2str(1) '.mat'],'tag')

%% Opening and Saving multiple consecutive experiments
% ExpNo = 1100;
% UserNo = '02'; % or the scenario number
% for i = ExpNo+1:2:ExpNo+11
%     tag = load(['C:\Users\alikamr3\Desktop\RFIDs_Research\FinalSinglePersonImaging\FinalExpSinglePersonTrainProd_' num2str(i-ExpNo) '_' num2str(i+1-ExpNo) '_' UserNo '.csv']);
%     save(['Exp' num2str(i) '.mat'],'tag')
% end

%% Reading tag values
MovingAvgWindow = [10/5 15/5 15/5]*5;           %just an array of [10 15 15]. Don't know why there's scalar mult. Doesn't seem like types are diff
%DopplerOn = 0;                                  %XXX
MedianFilterWindow = [1/5 12/5 15/5]*5;         % NO median filtering actually smoother for phase, and may be more finegrained/precise as well. Might not be accurate
PhaseZeroErrorDelta = pi/15.7;                  % 1/5 * pi
PhaseZeroError = 3.132;                         % To reduce some abrupt values and pi shifts, we constraint that person can't move about 7cm in just 1/300 seconds
ZeroErrorFactor = 15.7;                         %XXX
DopplerZeroError = 0.04*10;  %OG 0.04*10        % Kind of like subtracting mean of the noise, mean determined empirically.
PhaseBasedDopplerNoiseReductionFactor = 0.01; %OG .01  % If distance moved according to phase is at least this.
ConversionConstant = 39.3701;                   % meters to inches
TimeDiffThresh = 1000;                          % 1 second
PhaseDistanceMultiplier = 1/2;                  % 360 tag signal phase changes every lambda/2 radial tag movement.
DopplerDistanceMultiplier = 2/1;                % Because wave travels back and forth, so distance travelled to target back to
                                                %   radar/transmitter is 2R, distance moved 2*delta(R)

% the different fields of measurement the reader picked up
% radar/transmitter is 2R, distance moved 2*delta(R)
keySet =   {'RSSI', 'AbsPhase', 'Doppler', 'AntennaID', 'Wavelength','TagID','Timestamp'};
valueSet = 1:length(keySet);
mapObj = containers.Map(keySet,valueSet);
PhaseAngle = struct;
START = 1; END = -0;
%initialize length of the minimum phase output vector
minPhase_out = 0;
%get a unique array of tag ids to determine which tags to iterate over
tag_id_arr = unique(tag(:,6));
length_counter = 1;
%figure;
%loop through each element of the tag id array (so that each tag gets hit)
for idx = 1:numel(tag_id_arr)
    %assign actual tag number to the TagID
    TagID=tag_id_arr(idx);
    %loop through the number of antennas given
    for AntennaID = 1:ant_num
        %subselect just values that have right AntennaID
        Antenna_1 = tag(tag(:,mapObj('AntennaID'))==AntennaID,:);
        %get subselect of time stamps
        Timestamps_1 = Antenna_1(:,mapObj('Timestamp'));
        TimeStampsOriginal = Antenna_1(:,mapObj('Timestamp'));          % Keeping originals in separate vector
        %get initial time
        init_time=TimeStampsOriginal(1,1);
        %subtract initial time from all timestamps
        TimeStampsOriginal=TimeStampsOriginal-init_time;
        %correct time stamp units
        TimeStampsOriginal = TimeStampsOriginal(:,1)./1000000;
        
        % If we have multiple readers, then we need to tell reader ID as well.
        Antenna_1 = Antenna_1(Antenna_1(:,mapObj('TagID'))==TagID,:);
        % JH - TBD, not sure what this does
        Antenna_1 = Antenna_1(START:end+END,:);
        %get subselect of wavelengths
        Wavelengths = Antenna_1(:,mapObj('Wavelength'));
        %first order wavelength differential
        WavelengthDiff = diff(Wavelengths);
        %calculate frequencies in correct units
        Frequencies = 100*(3*10^2)./(Wavelengths);
        %first order frequency differential
        FrequenciesDiff = diff(Frequencies);
        Timestamps_1 = Timestamps_1(WavelengthDiff==0);
        TimestampsDiff_1 = diff(Timestamps_1)*0.000001; % seconds
        
        AbsPhase_1 = Antenna_1(:,mapObj('AbsPhase'));
        [CorrectPhase1] = PhaseCorrector(AbsPhase_1, WavelengthDiff,PhaseZeroError,PhaseZeroErrorDelta);
        [CorrectPhase2] = PhaseCorrector(flip(CorrectPhase1), flip(WavelengthDiff),PhaseZeroError,PhaseZeroErrorDelta);
        [PhaseAngle.(['A' num2str(AntennaID) 'T' num2str(TagID)])] = ChannelCorrector(flip(CorrectPhase2),WavelengthDiff,Wavelengths);
        PhaseAngle.(['A' num2str(AntennaID) 'T' num2str(TagID)]) = flip(CorrectPhase2);%CorrectPhase1;%flip(CorrectPhase2);
        AbsPhase_1 = PhaseAngle.(['A' num2str(AntennaID) 'T' num2str(TagID)]);
            %AbsPhase_1 = unwrap(AbsPhase_1,2*3.1); % In case motion above wavelength/2 occurs
        AbsPhase_1 = diff(AbsPhase_1); % Phase difference
        AbsPhase_1 = AbsPhase_1((WavelengthDiff==0));
        AbsPhase_1 = AbsPhase_1((TimestampsDiff_1 < TimeDiffThresh));
            %AbsPhase_1(AbsPhase_1 > PhaseZeroError) = AbsPhase_1(AbsPhase_1 > PhaseZeroError) - 3.14; % subtracting
            %AbsPhase_1(AbsPhase_1 < -PhaseZeroError) = AbsPhase_1(AbsPhase_1 < -PhaseZeroError) + 3.14; % adding
            %AbsPhase_1(AbsPhase_1>1.95*(PhaseZeroError)) = AbsPhase_1(AbsPhase_1>1.95*(PhaseZeroError)) - 2*PhaseZeroError;
            %AbsPhase_1(AbsPhase_1<-1.95*(PhaseZeroError)) = AbsPhase_1(AbsPhase_1<-1.95*(PhaseZeroError)) + 2*PhaseZeroError;
        AbsPhase_1 = PhaseDiffCorrector(AbsPhase_1,PhaseZeroError,ZeroErrorFactor);
        AP_1 = AbsPhase_1(AbsPhase_1 <= PhaseZeroError & AbsPhase_1 >= -PhaseZeroError);
        AP_1 = medfilt1(AP_1,MedianFilterWindow(1));
            %AbsPhase_1 = filter(ones(1,length(Antenna_1)),[1],AbsPhase_1);
        AP_1 = conv(AP_1, ones(1,MovingAvgWindow(1))/MovingAvgWindow(1), 'valid');
        
        Timestamps_1 = Antenna_1(:,mapObj('Timestamp'));
        Timestamps_1 = Timestamps_1(WavelengthDiff==0);
        Timestamps_1 = Timestamps_1(AbsPhase_1 <= PhaseZeroError & AbsPhase_1 >= -PhaseZeroError);
        TimestampsDiff_1 = diff(Timestamps_1)*0.000001 % seconds
        
        WavelengthsPruned_1 = Wavelengths(AbsPhase_1 <= PhaseZeroError & AbsPhase_1 >= -PhaseZeroError);
        DistanceMoved_1 = PhaseDistanceMultiplier*0.01*(WavelengthsPruned_1(1:end).*[zeros(MovingAvgWindow(1)-1,1); AP_1])/(2*pi); % meters
        DistanceMovedPerPhaseChange_1 = DistanceMoved_1;
            %DistanceMovedPerPhaseChange_1 = DistanceMovedPerPhaseChange_1(DistanceMovedPerPhaseChange_1 > 0.001 | DistanceMovedPerPhaseChange_1 < -0.001);
        DistanceMovedRunSum_1 = filter(1,[1 -1],abs(DistanceMovedPerPhaseChange_1(:))); %DistanceMoved > 0.0005
        DistanceMoved_1 = filter(1,[1 -1],DistanceMovedPerPhaseChange_1(:));
        
        DopplerVector = Antenna_1(:,mapObj('Doppler'));
        DopplerVector = DopplerVector(AbsPhase_1 <= PhaseZeroError & AbsPhase_1 >= -PhaseZeroError);
        DopplerVector = tsmovavg(DopplerVector,'w',gausswin(3*MedianFilterWindow(2))',1);
        DopplerVector = [DopplerVector(~isnan(DopplerVector)); zeros(sum(isnan(DopplerVector)),1)];
            %=DopplerVector = tsmovavg(DopplerVector,'e',20,1);
        DopplerVelocityChange_1 = 0.01*DopplerVector.*WavelengthsPruned_1/2; % meter/s
        DopplerVelocityChange_1 = conv(DopplerVelocityChange_1, ones(1,MovingAvgWindow(2))/MovingAvgWindow(2), 'valid');
        TEMP = [zeros(MovingAvgWindow(2)-1,1); (DopplerVelocityChange_1(:)-DopplerZeroError)];
        VfSquareMinusViSquare = (TEMP(2:end).^2) - (TEMP(1:end-1).^2)
        AccCalculation_1 = diff(TEMP)./TimestampsDiff_1;
        AccCalculation_1(isnan(AccCalculation_1)) = 0;
        DistanceCalculation_1 = DopplerDistanceMultiplier*VfSquareMinusViSquare./(2*AccCalculation_1);
        %JH - Kill off invalid distance calculations
        DistanceCalculation_1(isnan(DistanceCalculation_1)) = 0;
        %JH - kill off distance calculations that are outside noise bounds
        DistanceCalculation_1((DistanceMoved_1 < PhaseBasedDopplerNoiseReductionFactor & DistanceMoved_1 > -PhaseBasedDopplerNoiseReductionFactor)) = 0;
        %JH - 2-D median filtering, using the MedianFilterWindow
        DistanceCalculation_1 = medfilt1(DistanceCalculation_1,MedianFilterWindow(2));
            %DistanceCalculation_1 = DistanceCalculation_1(DistanceCalculation_1 > 0.0004 | DistanceCalculation_1 < -0.0004);
        DopplerDistanceMoved_1 = filter(1,[1 -1],DistanceCalculation_1);
        DopplerDistanceMoved_1RunSum_1 = filter(1,[1 -1],abs(DistanceCalculation_1));
        
        RSSI_1 = Antenna_1(:,mapObj('RSSI'));
        RSSI_1 = RSSI_1(AbsPhase_1 <= PhaseZeroError & AbsPhase_1 >= -PhaseZeroError);
        RSSI_1 = medfilt1(RSSI_1,MedianFilterWindow(3));
        RSSI_1 = conv(RSSI_1, ones(1,MovingAvgWindow(3))/MovingAvgWindow(3), 'valid');
        
            % Alligned plot Timestamps_1(2:end), DistanceMovedPerPhaseChange_1 = DistanceMoved_1
            % (ConversionConstant*DistanceMovedPerPhaseChange_1(DistanceMovedPerPhaseChange_1 > PhaseBasedDopplerNoiseReductionFactor | DistanceMovedPerPhaseChange_1 < -PhaseBasedDopplerNoiseReductionFactor))'
            %plot(diff(ConversionConstant*DistanceMoved_1));%ConversionConstant*DistanceMoved_1 % AP_1 - filtered phase, RSSI

        %% this section trims the phase output vectors to make sure they aren't the max size
        %copy the final phase-output values
        DistMoved_Out = diff(ConversionConstant*DistanceMoved_1);
        %store temp copy bc I am not sure how DistMoved_Out is used
        tempDistMoved_Out = DistMoved_Out;
        %get the temp candidate for minimum length
        tempminPhase_out = length(DistMoved_Out);
        %store the phase out length so I can debug
        MinPhaseOutLength(length_counter)=tempminPhase_out;
        %if this is the first time, I am the shortest, add me
        if idx==1 && AntennaID==1
            %this is the shortest (and only), make it the new minimum
            minPhase_out=tempminPhase_out;
            OutDist_Phase_mat=tempDistMoved_Out;
        %if this is the shortest phase vector to come through, trim others
        %in the matrix and add the current shortest
        elseif tempminPhase_out<minPhase_out
            %this is the shortest, make it the new minimum
            minPhase_out=tempminPhase_out;
            %delete the last portions of the matrix
            OutDist_Phase_mat((minPhase_out+1):end,:)=[];
            %add the latest dist calc to the final matrix
            OutDist_Phase_mat = [OutDist_Phase_mat tempDistMoved_Out];
        %if this is not the shortest, trim it to the shortest, add to the
        %matrix
        else
            %trim the temp to the length of the shortest (and size) of
            %matrix
            tempDistMoved_Out((minPhase_out+1):end,:)=[];
            %add the latest dist calc to the final matrix
            OutDist_Phase_mat = [OutDist_Phase_mat tempDistMoved_Out];
        end
        
            %plot(Timestamps_1, ConversionConstant*DistanceMoved_1)
    %         hold on;
    %         if DopplerOn
    %             plot((ConversionConstant*DistanceCalculation_1)) %Doppler velocity ./TimestampsDiff_1
    %             legend('Calculated Using Phase', 'Calculated Using Doppler')

            %plot((ConversionConstant*DistanceCalculation_1)) %Doppler velocity ./TimestampsDiff_1
            %plot(diff(ConversionConstant*DistanceCalculation_1), 'r') %Doppler velocity ./TimestampsDiff_1
            %% 
            %copy the final doppler-output values
            Dop_Out = ConversionConstant*DistanceCalculation_1;
            tempDop_Out = Dop_Out;
            tempminDop_out = length(Dop_Out);
            MinDopOutLength(length_counter)=tempminDop_out;
            %if this is the first time, I am the shortest, add me
            if idx==1 && AntennaID==1
                minDop_out=tempminDop_out;
                OutDist_Dop_mat=tempDop_Out;
            %if this is the shortest phase vector to come through, trim others
            %in the matrix and add the current shortest
            elseif tempminDop_out<minDop_out
                minDop_out=tempminDop_out;
                OutDist_Dop_mat((minDop_out+1):end,:)=[];
                OutDist_Dop_mat = [OutDist_Dop_mat tempDop_Out];
            %if this is not the shortest, trim it to the shortest, add to the
            %matrix
            else
                %shortest = size(OutDist_Phase_mat,1);
                tempDop_Out((minDop_out+1):end,:)=[];
                OutDist_Dop_mat = [OutDist_Dop_mat tempDop_Out];
            end
            length_counter=length_counter+1;
%             if combine==0
%                 legend('Phase', 'Doppler');%'First order diff Doppler')
%             end
            
% %         end

        
%         DistMoved_Out = diff(ConversionConstant*DistanceMoved_1);
%         OutDist_Phase=DistMoved_Out;
%         DopplerDist_Out = ConversionConstant*DistanceCalculation_1;
%         OutDist_Dop = DopplerDist_Out;
%         tempmaxDop = max(DopplerDist_Out);
%         tempminDop = min(DopplerDist_Out);
%         tempmaxPhase = max(DistMoved_Out);
%         tempminPhase = min(DistMoved_Out);
%         if tempmaxDop>maxDop
%             maxDop=tempmaxDop;
%         end
%         if tempmaxPhase>maxPhase
%             maxPhase=tempmaxPhase;
%         end
%         if tempminDop<minDop
%             minDop=tempminDop;
%         end
%         if tempminPhase>minPhase
%             minPhase=tempminPhase;
%         end
%         if plot_taps==1
%             hold on
%             for i=1:length(tap_arr)
%                 line([tap_arr(i) tap_arr(i)], [minPhase maxPhase]);
%                 line([tap_arr(i) tap_arr(i)], [minDop maxDop]);
%             end
%             hold off
%         end
%         line([tap1 tap1], [minPhase maxPhase]);
%         line([tap1 tap1], [minDop maxDop]);
%         line([tap2 tap2], [minPhase maxPhase]);
%         line([tap2 tap2], [minDop maxDop]);
%         line([tap3 tap3], [minPhase maxPhase]);
%         line([tap3 tap3], [minDop maxDop]);
%         line([tap4 tap4], [minPhase maxPhase]);
%         line([tap4 tap4], [minDop maxDop]);
%         line([tap5 tap5], [minPhase maxPhase]);
%         line([tap5 tap5], [minDop maxDop]);
%         %% 
%         
%         X = DistMoved_Out;
%         coeff = pca(X);
%         Itransformed = X*coeff;
%         figure; 
%         hold on
%         plot(Itransformed)
%         title('Wonky PCA');
%         hold off
        %% 


        %%plot(conv(AbsPhase_1, ones(1,window_size)/window_size, 'valid'))
        %%quarticMA = sgolayfilt(AbsPhase_1, 4, 7);
        
        %         subplot(4,1,1);plot(ConversionConstant*DistanceMovedRunSum_1);title('DistanceMovedRunSum - Phase & Doppler (inches)')
        %         hold on
        %         if DopplerOn
        %             plot(ConversionConstant*DopplerDistanceMoved_1RunSum_1);
        %             hold on
        %             if (AntennaID == 2)
        %                 legend('Ant-1-Using-Phase', 'Ant-1-Using-Doppler','Ant-2-Using-Phase', 'Ant-2-Using-Doppler','Location','northwest')
        %             end
        %         end
        %         if (AntennaID == 2)
        %             legend('Ant-1-Using-Phase','Ant-2-Using-Phase','Location','northwest')
        %         end
        %
        %         subplot(4,1,2);
        %         plot(ConversionConstant*DistanceMoved_1);title('DirectionalDistanceMoved - Phase & Doppler (inches)')
        %         hold on
        %         if DopplerOn
        %             plot(ConversionConstant*DopplerDistanceMoved_1);
        %             hold on
        %             if (AntennaID == 2)
        %                 legend('Ant-1-Using-Phase', 'Ant-1-Using-Doppler','Ant-2-Using-Phase', 'Ant-2-Using-Doppler','Location','northwest')
        %             end
        %         end
        %         if (AntennaID == 2)
        %             legend('Ant-1-Using-Phase','Ant-2-Using-Phase','Location','northwest')
        %         end
        %
        %         subplot(4,1,3);plot(ConversionConstant*DistanceMovedPerPhaseChange_1);title('DistanceMovedPerPhaseChange (inches)')
        %         hold on
        %         subplot(4,1,4);plot(ConversionConstant*DistanceCalculation_1);title('DopplerMovedPerDopplerChange (inches)')
        %         hold on
        
    end
end
% if combine==1
%         title('Radial-Dist, Tag Num: 1-3');
%         ylabel('inches'); xlabel('timestamps');
%          legend('Phase 1', 'Phase 2', 'Phase 3');
% end
% % maxDop = 0;
% % minDop = 0;
% maxPhase = .05;
% minPhase = -.05;
%  
%         line([tap1 tap1], [minPhase maxPhase]);
% %         line([tap1 tap1], [minDop maxDop]);
%         line([tap2 tap2], [minPhase maxPhase]);
% %         line([tap2 tap2], [minDop maxDop]);
%         line([tap3 tap3], [minPhase maxPhase]);
% %         line([tap3 tap3], [minDop maxDop]);
%         line([tap4 tap4], [minPhase maxPhase]);
% %         line([tap4 tap4], [minDop maxDop]);
%         line([tap5 tap5], [minPhase maxPhase]);
% %         line([tap5 tap5], [minDop maxDop]);
% R1 = ((ConversionConstant*DistanceMovedPerPhaseChange_1)./TimestampsDiff_1);
% R2 = ((ConversionConstant*DistanceCalculation_1)./TimestampsDiff_1);
% %% Angle of Arrival
% ANTENNA = struct;
% figure;
% for TagID = [999]
%     for AntennaID = 1
%         plot(PhaseAngle.(['A' num2str(AntennaID) 'T' num2str(TagID)]));
%         hold on
%     end
% end
%% our moving average filter has a delay.
% Any symmetric filter of length N will have a delay of (N-1)/2 samples. We can account for this delay manually.
end

