function [out] = DanFunc(filename, ant_num)
% Clearing
% clc
% clear all
% close all

% Hard code loading the file and saving experiment number
out=struct;
dist_str = 'Dist';
tag_str = 'Time';
tag = load(filename);
save(['Exp' num2str(1) '.mat'],'tag')

%% Reading tag values
MovingAvgWindow = [10/5 15/5 15/5]*5;           %just an array of [10 15 15]. Don't know why there's scalar mult. Doesn't seem like types are diff
%DopplerOn = 0;                                  %XXX
MedianFilterWindow = [1/5 12/5 15/5]*5;         % NO median filtering actually smoother for phase, and may be more finegrained/precise as well. Might not be accurate
PhaseZeroErrorDelta = pi/15.7;                  % 1/5 * pi
PhaseZeroError = 3.14;                         % To reduce some abrupt values and pi shifts, we constraint that person can't move about 7cm in just 1/300 seconds
ZeroErrorFactor = 15.7;                         %XXX
DopplerZeroError = 0.04*10;  %OG 0.04*10        % Kind of like subtracting mean of the noise, mean determined empirically.
PhaseBasedDopplerNoiseReductionFactor = 0.01;   %OG .01  % If distance moved according to phase is at least this.
ConversionConstant = 39.3701;                   % meters to inches
TimeDiffThresh = 1000;                          % 1 second
PhaseDistanceMultiplier = 1/2;                  % 360 tag signal phase changes every lambda/2 radial tag movement.
DopplerDistanceMultiplier = 2/1;                % Because wave travels back and forth, so distance travelled to target back to
                                                %   radar/transmitter is 2R, distance moved 2*delta(R)
OutDist_Phase_mat = [];
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

Timestamp_vec = tag(:,mapObj('Timestamp'));
init_time=Timestamp_vec(1,1);
Timestamp_vec=Timestamp_vec-init_time;
Timestamp_vec_true = Timestamp_vec(:,1)./1000000;
display('hello')
%% loop through the number of antennas given
for AntennaID = 1:ant_num
    
     Data = tag(tag(:,mapObj('AntennaID'))==AntennaID,:);
    %% loop through each element of the tag id array (so that each tag gets hit)
    for idx = 1:numel(tag_id_arr)
        %assign actual tag number to the TagID
        TagID=tag_id_arr(idx);
        
        % get only data for current tag and antenna
        Data = tag;
        Data = Data(Data(:,mapObj('TagID'))==TagID,:);
        Timestamp_vec = Timestamp_vec_true;
        % select only ID, Antenna, Phase, and Timestamp
        ID_vec = Data(:,mapObj('TagID'));
        Antenna_vec = Data(:,mapObj('AntennaID'));
        PhaseAngle_vec = Data(:,mapObj('AbsPhase'));
%         Timestamp_vec = Data(:,mapObj('Timestamp'));
        %123321
        % Convert to seconds since runtime
%         init_time=Timestamp_vec(1,1);
%         Timestamp_vec=Timestamp_vec-init_time;
%         Timestamp_vec = Timestamp_vec(:,1)./1000000;

        %get subselect of wavelengths
        Wavelength_vec = Data(:,mapObj('Wavelength'));
        %first order wavelength differential
        WavelengthDiff = diff(Wavelength_vec);
        Timestamps_1 = Timestamp_vec(WavelengthDiff==0)
        Antenna_vec = Antenna_vec(WavelengthDiff==0);
        ID_vec=ID_vec(WavelengthDiff==0);
        TimestampsDiff_1 = diff(Timestamps_1)*0.000001; % seconds
        
        CorrectPhaseAngle_vec = PhaseCorrectorTrim(PhaseAngle_vec, WavelengthDiff,PhaseZeroError,PhaseZeroErrorDelta);
        CorrectPhaseAngle_vec = PhaseCorrectorTrim(flip(CorrectPhaseAngle_vec), flip(WavelengthDiff),PhaseZeroError,PhaseZeroErrorDelta);
        CorrectPhase2=CorrectPhaseAngle_vec;
        [PhaseAngle.(['A' num2str(AntennaID) 'T' num2str(TagID)])] = ChannelCorrector(flip(CorrectPhase2),WavelengthDiff,Wavelength_vec);
        PhaseAngle.(['A' num2str(AntennaID) 'T' num2str(TagID)]) = flip(CorrectPhase2);%CorrectPhase1;%flip(CorrectPhase2);
        AbsPhase_1 = PhaseAngle.(['A' num2str(AntennaID) 'T' num2str(TagID)]);
        AbsPhase_1 = diff(AbsPhase_1); % Phase difference
        AbsPhase_1 = AbsPhase_1((WavelengthDiff==0));
        AbsPhase_1 = AbsPhase_1((TimestampsDiff_1 < TimeDiffThresh));
        AbsPhase_1 = PhaseDiffCorrector(AbsPhase_1,PhaseZeroError,ZeroErrorFactor);
        AP_1 = AbsPhase_1(AbsPhase_1 <= PhaseZeroError & AbsPhase_1 >= -PhaseZeroError);
        Timestamps_1 = Timestamps_1(AbsPhase_1 <= PhaseZeroError & AbsPhase_1 >= -PhaseZeroError);
        WavelengthsPruned_1 = Wavelength_vec(AbsPhase_1 <= PhaseZeroError & AbsPhase_1 >= -PhaseZeroError);
        AP_1 = medfilt1(AP_1,MedianFilterWindow(1));
        AP_1 = conv(AP_1, ones(1,MovingAvgWindow(1))/MovingAvgWindow(1), 'valid');

        DistanceMoved_1 = PhaseDistanceMultiplier*0.01*(WavelengthsPruned_1(1:end).*[zeros(MovingAvgWindow(1)-1,1); AP_1])/(2*pi); % meters
        DistanceMovedPerPhaseChange_1 = DistanceMoved_1;
        DistanceMoved_1 = filter(1,[1 -1],DistanceMovedPerPhaseChange_1(:));

        %copy the final phase-output values
        DistMoved_Out = diff(ConversionConstant*DistanceMoved_1);
        TagIDstr_temp=num2str(TagID);
        Dist_str_temp = strcat(dist_str, TagIDstr_temp);
        TagIDstr_temp = strcat(tag_str, TagIDstr_temp);
        length(Timestamps_1)
        length(DistMoved_Out)
        out.(TagIDstr_temp) = Timestamps_1;
        out.(Dist_str_temp) = DistMoved_Out;
%         length(ID_vec);
%         length(Antenna_vec);
%         length(DistMoved_Out)
%         length(Timestamps_1)
%         
%         %% this section trims the phase output vectors to make sure they aren't the max size
%         %store temp copy bc I am not sure how DistMoved_Out is used
%         tempDistMoved_Out = DistMoved_Out;
%         %get the temp candidate for minimum length
%         tempminPhase_out = length(DistMoved_Out);
%         %if this is the first time, I am the shortest, add me
%         if idx==1 && AntennaID==1
%             %this is the shortest (and only), make it the new minimum
%             minPhase_out=tempminPhase_out;
%             OutDist_Phase_mat=tempDistMoved_Out;
%         %if this is the shortest phase vector to come through, trim others
%         %in the matrix and add the current shortest
%         elseif tempminPhase_out<minPhase_out
%             %this is the shortest, make it the new minimum
%             minPhase_out=tempminPhase_out;
%             %delete the last portions of the matrix
%             OutDist_Phase_mat((minPhase_out+1):end,:)=[];
%             %add the latest dist calc to the final matrix
%             OutDist_Phase_mat = [OutDist_Phase_mat tempDistMoved_Out];
%         %if this is not the shortest, trim it to the shortest, add to the
%         %matrix
%         else
%             %trim the temp to the length of the shortest (and size) of
%             %matrix
%             tempDistMoved_Out((minPhase_out+1):end,:)=[];
%             %add the latest dist calc to the final matrix
%             OutDist_Phase_mat = [OutDist_Phase_mat tempDistMoved_Out];
%         end
%         
%                 %% this section trims the phase output vectors to make sure they aren't the max size
%         %store temp copy bc I am not sure how DistMoved_Out is used
%         tempTime = Timestamps_1;
%         %get the temp candidate for minimum length
%         tempminTime = length(Timestamps_1);
%         %if this is the first time, I am the shortest, add me
%         if idx==1 && AntennaID==1
%             %this is the shortest (and only), make it the new minimum
%             minTime_out=tempminTime;
%             Out_Time_mat=tempTime;
%         %if this is the shortest phase vector to come through, trim others
%         %in the matrix and add the current shortest
%         elseif tempminTime<minTime_out
%             %this is the shortest, make it the new minimum
%             minTime_out=tempminTime;
%             %delete the last portions of the matrix
%             Out_Time_mat((minTime_out+1):end,:)=[];
%             %add the latest dist calc to the final matrix
%             Out_Time_mat = [Out_Time_mat tempTime];
%         %if this is not the shortest, trim it to the shortest, add to the
%         %matrix
%         else
%             %trim the temp to the length of the shortest (and size) of
%             %matrix
%             tempTime((minTime_out+1):end,:)=[];
%             %add the latest dist calc to the final matrix
%             Out_Time_mat = [Out_Time_mat tempTime];
%         end
        
%         SubData = [ID_vec Antenna_vec DistMoved_Out Timestamps_1];
        %SubData = [ID_vec Antenna_vec CorrectPhaseAngle_vec Timestamp_vec];
        
        
%         OutDist_Phase_mat = [OutDist_Phase_mat; SubData];
    end
            
end

