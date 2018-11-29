function [out] = DanFunc(filename, ant_num)
% Clearing
% clc
% clear all
% close all

% Hard code loading the file and saving experiment number
out=struct;
dist_str = 'Dist';
tag_str = 'Time';
ant_str = 'Ant';
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

%Replace Timestamps with time since run began
Timestamp_replace = tag(:,mapObj('Timestamp'));
init_time= Timestamp_replace(1,1);
Timestamp_replace=Timestamp_replace-init_time;
Timestamp_replace = Timestamp_replace(:,1)./1000000;
tag(:,mapObj('Timestamp')) = Timestamp_replace;

%% loop through the number of antennas given
for AntennaID = 1:ant_num
    
    %% loop through each element of the tag id array (so that each tag gets hit)
    for idx = 1:numel(tag_id_arr)
        %assign actual tag number to the TagID
        TagID=tag_id_arr(idx);
        
        % get only data for current tag and antenna
        Data = tag(tag(:,mapObj('AntennaID'))==AntennaID,:);
        Data = Data(Data(:,mapObj('TagID'))==TagID,:);
        % select only ID, Antenna, Phase, and Timestamp
        ID_vec = Data(:,mapObj('TagID'));
        Antenna_vec = Data(:,mapObj('AntennaID'));
        PhaseAngle_vec = Data(:,mapObj('AbsPhase'));
        Timestamp_vec = Data(:,mapObj('Timestamp'));

        %get subselect of wavelengths
        Wavelength_vec = Data(:,mapObj('Wavelength'));
        %first order wavelength differential
        WavelengthDiff = diff(Wavelength_vec);
        Timestamps_1 = Timestamp_vec(WavelengthDiff==0);
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
        %DistMoved_Out = (ConversionConstant*DistanceMoved_1);
        TagIDstr_temp=num2str(TagID);
        AntennaIDstr_temp=num2str(AntennaID);
        Dist_str_temp = strcat(ant_str, AntennaIDstr_temp, dist_str, TagIDstr_temp);
        TagIDstr_temp = strcat(ant_str, AntennaIDstr_temp, tag_str, TagIDstr_temp);
        out.(TagIDstr_temp) = Timestamps_1;
        out.(Dist_str_temp) = DistMoved_Out;

    end
            
end

