function [OutDist_Phase_mat] = DanFunc(filename, ant_num)
% Clearing
% clc
% clear all
% close all

% Hard code loading the file and saving experiment number
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
        %123321
        % Convert to seconds since runtime
        init_time=Timestamp_vec(1,1);
        Timestamp_vec=Timestamp_vec-init_time;
        Timestamp_vec = Timestamp_vec(:,1)./1000000;

        %get subselect of wavelengths
        Wavelength_vec = Data(:,mapObj('Wavelength'));
        %first order wavelength differential
        WavelengthDiff = diff(Wavelength_vec);
        Timestamp_1 = Timestamp_vec(WavelengthDiff==0);
        TimestampsDiff_1 = diff(Timestamp_1)*0.000001; % seconds
        
        CorrectPhaseAngle_vec = PhaseCorrectorTrim(PhaseAngle_vec, WavelengthDiff,PhaseZeroError,PhaseZeroErrorDelta);
        CorrectPhaseAngle_vec = PhaseCorrectorTrim(flip(CorrectPhaseAngle_vec), flip(WavelengthDiff),PhaseZeroError,PhaseZeroErrorDelta);
        
        SubData = [ID_vec Antenna_vec CorrectPhaseAngle_vec Timestamp_vec];
        
        
        OutDist_Phase_mat = [OutDist_Phase_mat; SubData];
    end
            
end

