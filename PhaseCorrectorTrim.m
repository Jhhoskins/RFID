% Params:
%   vector InputPhase -     The phase angles recorded
%   vector InputWavelengthDiff -    The change in wavelength measured from read to read
%                                   Stores the locations where the RFID reader changed wavelengths
%   double PhaseZeroError - 
%   double PhaseZeroErrorDelta - 
% returns:
%   vector CorrectPhase - 
function [CorrectPhase] = PhaseCorrectorTrim(InputPhase,InputWavelengthDiff,PhaseZeroError,PhaseZeroErrorDelta)

PiSubtraction = 3.135; % Value of pi to be subtracted.
SP = InputPhase(1); % The last value seen in the array
%WP_1 = find(InputWavelengthDiff~=0); % find locations of all non-zero elements
WP = [SP; InputWavelengthDiff];
TEMPPhaseAngle = zeros(length(InputPhase),1);   %Preallocate array for time savings
baseline = 0;

% Sometimes the phase jumps pi radians, correct that
for P = 1:length(InputPhase)
    % RFID Reader has switched to a new wavelength. Don't need to do anything this time
    if (WP(P))
        baseline = InputPhase(P);
        TEMPPhaseAngle(P) = InputPhase(P) - baseline;
        % When a new wavelength comes, the phase value jumps. Account for that
        % by returning to baseline
        SP = InputPhase(P);
    % RFID reader is on the same wavelength as last time
    else
        
        phaseDiff = InputPhase(P) - SP;     % difference between last value seen and current value
        
        
        % if phaseDiff is near pi, subtract pi from value
        if (abs(phaseDiff - PhaseZeroError) < PhaseZeroErrorDelta)
            TEMPPhaseAngle(P) =  InputPhase(P) - PiSubtraction;
            
        % if phaseDiff is near -pi, add pi to value
        elseif (abs(phaseDiff + PhaseZeroError) < PhaseZeroErrorDelta)
            TEMPPhaseAngle(P) = InputPhase(P) + PiSubtraction;
        
        %2pi
        elseif (abs(phaseDiff - 2*PhaseZeroError) < PhaseZeroErrorDelta)
            TEMPPhaseAngle(P) =  InputPhase(P) - 2*PiSubtraction;
        
        elseif (abs(phaseDiff + 2*PhaseZeroError) < PhaseZeroErrorDelta)
            TEMPPhaseAngle(P) =  InputPhase(P) + 2*PiSubtraction;
        
        % Who knows, but is rarely called if ever
        elseif (SP > InputPhase(P))...      % Taking into account the 2*pi wraps.
                && (abs(InputPhase(P)-mod(SP+PiSubtraction,PiSubtraction*1.9))...
                    < PhaseZeroErrorDelta+0.1) 
            TEMPPhaseAngle(P) = mod(InputPhase(P)-PiSubtraction,PiSubtraction*1.9);%+3.14]; % Just decrease previous one.
        
        % Otherwise there wasn't a spike, record normally
        else
            TEMPPhaseAngle(P) = InputPhase(P);
        end

        % We saw a new value, SP = this value
        SP = TEMPPhaseAngle(P);
        TEMPPhaseAngle(P) = TEMPPhaseAngle(P) - baseline;
    end
    
    
    
end
% TEMPPhaseAngle(TEMPPhaseAngle>PiSubtraction) = TEMPPhaseAngle(TEMPPhaseAngle>PiSubtraction)-PiSubtraction;
% Subtracting pi, as PLL reports (Phase + pi), why not just subtract pi :),
% % I will help in AoA measurements.
% Might need to add pi and take mod ?? == > Exactly similar ? YES
% But it will not be accurate in cases where there are readings with 2pi
% wrap. E.g. 5.1416 and then 2 (5.1416+pi mod 2*pi). Don't know which to chose
% TEMPPhaseAngle(TEMPPhaseAngle>PiSubtraction) = mod(TEMPPhaseAngle(TEMPPhaseAngle>PiSubtraction)+PiSubtraction,2*PiSubtraction);
% TEMPPhaseAngle(TEMPPhaseAngle<-PiSubtraction) = mod(TEMPPhaseAngle(TEMPPhaseAngle<-PiSubtraction)-PiSubtraction,-2*PiSubtraction);
%wpDiff = diff(WP_1);
CorrectPhase = TEMPPhaseAngle;
% figure;
% hold on
% plot(CorrectPhase)
% title('Tag 2')
% hold off
end