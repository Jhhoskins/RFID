function [CorrectPhase] = PhaseCorrector(InputPhase,InputWavelengthDiff,PhaseZeroError,PhaseZeroErrorDelta)

PiSubtraction = 3.135; % Value of pi to be subtracted.
PhaseFactor = 1.85; % Multiplier for checking a 2*pi phase diff., 1.95
SP = InputPhase(1); % Input raw phase matrix.
WP = find(InputWavelengthDiff~=0); % Input raw wavelength diff matrix.
TEMPPhaseAngle = [];
WPLatest = 1;
for P = 1:length(InputPhase)
    if (~isempty(find((WP+1)==P))) % The New Wavelength, first index. Don't need to do anything with that.
        TEMPPhaseAngle = [TEMPPhaseAngle; InputPhase(P)];
        SP = TEMPPhaseAngle(P);
        WPLatest = P;
    else
        if (InputPhase(P) - SP > (PhaseZeroError-PhaseZeroErrorDelta)) && (InputPhase(P) - SP < (PhaseZeroError+PhaseZeroErrorDelta))
            TEMPPhaseAngle = [TEMPPhaseAngle; InputPhase(P) - PiSubtraction];
            SP = TEMPPhaseAngle(P);
        elseif (InputPhase(P) - SP < (-PhaseZeroError+PhaseZeroErrorDelta)) && (InputPhase(P) - SP > -(PhaseZeroError+PhaseZeroErrorDelta))
            TEMPPhaseAngle = [TEMPPhaseAngle; InputPhase(P)];%+3.14]; % Just decrease previous one.
            TEMPPhaseAngle(P-1) = TEMPPhaseAngle(P-1) - PiSubtraction;
            SP = TEMPPhaseAngle(P);
        elseif (SP > InputPhase(P)) && (abs(InputPhase(P)-mod(SP+PiSubtraction,PiSubtraction*1.9)) < PhaseZeroErrorDelta+0.1) % Taking into account the 2*pi wraps.
            TEMPPhaseAngle = [TEMPPhaseAngle; mod(InputPhase(P)-PiSubtraction,PiSubtraction*1.9)];%+3.14]; % Just decrease previous one.
            SP = TEMPPhaseAngle(P);
        else
            TEMPPhaseAngle = [TEMPPhaseAngle; InputPhase(P)];
            SP = TEMPPhaseAngle(P);
        end
    end
    
    % For the consecutive values which turn out to be very close to 2pi
    % wrap. A little phase change and the jump occurs. This is to avoid
    % those jumps.
    if P > 1 && WPLatest ~= P % Don't want to do this for the first value of phase when channel changes.
        if TEMPPhaseAngle(P-1)-SP >= PhaseFactor*(PhaseZeroError)
            TEMPPhaseAngle(P-1) = TEMPPhaseAngle(P-1) - 1.9*PiSubtraction;
        elseif TEMPPhaseAngle(P-1)-SP <= -PhaseFactor*(PhaseZeroError)
            TEMPPhaseAngle(P) = TEMPPhaseAngle(P) - 1.9*PiSubtraction;
            SP = TEMPPhaseAngle(P);
        end
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
CorrectPhase = TEMPPhaseAngle;
end