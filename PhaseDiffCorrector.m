function [CorrectPhaseDiff] = PhaseDiffCorrector(PhaseDiff,PhaseZeroError,ZeroErrorFactor)

TEMPPhaseDiff = zeros(length(PhaseDiff),1);
SP = 0;
for P = 1:length(PhaseDiff)
    if PhaseDiff(P) <= PhaseZeroError/ZeroErrorFactor && PhaseDiff(P) >= -PhaseZeroError/ZeroErrorFactor
        TEMPPhaseDiff(P) = PhaseDiff(P);
        SP = PhaseDiff(P);
    else
        TEMPPhaseDiff(P) = SP;
    end
end
CorrectPhaseDiff = TEMPPhaseDiff;
end