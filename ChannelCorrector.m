function [CorrectPhase] = ChannelCorrector(InputPhase,InputWavelengthDiff,Wavelengths)

Frequencies = 100*(3*10^2)./(Wavelengths);
WP = find(InputWavelengthDiff~=0); % Input raw wavelength diff matrix.
TEMPPhaseAngle = [];
for P = 1:length(InputPhase);
    if ~isempty(find((WP+1)==P)) % The New Wavelength, first index. Don't need to do anything with that.
        SP = InputPhase(P); % -(4*pi*100*(Wavelengths(P-1)-Wavelengths(P))/(Wavelengths(P-1)*Wavelengths(P)));
        % SP = SP - 0.142*(round(Wavelengths(P-1)-Wavelengths(P),2))/0.02;
        TEMPPhaseAngle = [TEMPPhaseAngle; SP];
    else
        TEMPPhaseAngle = [TEMPPhaseAngle; InputPhase(P)];
    end
end
CorrectPhase = TEMPPhaseAngle;

end