function [  ] = cse825_plotRFIDinfo( wavelength, ant, filename, Dataset)
%pass in a 0 for wavelength, ant, or dataset to ignore variable

%if a dataset is passed in, use it. If not, read the file
if size(Dataset,1)>1
    Data=Dataset;
    filename='Passed In';
else
    Data = csvread(filename);
end
Data_c=Data;
%select subset of data, all rows of wavelength=wavelength
if wavelength==0
    display('No wavelength given, plotting all wavelengths');
else
    Data = Data(Data(:,5)==wavelength, :);
end

if ant==0
    display('No antenna given, plotting all antennas');
else
    Data = Data(Data(:,4)==ant, :);
end

Length1 =size(Data, 1)
time1 = linspace(0, 100, Length1);
% Length1 =size(Data_c, 1)
% time1 = linspace(0, 100, Length1);

figure;
hold on;
%plot(time1, Data(:,2), 'b+');
plot(Data(:,7), Data(:,2), 'b+');
xlabel('Time')
ylabel('Phase Angle');
title(['File: ', filename, ' Phase, wav:', num2str(wavelength), ' ant: ', num2str(ant)]);
hold off;
% 
% figure;
% hold on;
% plot(time1, Data(:,1), 'b+');
% xlabel('Time (made up)')
% ylabel('RSSI');
% title(['File: ', filename, ' RSSI, wav:', num2str(wavelength), ' ant: ', num2str(ant)]);
% hold off;


end

