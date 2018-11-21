% CSE 825 RFID Project
% Jesse and Dan
% Main Script

%kkcontrol
%Columns
% Peak RSSI, Phase Angle (rad), RF Doppler Freq, Antenna Port Num, Wave
% length, ID num, Last Seen Time
%%

%cse825_plotRFIDinfo( 32.493908, 2, 'recreateSetup1ChipTouch.csv', 0);
%Data_kk = csvread('recreateControl2_0.csv');
%cse825_plotRFIDinfo( 32.493908, 2, 'Control', Data_kk);

% To do
% 1. Figure out our own (or understand Kamran's) preprocessing module(?)
% 2. Take a multi chip dataset (known chip location)
% 3. Take full test dataset
% 4. Manual analysis
%   A. Visual
%   B. Data mining techniques (PCA, modeling/SVM)
% 5. Automate as much analysis as possible
% 6. Write paper


addpath('..\KamRFID\rfidtouch\Matlab_Processing')
%% Playing with Kam's code
clear all
close all
clc
% [OutDist_Phase_mat, OutDist_Dop_mat, MinPhaseOutLength, MinDopOutLength] = KamFunc(filename, ant_num)
[phase_out, dop_out, phaselength, doplength] = KamFunc('258Pin_2Read.csv', 2);


mat1 = dop_out;

figure;
hold on
plot(mat1(:,1), 'g')
%plot(mat1(:,2), 'g')
plot(mat1(:,3), 'r')
%plot(mat1(:,4), 'r')
plot(mat1(:,5), 'b')
%plot(mat1(:,6), 'b')
title('Kam Filtered Dop - Antenna 1')
legend('Tag 2', 'Tag 5', 'Tag 8');
%legend('Tag 2 Ant 1', 'Tag 2 Ant 2', 'Tag 5 Ant 1', 'Tag 5 Ant 2', 'Tag 8 Ant 1', 'Tag 8 Ant 2');
hold off


figure;
hold on
plot(mat1(:,2), 'g')
%plot(mat1(:,2), 'g')
plot(mat1(:,4), 'r')
%plot(mat1(:,4), 'r')
plot(mat1(:,6), 'b')
%plot(mat1(:,6), 'b')
title('Kam Filtered Dop - Antenna 2')
legend('Tag 2', 'Tag 5', 'Tag 8');
hold off


mat1 = phase_out;

figure;
hold on
plot(mat1(:,1), 'g')
%plot(mat1(:,2), 'g')
plot(mat1(:,3), 'r')
%plot(mat1(:,4), 'r')
plot(mat1(:,5), 'b')
%plot(mat1(:,6), 'b')
title('Pre Filtered Phase - Antenna 1')
legend('Tag 2', 'Tag 5', 'Tag 8');
%legend('Tag 2 Ant 1', 'Tag 2 Ant 2', 'Tag 5 Ant 1', 'Tag 5 Ant 2', 'Tag 8 Ant 1', 'Tag 8 Ant 2');
hold off

figure;
hold on
plot(mat1(:,2), 'g')
%plot(mat1(:,2), 'g')
plot(mat1(:,4), 'r')
%plot(mat1(:,4), 'r')
plot(mat1(:,6), 'b')
%plot(mat1(:,6), 'b')
title('Kam Filtered Phase - Antenna 2')
legend('Tag 2', 'Tag 5', 'Tag 8');
hold off

%%

% outright spikes
% large changes (swing from neg to pos, large amplitude)
% other antenna reading?
% other tags?
% phys proximity

%find maxes
selector = phase_out(:,6);
%find nummax+nummin abs(biggest) components
nummax=50;
nummin=50;
[maxval, idx_max]=maxk(selector,nummax);
max_mat = [idx_max maxval];
[minval, idx_min]=mink(selector,nummin);
min_mat = [idx_min minval];
tokillidx=find(maxval<(mean(maxval)+std(maxval)));
max_mat(tokillidx,:)=[];
dummy_phase=zeros(length(selector));
dummy_phase(max_mat(:,1))=max_mat(:,2);

tokillidx_min=find(minval>(mean(minval)-std(minval)));
min_mat(tokillidx_min,:)=[];
dummy_phase(min_mat(:,1))=min_mat(:,2);

%dummy_phase(idx_min)=minval;

% dummy_phase=zeros(length(selector));
% dummy_phase(idx_max)=maxval;
% dummy_phase(idx_min)=minval;
figure;
hold on
plot(dummy_phase, 'b')
%plot(dummy_phase, 'b')
title('Antenna 1 Tag 2');
hold off

% tokillidx=find(maxval<(mean(maxval)+std(maxval)));
% max_mat(tokillidx,:)=[];

%%

windowSize = 8; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
mat1_f = filter(b,a,mat1(:,1));
mat2_f = filter(b,a,mat1(:,2));
mat3_f = filter(b,a,mat1(:,3));
mat4_f = filter(b,a,mat1(:,4));
mat5_f = filter(b,a,mat1(:,5));
mat6_f = filter(b,a,mat1(:,6));

figure;
hold on
plot(mat1_f, 'g')
plot(mat3_f, 'r')
plot(mat5_f, 'b')
title('Post Filtered Phase - Antenna 1')
legend('Tag 2', 'Tag 5', 'Tag 8');
hold off


figure;
hold on
plot(mat2_f, 'g')
plot(mat4_f, 'r')
plot(mat6_f, 'b')
title('Post Filtered Phase - Antenna 2')
legend('Tag 2', 'Tag 5', 'Tag 8');
hold off

%% 
%close all
Data_A = csvread('smallTapTest.csv');
%Data_A = csvread('kkcontrol.csv');
Data_kk = csvread('bigControl.csv');
%Data_A=Data_kk;
%Normalize time
init_time=Data_A(1,7);
phase_sev=Data_A(:,7)-init_time;
phase_sev = phase_sev(:,1)./1000000;
Data_A(:,7)=phase_sev;

linspace(0,5057,1);
%Data_A(Data_A(:,2)<.5)=6.2832;

%Find all wavelength's calculated
wave_arr = unique(Data_A(:, 5));
% iteratively plot every wavelength in given file
%thesize = size(wave_arr,1);
% iteratively plot the first 10 unique wavelengths in given file
% thesize=10;
for i=1:thesize
    cse825_plotRFIDinfo(wave_arr(i), 1, 'recreateSetup1ChipTouch.csv', Data_A);
end
%% 
%tic
%Normalize time
phase=Data_A;
% init_time=phase(1,7);
% phase_sev=phase(:,7)-init_time;
% phase_sev = phase_sev(:,1)./1000000;
% phase(:,7)=phase_sev;
%select just antenna 2 readings
phase = phase(phase(:,4)==1, :);



 
pca_in = zeros(50,length(phase(:,7)));
wave_arr = unique(phase(:, 5));
%toc
%temp_size = zeros(50,1);
for j=1:size(phase,1)
    wavelength_curr = phase(j,5);
    wave_vec_idx = find(wave_arr==wavelength_curr);
    %time_curr = phase(j,7);
    %propogate previous values forward
    if j>1
        pca_in(:,j) = pca_in(:,j-1);
    end
    %build pca_in, take current phase value and add to pca_in matrix
    pca_in(wave_vec_idx, j) = phase(j, 2);
end
%toc

figure;
hold on;
plot(phase(:,7) , pca_in(1,:), 'b+');
xlabel('Time')
ylabel('Phase Angle');
%title(['File: ', filename, ' Phase, wav:', num2str(wavelength), ' ant: ', num2str(ant)]);
hold off;




%% 
%select the minimum size from wavelengths
%transpose pca_in
pca_in=pca_in';
pca_in = [zeros(1,size(pca_in,2)); diff(pca_in)];
[Coeff, Score] = pca(pca_in);

downselect= pca_in*Coeff;

%% 
% figure;
% hold on;
% 
% for k=1:comp
%     plot
% end
% 
% hold off;


figure;
hold on;
plot(phase(:,7) , downselect(:,1), 'b+');
xlabel('Time')
ylabel('Phase Angle');
title('PCA Output 1st Component');

hold off;

figure;
hold on;
plot(phase(:,7) , downselect(:,1), 'b+');
plot(phase(:,7) , downselect(:,2), 'r+');
plot(phase(:,7) , downselect(:,3), 'k+');
plot(phase(:,7) , downselect(:,4), 'g+');
plot(phase(:,7) , downselect(:,5), 'c+');
xlabel('Time')
ylabel('Phase Angle');
title('PCA Output of Downselect');
legend('First', 'Second', 'Third', 'Fourth', 'Fifth');
hold off;


%% 
Length1 =size(Data_kk, 1);
time1 = linspace(0, 100, Length1);


Length =size(Data_A, 1);
time = linspace(0, 100, Length);

%Select only Antenna 1 data
Ant1=Data_A(Data_A(:,4)==1,:);
Ant1L =size(Ant1, 1);
Ant1T = linspace(0, 100, Ant1L);

%Select only Antenna 2 data
Ant2=Data_A(Data_A(:,4)==2,:);
Ant2L =size(Ant2, 1);
Ant2T = linspace(0, 100, Ant2L);
% 
% figure;
% hold on;
% plot(Ant1T, Ant1(:,1), 'r+', Ant2T, Ant2(:,1), 'b+');
% xlabel('Time (made up)')
% ylabel('RSSI Mag');
% legend('Ant 1', 'Ant 2');
% title('Raw Peak RSSI vs Time');
% hold off;
% 
% figure;
% hold on;
% plot(Ant1T, Ant1(:,2), 'r+', Ant2T, Ant2(:,2), 'b+');
% xlabel('Time (made up)')
% ylabel('Phase Angle');
% legend('Ant 1', 'Ant 2');
% title('Raw Phase Angle vs Time');
% hold off;

%filter data, apply a moving average filter of window=windowSize
windowSize = 50; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
Ant1_RSSI_F = filter(b,a,Ant1(:,1));
Ant1_P_F = filter(b,a,Ant1(:,2));
Ant2_RSSI_F = filter(b,a,Ant2(:,1));
Ant2_P_F = filter(b,a,Ant2(:,2));


%% Out of data/not currently useful
%Plot the filtered RSSI and Phase per Antenna
% figure;
% hold on;
% plot(Ant1T(1, 50:end), Ant1_RSSI_F(50:end, 1), 'r+', Ant2T(1, 50:end), Ant2_RSSI_F(50:end, 1), 'b+');
% xlabel('Time (made up)')
% ylabel('RSSI Mag');
% legend('Ant 1', 'Ant 2');
% title('Filtered Peak RSSI vs Time');
% hold off;

% figure;
% hold on;
% plot(Ant1T, Ant1_P_F, 'r+', Ant2T, Ant2_P_F, 'b+');
% xlabel('Time (made up)')
% ylabel('Phase Angle');
% legend('Ant 1', 'Ant 2');
% title('Filtered Phase Angle vs Time');
% hold off;
% 

%Plot the difference between antennas (filtered)
% L1 = size(Ant1_P_F, 1);
% L2 = size(Ant2_P_F, 1);
% L3 = abs(L2-L1)+1;
% if L1>L2
%     deltaPhase = abs(Ant2_P_F(:,1) - Ant1_P_F(L3:end,1));
%     AntPlt=Ant2T;
% elseif L2>L1
%     deltaPhase = abs(Ant2_P_F(L3:end,1) - Ant1_P_F(:,1));    
%     AntPlt=Ant1T;
% end
%     
% figure;
% hold on;
% plot(AntPlot, deltaPhase, 'r+');
% xlabel('Time (made up)')
% ylabel('Phase Angle');
% title('Diff Filtered Phase Angle vs Time');
% hold off;

%Normalize the RSSI
% degA = subA*(180/pi);
% RSSI_norm = abs(Data_A(:,1))./mean(abs(Data_A(:,1)));
% y_norm = y./max(y);

