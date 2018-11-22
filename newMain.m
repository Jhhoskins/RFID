% CSE 825 RFID Project
% Jesse and Dan
% Main Script

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
%close all
clc

file = '123Pin_1.csv';
%file = '../tagTap/test/123Pin_0.csv';

% [OutDist_Phase_mat, OutDist_Dop_mat, MinPhaseOutLength, MinDopOutLength] = KamFunc(filename, ant_num)
%[phase_out, dop_out, phaselength, doplength] = KamFunc('../tagTap/test/123Pin_control.csv', 2);
%[phase_out, dop_out, phaselength, doplength] = KamFunc('../tagTap/test/123Pin_1.csv', 1);
%[phase_out1, dop_out, phaselength1, doplength] = KamFunc(file, 1);
%[phase_out2, phaselength2] = KamFuncTrim(file, 1);
[phase_out] = DanFunc(file, 1);
mat1=phase_out;
figure;
hold on
plot(mat1(:,1), 'r')
title('Pre Filtered Phase - Tag 1 Antenna 1')
hold off

figure;
hold on
plot(mat1(:,2), 'b')
title('Pre Filtered Phase - Tag 2 Antenna 1')
hold off

figure;
hold on
plot(mat1(:,3), 'g')
title('Pre Filtered Phase - Tag 3 Antenna 1')
hold off

% error check, lots of values are about x*pi away from line. Fix them
%subtract pi from values in array
% piVal = 2.5;
% for i = 1:length(phase_out)
%     phaseDiff = phase_out(i,3);
%     %above
%     if (phase_out(i,3) > 2*piVal)
%         phase_out(i,3) = phase_out(i,3) - 2*pi;
%     elseif (phase_out(i,3) > piVal)
%         phase_out(i,3) = phase_out(i,3) - pi;
%     %below
%     elseif (phase_out(i,3) < -2*piVal)
%         phase_out(i,3) = phase_out(i,3) + 2*pi;
%     elseif (phase_out(i,3) < -piVal)
%         phase_out(i,3) = phase_out(i,3) + pi;
%     end
% end




mat1 = phase_out(phase_out(:,1)==1,:);
mat2 = phase_out(phase_out(:,1)==2,:);
mat3 = phase_out(phase_out(:,1)==3,:);



% windowSize = 5;
% b = (1/windowSize)*(ones(1,windowSize));
% a = 5;
% 
% mat1 = filter(b,a,mat1);
% mat2 = filter(b,a,mat2);
% mat3 = filter(b,a,mat3);


figure;
hold on
plot(mat1(:,4), mat1(:,3), 'r')
plot(mat2(:,4), mat2(:,3), 'g')
plot(mat3(:,4), mat3(:,3), 'b')
%plot(mat1(:,4), 'r')
%plot(mat1(:,5), 'b')
%plot(mat1(:,6), 'b')
title('Plot by Timestamp')
legend('Tag 1', 'Tag 2', 'Tag 3');
%legend('Tag 2 Ant 1', 'Tag 2 Ant 2', 'Tag 5 Ant 1', 'Tag 5 Ant 2', 'Tag 8 Ant 1', 'Tag 8 Ant 2');
hold off

mat11 = diff(mat1);
mat22 = diff(mat2);
mat33 = diff(mat3);

figure;
hold on
plot(mat1(2:end,4), mat11(:,3), 'r')
plot(mat2(2:end,4), mat22(:,3), 'g')
plot(mat3(2:end,4), mat33(:,3), 'b')
title('Plot by Timestamp - Diff')
legend('Tag 1', 'Tag 2', 'Tag 3');
%legend('Tag 2 Ant 1', 'Tag 2 Ant 2', 'Tag 5 Ant 1', 'Tag 5 Ant 2', 'Tag 8 Ant 1', 'Tag 8 Ant 2');
hold off



%every time the wavelength shifts,
%subtract new baseline from following
%values
 
figure;
hold on
plot(mat1(:,4), mat1(:,3), 'r')
title('Tag 1')
hold off

figure;
hold on
plot(mat2(:,4), mat2(:,3), 'g')
title('Tag 2')
hold off

figure;
hold on
plot(mat3(:,4), mat3(:,3), 'b')
title('Tag 3')
hold off




%% Sidenotes
%%% plotting phase over time
% col1 = TagNum, col2 = Value, col3 = Timestamp
% vec = [1 1 0;
%        1 2 2; 
%        1 1 6;
%        1 1 9; 
%        1 3 13; 
%        1 2 15;
%        2 1 3;
%        2 2 4;
%        2 3 8;
%        2 2 16;
%        2 1 17;
%        2 1 18;
%        2 1 19; 
%        2 2 23;
%        2 1 26;
%        3 3 7;
%        3 3 10;
%        3 1 11;
%        3 3 12;
%        3 2 14];

% a = vec(vec(:,1) == 1,:);
% b = vec(vec(:,1) == 2,:);
% c = vec(vec(:,1) == 3,:);

% figure;
% hold on
% plot(a(:,3),a(:,2))
% plot(b(:,3),b(:,2))
% plot(c(:,3),c(:,2))
% xlim([0 30])
% ylim([0 4])
