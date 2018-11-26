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



%% Playing with Kam's code
clear all
%close all
clc

file = 'fullPin_1.csv';
%call the function to calc radial dist
[output1] = DanFunc(file, 1);
%delete last timestamp due to the diff on dist
color_arr = [1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 .5; 0.8500 0.3250 0.0980; 0 .5 0];
plotRFIDTry2(output1, color_arr)
legend('Tag 1', 'Tag 2', 'Tag 3', 'Tag 4', 'Tag 5' , 'Tag 6', 'Tag 8', 'Tag 9');


mat1=phase_out;
figure;
hold on
plot(time_out(:,1), mat1(:,1), 'r')
plot(time_out(:,2), mat1(:,2), 'b')
plot(time_out(:,3), mat1(:,3), 'g')
plot(time_out(:,4), mat1(:,4), 'k')
plot(time_out(:,5), mat1(:,5), 'c')
plot(time_out(:,6), mat1(:,6), 'm')
plot(time_out(:,7), mat1(:,7), 'Color', [0.8500, 0.3250, 0.0980])
plot(time_out(:,8), mat1(:,8), 'Color', [0, 0.5, 0])

title('Radial Dist - All Tags')
legend('Tag 1', 'Tag 2', 'Tag 3', 'Tag 4', 'Tag 5' , 'Tag 6', 'Tag 8', 'Tag 9');
hold off

%% Plot the Radial dist values
%make a copy to avoid messing with output
mat1=phase_out;
figure;
hold on
plot(time_out(:,1), mat1(:,1), 'r')
title('Radial Dist - Tag 1')
hold off

figure;
hold on
plot(time_out(:,2), mat1(:,2), 'b')
title('Radial Dist - Tag 2')
hold off

figure;
hold on
plot(time_out(:,3), mat1(:,3), 'g')
title('Radial Dist - Tag 3')
hold off


figure;
hold on
plot(time_out(:,1), mat1(:,1), 'r')
plot(time_out(:,2), mat1(:,2), 'b')
plot(time_out(:,3), mat1(:,3), 'g')
title('Radial Dist - All Tags')
legend('Tag 1', 'Tag 2', 'Tag 3');
hold off


%% Perform keystroke extraction
out1 = keystrokeExtraction(phase_out, time_out, 70, 70, 1);
out2 = keystrokeExtraction(phase_out, time_out, 70, 70, 2);
out3 = keystrokeExtraction(phase_out, time_out, 70, 70, 3);
figure;
hold on
plot(time_out(:,1), out1, 'r')
plot(time_out(:,2), out2, 'b')
plot(time_out(:,3), out3, 'g')
legend('Tag 1', 'Tag 2', 'Tag 3');
title('Antenna 1');
hold off

%% Tester
data = phase_out;
nummax = 70;
nummin = 70;
col = 1;
%find maxes
selector = data(:,col);
length(selector)
%find nummax+nummin abs(biggest) components
[maxval, idx_max]=maxk(selector,nummax);
max_mat = [idx_max maxval];
[minval, idx_min]=mink(selector,nummin);
min_mat = [idx_min minval];
tokillidx_max=find(maxval<(mean(maxval)+std(maxval)));
tokillidx_min=find(minval>(mean(minval)-std(minval)));
min_mat(tokillidx_min,:)=[];
max_mat(tokillidx_max,:)=[];

dummy_phase=zeros(length(selector),1);
dummy_phase(max_mat(:,1))=max_mat(:,2);
dummy_phase(min_mat(:,1))=min_mat(:,2);

dummy_phase1= dummy_phase(dummy_phase~=0);
figure; plot(dummy_phase1)
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
