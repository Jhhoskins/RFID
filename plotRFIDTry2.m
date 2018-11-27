function [out_cell] = plotRFIDTry2(out_struct, color_arr)

size=length(fieldnames(out_struct));
out_cell= struct2cell(out_struct);
figure;
hold on;
color_count=1;
%iterate i being every other
for i=1:2:size

    % delete the last time value for each part
    out_cell{i}(end)=[];
    time = out_cell{i};
    dist = out_cell{i+1};
    plot(time,dist, 'Color', color_arr(color_count,:));
    color_count = color_count+1;
end

end

