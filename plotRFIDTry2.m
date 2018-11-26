function [] = plotRFIDTry2(out_struct, color_arr)

size=length(fieldnames(out_struct));
out_cell= struct2cell(out_struct);
figure;
hold on;
color_count=1;
for i=1:2:size
    
%     time = cell2mat(out_cell(i));
%     time(end)=[]
%     val=i+1;
%     dist = cell2mat(out_cell(val));

    time = out_cell{i};
    time(end)=[]
    val=i+1;
    dist = out_cell{val};
    plot(time,dist, 'Color', color_arr(color_count,:))
    color_count = color_count+1;
    
end

end

