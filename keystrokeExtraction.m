function [dummy_phase] = keystrokeExtraction(data, time, nummax, nummin, col)
% data is the data input

%find maxes
selector = data(:,col);
time_sel = time(:,col);
%find nummax+nummin abs(biggest) components
[maxval, idx_max]=maxk(selector,nummax);
max_mat = [idx_max maxval];
[minval, idx_min]=mink(selector,nummin);
min_mat = [idx_min minval];
tokillidx_max=find(maxval<(mean(maxval)+std(maxval)));
tokillidx_min=find(minval>(mean(minval)-std(minval)));
min_mat(tokillidx_min,:)=[]
max_mat(tokillidx_max,:)=[]

dummy_phase=zeros(length(selector),1);
dummy_phase(max_mat(:,1))=max_mat(:,2);
dummy_phase(min_mat(:,1))=min_mat(:,2);
        
end

