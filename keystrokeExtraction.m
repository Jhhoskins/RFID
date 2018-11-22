function [dummy_phase] = keystrokeExtraction(data, nummax, nummin, col)
% data is the data input

%find maxes
selector = data(:,col);
%find nummax+nummin abs(biggest) components
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

end

