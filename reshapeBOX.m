function OUT = reshapeBOX(IN,subset)
%Reshapes cell array for box plot function.
%
%IN must be a cell array. Subset is the subset of arrays that should be
%included in the reshaping.
%
%OUT is a (:,2) double, where the second column is the column of the
%cell array, starting at 1.

numb = num2cell(1:size(IN(:),1));

out = cellfun(@(x,y) [x(:) y*ones(size(x(:)))],IN,numb,'uniformoutput',0);
OUT = vertcat(out{subset});
