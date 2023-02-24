%Converts genenames to BSU names
function [operons_bsu] = gene2bsu_op(x)
%[~, ~, readin]= xlsread('\\AGMAIERS2\Melih\NGS\Gene_names.xlsx','1');
%[~, ~, readin]= xlsread('\\AGMAIERS2\Melih\NGS\Gene_names_2020.xlsx','1');
[~, ~, readin]= xlsread('C:\Users\ESN\Documents\MATLAB\Matlab\EvoW23\Gene_names_2020.xlsx','1'); 

%BSU Numbers to gene name
operons_bsu = cell(size(x,1),size(x,2));
for n=1:size(x,2)
    for i= 1:length(x)
        for z= 1:size(readin,2)
        index = find(strcmp(x(i,n),readin(:,z)));
            if not(isempty(index))
             operons_bsu(i,n) = readin(index(1),1);
            end
        end
    end
end
