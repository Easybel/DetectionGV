% What the script does:
% Loads the operons list (March 2020) --> operons, converts genes to BSU names with gene2bsu function -->
% operons_bsu, with bed file it find the StartBSU (first BSU in operon), LastBSU (last BSU in operon),
% StartPos and EndPos of the operon --> operons_pos, sorts by start position --> operons_pos_sorted.

% Input: 
% files : operons.csv, Gene_names_2020.xlsx, BsubNC_000964wt_newRef.bed
% functions: gene2bsu_op.m Load_operons.m

clear all; close all

Load_operons 
operons_bsu=gene2bsu_op(operons); % Converting gene names in operons to BSU
recipbed = 'C:\Users\ESN\Documents\MATLAB\Matlab\EvoW23\BsubNC_000964wt_newRef.bed.txt';

% v is a vector with number of genes in operons
v = zeros(length(operons_bsu),1);
for i = 1:size(operons_bsu,1)
    x = find(~cellfun('isempty', operons_bsu(i,:)));
    v(i,1) = length (x);
end
clear x i 


% operons_large = find (vv(:) >= 5);


%Look at operons with >=2 genes
no = 2;
idx = v(:) >= no;
D = operons_bsu(idx,:);
clear idx operons_bsu no v 
operons_bsu = D;

%Load gene names, BSU and positions from bed file 
fid = fopen(recipbed); bed = textscan(fid,'%s %f %f %f %f %s','delimiter',' ');
fclose(fid);
gene168.GN=bed{1}; gene168.S=bed{2}; gene168.E=bed{3}; gene168.L=bed{4};
gene168.BSU=bed{6}; clear bed fid ans


operons_pos = struct('StartBSU',[],'LastBSU',[],'StartPos',[],'EndPos',[]);

%Finds first gene from operon in gene168 structure, if not (unknown gene
%like new_...) then second
for i=1:length(operons_bsu)
    idx = find(strcmp(operons_bsu(i,1),gene168.BSU));
    if not(isempty(idx))
        x = gene168.S(idx);
    else
        idx = find(strcmp(operons_bsu(i,v+1),gene168.BSU));
        x = gene168.S(idx);
    end
    
    % Finds last gene, if not last - 1 gene
    v = length(find(~cellfun('isempty', operons_bsu(i,:))));
    idxx = find(strcmp(operons_bsu(i,v),gene168.BSU));
    if not(isempty(idxx))
        y = gene168.S(idxx);
    else
        idxx = find(strcmp(operons_bsu(i,v-1),gene168.BSU));
        y = gene168.S(idxx);
    end
    
    % If Position od first gene is lower the order stays if not the order
    % is flipped 
    if x < y
        operons_pos(i).StartBSU = gene168.BSU(idx);
        operons_pos(i).LastBSU = gene168.BSU(idxx);
        operons_pos(i).StartPos = gene168.S(idx);
        operons_pos(i).EndPos = gene168.E(idxx);
    else
        operons_pos(i).StartBSU = gene168.BSU(idxx);
        operons_pos(i).LastBSU = gene168.BSU(idx);
        operons_pos(i).StartPos = gene168.S(idxx);
        operons_pos(i).EndPos = gene168.E(idx);
    end
end
        
clear i idx idxx v x y 

%Clears all empty rows
x = operons_pos(arrayfun(@(operonpos1) ~isempty(operonpos1.StartBSU),operons_pos));
y = x(arrayfun(@(D) ~isempty(D.LastBSU),x));
%clear operons_pos;
operons_pos = y;
clear x y 


%Sort the structure
T = struct2table(operons_pos); 
T_sorted = sortrows(T, 'StartBSU'); 
operons_pos_sorted = (table2struct(T_sorted)).';
clear T T_sorted
