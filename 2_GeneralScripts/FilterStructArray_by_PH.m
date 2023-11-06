%
%
%%%%%%%%%%%%%%%%   FilterStructArray created by Philipp  %%%%%%%%%%%%%%%%%%%%%%%%
%
%
% What this script does: It loads a struct saved in .mat file and
% filters/removes each row for given criteria. Can save to a new .mat file.
% Example usage: filter AccMM2Genes.mat for A2_List2Genes
%
% Also, one could use it for existing workspace variables by replacing all
% occurencies of loadedStruct.(nameOfStruct) with the variable name, BUT
% it should be easier just to save it with right-click in workspace view
% and use as intended.
% 
%    input:
%%%%     -- .mat file with struct to load
%%%%     -- parameters
%%%%     -- criteria for rows to filter out in line 66
%    output:
%%%%     -- filtered .mat file with rows hitting criteria missing

%% settings %%

% the .mat file with a struct to load
pathOfMatFile = "C:\Users\[username]\...\file.mat";
% whether to save filtered .mat file
shallFileBeSaved = true;
% the path where to save the filtered .mat file
pathToSave = "C:\Users\[username]\...\file_filtered.mat";
% name of the struct that was saved/should be filtered
nameOfStruct = "xxx"; % e.g. "MultiHitStat"
% field name / column to apply critera on
fieldToCheck = "xxx"; %e.g. "FracMean"

%%%% CHANGE CRITERIA BELOW %%%%
%%%% you just have to edit the criteria below this line %%%%


%% load struct from .mat file into variable loadedStruct
try
    loadedStruct = load(pathOfMatFile);
% on errors
catch ME
    % show error
    disp(ME);
    % handle wrong file format
    if ME.identifier == "MATLAB:load:numColumnsNotSame"
        error("Error loading file. Is it a valid .mat file?");
    % handle file missing
    elseif ME.identifier == "MATLAB:load:couldNotReadFile"
        error('Error loading file. The path correct? It should end with ".mat"');
    end
end
disp("The file " + pathOfMatFile + " was loaded.")

%% initialize list of rows to delete
rowsToDelete = [];

%% go through the whole struct (2nd dimenstion of struct size)
% and check criteria
for i=1:size(loadedStruct.(nameOfStruct),2);
    % if this criteria is hit, then mark for deletion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% change this criterion as you need it %%%%%%%%%%%%%%%
    %%%%%%%%      here it is: "equals zero"      %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if loadedStruct.(nameOfStruct)(i).(fieldToCheck) == 0
        % mark for deletion
        rowsToDelete = [rowsToDelete, uint32(i)];
    end
end

%% delete (fill with empty) the marked rows of the struct
loadedStruct.(nameOfStruct)(rowsToDelete) = [];
%% print how many rows were deleted
disp("From " + string(structLengths(2)) + " rows, " ...
    + string(size(rowsToDelete,2)) + " were deleted.");

%% save new .mat if enabled
if shallFileBeSaved
    try
        save(pathToSave, "-struct", "loadedStruct");
    % on errors
    catch ME
        disp(ME);
        % handle parent directory missing
        if ME.identifier == "MATLAB:save:noParentDir"
            error("Error, the given directory to save the file " + ...
                "seems to be missing or incorrect.");
        end
    end
    % display if successful
    disp("Filtered .mat was saved to " + pathToSave);
else
    disp("The filtered struct was not saved. You may do it on your " + ...
        "own using the workspace.");
end