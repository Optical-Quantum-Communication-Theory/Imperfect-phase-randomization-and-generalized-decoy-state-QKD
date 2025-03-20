%% helper function for debuggin %%
% Converts from the labelling of (index1,index2) used in the code, where
% index2 stores number of clicks + 1, and index1 stores the specific click
% pattern (as ordered by nchoosek). This function outputs the exact click
% pattern.

function pattern = findClickPattern(index1, index2, numModes)
    clickPatternList = nchoosek([1:numModes],index2-1);

    pattern = clickPatternList(index1,:);
end

