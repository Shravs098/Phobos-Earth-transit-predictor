function mkdir_if_missing(folder)
% MKDIR_IF_MISSING  Create a folder only if it doesn't already exist.
if ~exist(folder,'dir')
    mkdir(folder);
end
end