function success = setupSTEAMPaths()
% function success = setupSTEAMPaths()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Helper script for setting up PATH variables to access the scripts and classes
% defined in this package
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Archit Gupta
% Date: January 19, 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name_of_this_file = 'setupSTEAMPaths.m';
    STEAM_PATHS_SET_VALUE = 'PATHS_SET';
    persistent STEAM_PATHS_SET; 

    if (strcmp(STEAM_PATHS_SET, STEAM_PATHS_SET_VALUE))
        disp('It looks like STEAM paths have already been set.')
        success = 0;
        return;
    else
        current_dir = strcat(pwd, '/');
        
        % Check that the function has been called from the directory in which
        % the .m file is saved. All paths will be used relative to the location
        % of this file.

        if (~exist(strcat(current_dir, name_of_this_file)))
            error('PATH setup script not called from its file location. Please CD to the correct location and try again.')
        else
            path_locations_to_add = {'src/common', ...
                'src/utils', ...
                'src/device-models/MOS'};
            for loc = 1 : length(path_locations_to_add)
                addpath(strcat(current_dir, path_locations_to_add{loc}));
            end
        end
    end
    STEAM_PATHS_SET = STEAM_PATHS_SET_VALUE;
    disp('PATHS have been set up succesfully for STEAM.')
    success = 1;
end
