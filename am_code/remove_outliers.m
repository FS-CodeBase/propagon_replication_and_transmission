% SCRIPT: REMOVE_OUTLIERS
% % Author: Fabian Santiago 
% % E-mail: FabianSantiago707@gmail.com
%
% DESCRIPTION
% %     Removes outliers from down sampled data, located in the folder
% %     propagon_data_down_sampled/, via IQR and saves the filtered data
% %     in the folder propagon_data_filtered_iqr/

% Locate propagon data
files_list = dir('propagon_data_down_sampled/');
files_list(1:2) = [];

% Save location
save_file = 'propagon_data_filtered_iqr/';

% Loop over each file
for File = 1:numel(files_list)
    % Load propagon data: propagon_data & sampling_times
    load([files_list(File).folder,'/',files_list(File).name])

    % Pre-allocate space to store data filtered for outliers.
    res = propagon_data;

    % Remove outliers
    for data = 1:numel(propagon_data)
        res{data} = rem_outliers_iqr(propagon_data{data});
    end
    propagon_data = res;
    % Save filtered data
    save([save_file,files_list(File).name(1:end-4),'_fil.mat'],...
            'propagon_data','sampling_times')
end