%% Load and Process Full Sleep EEG Recordings
% This script loads a full-night sleep EEG recording from a PhysioNet dataset.
% The EEG file is downloaded from https://physionet.org/content/capslpdb/1.0.0/.
% It extracts a specific channel and saves it as a MATLAB array.

%% Step 1: Specify the EDF File Path
edfFile = 'ins5.edf'; % Change filename if necessary

%% Step 2: Extract File Information
% Use edfinfo to retrieve metadata from the EDF file.
edfInfo = edfinfo(edfFile);

%% Step 3: Read EEG Data from the EDF File
% Use edfread to load the EEG recordings into a timetable
data = edfread(edfFile);

%% Step 4: Extract a Specific EEG Channel
% Define the column index corresponding to the desired EEG channel.
% Column index 4 corresponds to the C4-P4 channel.
columnIndex = 4; % Change this value to extract a different channel

% Select the specified column from the timetable (scalp channel)
selectedColumn = data{:, columnIndex};

% Concatenate all the 251x1 arrays in the selected column into a single array
new_X = vertcat(selectedColumn{:});

%% Step 5: Save the Extracted EEG Data
save('C4_P4.mat', 'new_X');

disp('EEG channel extraction complete. Data saved in C4_P4.mat.');
