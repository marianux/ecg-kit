function a2hbc_examples()
%%%%%%%%%%%%%%%%%%%%%%%%
% Examples using a2hbc %
%%%%%%%%%%%%%%%%%%%%%%%%
% 
a2hbc_root_path = [fileparts(mfilename('fullpath')) filesep];

% Here you have several examples using a2hbc from the command line.

%% MIT format

% To warm up you can start using the automatic mode.

a2hbc( ... 
    'recording_name', [ a2hbc_root_path '..' filesep 'example recordings' filesep '208.dat'], ...
    'recording_format', 'MIT', ... 
    'op_mode', 'auto');

keyboard

% then you can try the semi-assisted mode, the algorithm only ask you for
% help when it doubts about a cluster label.

a2hbc( ... 
    'recording_name', [ a2hbc_root_path '..' filesep 'example recordings' filesep '208.dat'], ...
    'recording_format', 'MIT', ... 
    'op_mode', 'slightly-assisted');   % <---- 

keyboard

% or the assisted mode. In this case, the algorithms always ask you to
% label each cluster found.

a2hbc( ... 
    'recording_name', [ a2hbc_root_path '..' filesep 'example recordings' filesep '208.dat'], ...
    'recording_format', 'MIT', ... 
    'op_mode', 'assisted');           % <---- 

keyboard

% if you want you can control some parameters of the algorithm, as the
% number of clusters to look for.

a2hbc( ... 
    'recording_name', [ a2hbc_root_path '..' filesep 'example recordings' filesep '208.dat'], ...
    'recording_format', 'MIT', ... 
    'NumOfClusters', 3, ...          % <---- 
    'op_mode', 'assisted');

keyboard

% After the processing if you want to continue tunning up some parameters
% you can try this in order to use the GUI.

a2hbc( ... 
    'recording_name', [ a2hbc_root_path '..' filesep 'example recordings' filesep '208.dat'], ...
    'recording_format', 'MIT', ... 
    'NumOfClusters', 3, ... 
    'InteractiveMode', true, ...      % <---- 
    'op_mode', 'assisted');

keyboard

%% Other formats, ISHNE


%% Other formats, AHA

