
%%%%%%%%%%%%%%%%%%%%%%%%
% Examples using a2hbc %
%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Here you have several examples using a2hbc from the command line.

%% MIT format

% To warm up you can start using the automatic mode.

a2hbc( ... 
    'recording_name', [ '.' filesep 'example recordings' filesep '208.dat'], ...
    'recording_format', 'MIT', ... 
    'op_mode', 'auto');

% then you can try the semi-assisted mode, the algorithm only ask you for
% help when it doubts about a cluster label.

a2hbc( ... 
    'recording_name', [ '.' filesep 'example recordings' filesep '208.dat'], ...
    'recording_format', 'MIT', ... 
    'op_mode', 'slightly-assisted');   % <---- 

% or the assisted mode. In this case, the algorithms always ask you to
% label each cluster found.

a2hbc( ... 
    'recording_name', [ '.' filesep 'example recordings' filesep '208.dat'], ...
    'recording_format', 'MIT', ... 
    'op_mode', 'assisted');           % <---- 

% if you want you can control some parameters of the algorithm, as the
% number of clusters to look for.

a2hbc( ... 
    'recording_name', [ '.' filesep 'example recordings' filesep '208.dat'], ...
    'recording_format', 'MIT', ... 
    'NumOfClusters', 3, ...          % <---- 
    'op_mode', 'assisted');

% After the processing if you want to continue tunning up some parameters
% you can try this in order to use the GUI.

a2hbc( ... 
    'recording_name', [ '.' filesep 'example recordings' filesep '208.dat'], ...
    'recording_format', 'MIT', ... 
    'NumOfClusters', 3, ... 
    'InteractiveMode', true, ...      % <---- 
    'op_mode', 'assisted');


%% Other formats, ISHNE


%% Other formats, AHA

