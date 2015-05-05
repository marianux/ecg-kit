%% A progress bar class example
%
% Example:

%% start of algorithm
clear

pb = progress_bar('Progress bar demo', 'Start of algorithm');

pause(2)

%% initialization code

pb = pb.checkpoint('Initialization');

pause(4)

%% some iteration known a priori

pb.Loops2Do = 10;
pb.Title = 'Iterations known a priori';

for ii = 1:10
    
    pb = pb.start_loop();

    pause(3+randn(1))
    
    pb = pb.checkpoint('Step 1');

    pause(3+randn(1))
    
    pb = pb.checkpoint('Step 2');

    pause(3+randn(1))
    
    pb = pb.checkpoint('Step 3');

    pause(3+randn(1))
    
    pb = pb.end_loop();
    
end

%% some iteration unknown a priori

pb = pb.reset();
pb.Title = 'Iterations Unknown a priori';

for ii = 1:round(8+2*rand(1))
    
    pb = pb.start_loop();

    pause(3+randn(1))
    
    pb = pb.checkpoint('Step 1');

    pause(3+randn(1))
    
    pb = pb.checkpoint('Step 2');

    pause(3+randn(1))
    
    pb = pb.checkpoint('Step 3');

    pause(3+randn(1))
    
    pb = pb.end_loop();
    
end

%this clear and close all.
clear pb
