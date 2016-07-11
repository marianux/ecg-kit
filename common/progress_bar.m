%% (Internal) A progress bar class for showing evolution of a process to users
%
% Description:
% A progress-bar object to show the progress of a process. If the process
% is a loop, and the iteration number is known a priori, the object can inform
% the remaining time to finish the loop. If the iteration number is
% unknown, the object informs the mean duration of a loop.
% If the process is linear, the object can indicate the advance in the
% execution. The interface is very simple, and this example show the way of
% using it.
%
% Example:
% 
% %% start of algorithm
% clear
% 
% pb = progress_bar('Progress bar demo', 'Start of algorithm');
% 
% pause(2)
% 
% %% initialization code
% 
% pb.checkpoint('Initialization');
% 
% pause(4)
% 
% %% some iteration known a priori
% 
% pb.Loops2Do = 10;
% pb.Title = 'Iterations known a priori';
% 
% for ii = 1:10
%     
%     pb.start_loop();
% 
%     pause(3+randn(1))
%     
%     pb.checkpoint('Step 1');
% 
%     pause(3+randn(1))
%     
%     pb.checkpoint('Step 2');
% 
%     pause(3+randn(1))
%     
%     pb.checkpoint('Step 3');
% 
%     pause(3+randn(1))
%     
%     pb.end_loop();
%     
% end
% 
% %% some iteration unknown a priori
% 
% pb.reset();
% pb.Title = 'Iterations Unknown a priori';
% 
% for ii = 1:round(8+2*rand(1))
%     
%     pb.start_loop();
% 
%     pause(3+randn(1))
%     
%     pb.checkpoint('Step 1');
% 
%     pause(3+randn(1))
%     
%     pb.checkpoint('Step 2');
% 
%     pause(3+randn(1))
%     
%     pb.checkpoint('Step 3');
% 
%     pause(3+randn(1))
%     
%     pb.end_loop();
%     
% end
% 
% %this clear and close all.
% clear pb
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate  : 23/8/2011
% Last update: 24/3/2014
% Copyright 2008-2015
% 
classdef progress_bar < handle

    properties(GetAccess = private, Constant)

        %constant for the waitbar
        bar_position = [0.05    0.3    0.9    0.25];
        default_evolution = 0.1;
        long_loop_in_sec = 10; % seconds
        update_msg_time = 10; % seconds
        loops2weight = 20; % last 20 loops
        
    end

    properties ( Access = private )
        bUIpresent = usejava('desktop');
        wb_handle = [];
        wb_axes_hdl = [];
        LoopTimes = [];
        counter = 0;
        obj_tic = [];
        Cleanup_hdl = [];
        bPBcreated = false;
        childs
        parent
        update_counter
        time_elapsed
        time_weighting
    end

    properties(SetAccess = private, GetAccess = public)
        LoopMeanTime = [];
    end
    
    properties
        Message = '';
        Title = '';
        Loops2Do = [];
        LoopsDone = 0;
    end
    
    methods 

        function obj = progress_bar(Title, Message, wb_handle)

            if( nargin > 0 && ischar(Title) )
                obj.Title = Title;
            end
                
            if( nargin > 1 && ischar(Message) )
                obj.Message = Message;
            end
            
            if( obj.bUIpresent )
                % log to a waitbar
                
                obj.update_counter = 0;
                obj.time_elapsed = 0;
                
                if( nargin < 3 || ~ishandle(wb_handle) )
                    obj.wb_handle = waitbar(0, obj.Message, 'name', obj.Title );
                else
                    obj.wb_handle = wb_handle;
                end
                set(obj.wb_handle, 'Tag', 'progress_bar_class');
                obj.wb_axes_hdl = findobj(obj.wb_handle,'Type','Axes');
                set(obj.wb_axes_hdl, 'units','normalized' );
                set(obj.wb_axes_hdl, 'position',  obj.bar_position);
                
            else
                % TODO:log to stdout
                
            end
            
            obj.time_weighting = linspace(1,0.1,obj.loops2weight);
            obj.time_weighting = obj.time_weighting / sum(obj.time_weighting);
            
            %to clean and delete the waitbar.
%             obj.Cleanup_hdl = onCleanup(@()DoPBHouseKeeping(, obj.wb_handle));
            
            obj.bPBcreated = true;
            
        end
        
        function hide(obj)
            set(obj.wb_handle, 'Visible', 'off');
        end
        
        function show(obj)
            set(obj.wb_handle, 'Visible', 'on');
        end
        
        function checkpoint(obj, Message )
        
            if(nargin < 2)
                Message = [];
            end
            
            if(  ischar(Message) )
                obj.Message = Message;
            end
            
            if( isempty(obj.LoopMeanTime) )
                %first loop, learning times. Just to show time evolution.
                obj.counter = obj.counter + obj.default_evolution;
                currTime = 0;
            else
                %estimate progress
                currTime = toc(obj.obj_tic);
                
                if( isempty(obj.parent) )
                    aux_mean_time = obj.LoopMeanTime;
                    aux_LoopsDone = obj.LoopsDone;
                    aux_Loop2do = obj.Loops2Do;
                else
                    aux_mean_time = obj.parent.LoopMeanTime + obj.LoopMeanTime;
                    aux_LoopsDone = obj.parent.LoopsDone * obj.Loops2Do + obj.LoopsDone;
                    aux_Loop2do = obj.parent.Loops2Do * obj.Loops2Do ;
                end

                if( ~isempty(aux_Loop2do) && aux_Loop2do > 0 )
                    obj.counter = aux_LoopsDone/aux_Loop2do;
                else
                    obj.counter = currTime/aux_mean_time;
                end                    
                
            end 

            % take care always a waitbar to draw.
            if( obj.bUIpresent && ~ishandle(obj.wb_handle) )
                obj.wb_handle = waitbar(0);
                set(obj.wb_handle, 'Tag', 'progress_bar');
                obj.wb_axes_hdl = findobj(obj.wb_handle,'Type','Axes');
                set(obj.wb_axes_hdl, 'units','normalized' );
                set(obj.wb_axes_hdl, 'position', obj.bar_position );
            end

            if( obj.bUIpresent )
                waitbar( obj.counter - fix(obj.counter), obj.wb_handle, obj.Message );
            end
            
            if( isempty(obj.parent) )
                aux_LoopMeanTime = obj.LoopMeanTime;
            else
                aux_LoopMeanTime = obj.parent.LoopMeanTime;
            end
            
            if( isempty(aux_LoopMeanTime) )
                if( obj.bUIpresent && ~isempty(obj.obj_tic) )
                    %Learning phase
                    set(obj.wb_handle, 'Name', [ obj.Title '. Learning loop time ...' ]);
                end
            else
               
                if( isempty(obj.parent) )
                    aux_Time2Finish = (obj.Loops2Do-obj.LoopsDone) * obj.LoopMeanTime;
                else
                    aux_Time2Finish = (obj.parent.Loops2Do - obj.parent.LoopsDone) * obj.parent.LoopMeanTime + (obj.Loops2Do-obj.LoopsDone) * obj.LoopMeanTime;
                end                
                
                if( obj.bUIpresent )
                
                    obj.update_counter = obj.update_counter + obj.LoopMeanTime;
                    
                    if( obj.update_counter > obj.update_msg_time )
                        obj.update_counter = 0;
                        if( isempty(obj.Loops2Do) )
                            set(obj.wb_handle, 'Name', [ adjust_string(obj.Title, 30) ' - [' Seconds2HMS( aux_LoopMeanTime ) ' s/loop]']);
                        else
                            set(obj.wb_handle, 'Name', [ adjust_string(obj.Title, 30) ' - Finishing in ' Seconds2HMS(aux_Time2Finish) ]);
                        end
                    end
                    
                end
            end
                
        end
        
        function start_loop(obj)
        
            %start of loop. Reset timers
            obj.obj_tic = tic;

            if( isempty(obj.parent) )
                aux_mean_time = obj.LoopMeanTime;
            else
                aux_mean_time = obj.parent.LoopMeanTime + obj.LoopMeanTime;
            end
            
            if( obj.bUIpresent && ~isempty(aux_mean_time) && aux_mean_time > obj.long_loop_in_sec )
                % long process: progress within a loop.
                waitbar( 0, obj.wb_handle, 'Start of loop.' );
%             else
                % short process: total progress                 
            end
            
            
        end
        
        function end_loop(obj)
            %end of loop. Calculate averages.
            obj.LoopTimes = [obj.LoopTimes; toc(obj.obj_tic)];
            obj.LoopsDone = obj.LoopsDone + 1;
            if( length(obj.loops2weight) >= obj.loops2weight )
                % last loop is the most weighted
                obj.LoopMeanTime = mean([ (rowvec(obj.LoopTimes((end-obj.loops2weight+1):end))*flipud(colvec(obj.time_weighting))) (obj.time_elapsed / obj.LoopsDone) ]);
            else
                obj.LoopMeanTime = mean(obj.LoopTimes);
            end
            
            obj.time_elapsed = obj.time_elapsed + obj.LoopTimes(end);
            
            if( obj.bUIpresent )

                if( obj.LoopMeanTime > obj.long_loop_in_sec )
                    % long process: progress within a loop.
                    waitbar( 1, obj.wb_handle, 'End of loop.' );
    %             else
                    % short process: total progress                 
                end
                
            end
            
        end
        
        function pb_inner = AddInerLoop(obj, strMessage)
            
            pb_inner = progress_bar(obj.Title, strMessage, obj.wb_handle);
            pb_inner.parent = obj;
            obj.childs = [ obj.childs; pb_inner ];
            
        end
        
        function reset(obj)

            obj.LoopTimes = [];
            obj.LoopMeanTime = [];
            obj.counter = 0;
            obj.LoopsDone = 0;
            obj.Loops2Do = [];
            obj.obj_tic = [];
            obj.Message = '';
            obj.update_counter = 0;
            obj.time_elapsed = 0;
            
            if( obj.bUIpresent )
                waitbar( 0, obj.wb_handle, obj.Message );
                set(obj.wb_handle, 'Name', obj.Title);
            end
            
        end
        
        function set.Message(obj,value)
            if( ischar(value) )
                obj.Message = value;
                if( obj.bUIpresent && obj.bPBcreated )
                    waitbar( obj.counter - fix(obj.counter), obj.wb_handle, obj.Message );
                end
            else
                warning('progress_bar:BadArg', 'Message must be a string.');
            end
        end
        
        function set.Title(obj,value)
            if( ischar(value) )
                obj.Title = value;
                if( obj.bUIpresent && obj.bPBcreated )
                    set(obj.wb_handle, 'Name', obj.Title);
                end
            else
                warning('progress_bar:BadArg', 'Title must be a string.');
            end
        end
        
        function set.Loops2Do(obj,value)
            if( isempty(value) || (isnumeric(value) && value > 0 ) )
                obj.Loops2Do = value;
            else
                warning('progress_bar:BadArg', 'Loops2Do must be a number > 1.');
            end
        end
        
        function delete(obj)
            if( obj.bUIpresent )
                if( ishandle(obj.wb_handle) ) 
                    % waitbar close
                    delete(obj.wb_handle)
                end
            else

            end
            
        end
        
    end

end


