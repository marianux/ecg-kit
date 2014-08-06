%PRWAITBAR Report PRTools progress by single waitbar
%
%   H = PRWAITBAR(N,M,TEXT)
%   H = PRWAITBAR(N,TEXT,FLAG)
%   S = PRWAITBAR
%
% INPUT
%   N      Integer, total number of steps in loop
%   M      Integer, progress in number of steps in loop
%   TEXT   Text to be displayed in waitbar
%   FLAG   Flag (0/1)   
%
% OUTPUT
%   H      Waitbar handle
%   S      Status PRWAITBAR ('on' or 'off')
%
% DESCRIPTION
% This routine may be used to report progress in PRTools experiments.
% It detects and integrates levels of loops. The following calls are
% supported:
%
%   PRWAITBAR(N,TEXT)        initialize loop
%   PRWAITBAR(N,TEXT,FLAG)   initialize loop if FLAG == 1
%   PRWAITBAR(N,M)           update progress
%   PRWAITBAR(N,M,TEXT)      update info
%   PRWAITBAR(0)             closes loop level
%   PRWAITBAR OFF            removes waitbar
%   PRWAITBAR ON             switches waitbar on again
%   PRWAITBAR REPORT         no waitbar, report progress in command window
%   PRWAITBAR                reset prwaitbar
%
% A typical sequence of calls is:
%    .....
%    prwaitbar(nfolds,'cross validation')
%    for j = 1:nfolds
%        prwaitbar(nfolds,j,['cross validation, fold ' int2str(j)])
%        ....
%    end
%    prwaitbar(0)
%    .....
% Calls to PRTWAITBAR may be nested and all progress is merged into a single
% waitbar. In this PRWAITBAR differs from Matlab's WAITBAR.
% A typical example can be visualised by:
%
%    prcrossval({gendatb,gendath},{svc,loglc,fisherc},5,2)
%
% PRWAITBAR may increase the computation time with 25%. For long running
% processes PRWAITBAR REPORT (e.g. in a screen window) is recommended.
%
% Some other, more high-level routines calling PRWAITBAR are
% PRWAITBARINIT, PRWAITBARNEXT, PRWAITBARONCE, PREIG, PRINV, PRPINV, PRRANK,
% PRSVD, PRCOV
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PRWAITBARINIT, PRWAITBARNEXT, PRWAITBARONCE

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function varargout = prwaitbar(n,m,text)

persistent N                % array, N(DEPTH) is the loop size at level DEPTH
persistent M                % array, M(DEPTH) is the present counter values at level DEPTH
persistent DEPTH            % nesting DEPTH
persistent WHANDLE          % handle of the waitbar
persistent MESS             % cell array with message to be displayed at level DEPTH
persistent STAT             % flag, if STAT = 0, PRWAITBAR is switched off
persistent WPOS             % position of the waitbar on the screen
persistent REP              % flag, if REP = 1, PRWAITBAR reports in the command window
persistent FNAME            % cell array with names of the calling routine at level DEPTH
persistent DOIT             % flag that may switch off tracking progress at deeper levels
persistent NSKIP            % level at which tracking progress has been switched off
persistent OLDPROG
persistent OLDPROG1         % previous progress
persistent OLDPROG2         % previous progress
persistent OLDPROG3         % previous progress for report

% disp(['prwaitbar ' callername])
% disp(nargin)
% if nargin > 0, disp(n),   end
% if nargin > 1, disp(m),   end
% if nargin > 2, disp(text),end
% disp('-------------')
if isempty(DEPTH)           % initialisation of persistent variables
	DEPTH = 0;
	MESS = cell(1,10); 
	FNAME = cell(1,10);
	STAT = 1;
	DOIT = 1;
  OLDPROG = 0;
end

decrease_depth = 0;
if nargin == 0                            % ---- call: prwaitbar ----
	if nargout == 0
		%prwaitbar off
		DEPTH = 0;
		STAT = 1;
    REP = 0;
		DOIT = 1;
		MESS = cell(1,10); 
		OLDPROG1 = 0;
		OLDPROG2 = 0;
		OLDPROG3 = 0;
		err.message = '##';
		lasterror(err);
% 		if ishandle(WHANDLE) 
% 			delete(WHANDLE);
% 		end
	else                                    % ---- call: out = prwaitbar ----
		if nargout == 1
			if STAT == 0, varargout = {'off'};
			else varargout = {'on'};
			end
		else
			varargout = {WHANDLE,N,M,DEPTH,MESS,STAT,REP,WPOS,FNAME,DOIT,NSKIP};
		end
		return
	end
	
elseif nargin == 1
	
	if isstr(n) & strcmp(n,'off')           % ---- call: prwaitbar off ----
		
		if ~isempty(WHANDLE) & ishandle(WHANDLE)
			WPOS = get(WHANDLE,'position');     % save position as user likes it
			close(WHANDLE);
		end
		N = []; M = []; DEPTH = 0;            % set variables at initial values
		WHANDLE = [];
		MESS = cell(1,10); 
		STAT = 0;
 		REP = 0;
		return
		
	elseif isstr(n) & strcmp(n,'on')        % ---- call: prwaitbar on ----
		STAT = 1;                            % restart waitbar
		REP = 0;
    DEPTH = [];
		return
		
  elseif isstr(n) & strcmp(n,'report')    % ---- call: prwaitbar report ----
		prwaitbar off;
    err.message = '##';
    lasterror(err);
		STAT = 1;                            
		REP = 1;
		DEPTH = 0;
		if ishandle(WHANDLE)
			delete(WHANDLE);
		end
		WHANDLE = [];
    OLDPROG3 = 0;
		return
    
	elseif n == 0 & STAT                    % ---- call: prwaitbar(0) ----
		if ~DOIT 
			if (NSKIP == DEPTH)                 % we are back at the right level, so
				DOIT = 1;                         % restart tracking progress
			end
			DEPTH = max(DEPTH - 1,1);           % counting problem !!! caused by catch-try loops
      return                              % let us continue with DEPTH = 1
		end
		if DEPTH < 1                        % wrong call, return
			return
		end
		M(DEPTH) = N(DEPTH);                % loop almost ready
		decrease_depth = 1;                 % decrease depth after plotting progress
	else
		return
	end
	
elseif ~STAT                            % next commands can be skipped if switched off
	return
	
elseif nargin == 2
	
	if isstr(m)                           % ---- call: prwaitbar(n,text) ----
			                                  % note: this is a typical call to
			                                  % start a waitbar cycle
		
		DEPTH = DEPTH + 1;                  % we nest one level deeper
		if ~DOIT                            % skip if we dont track progress
			return
		end 

		err = lasterror;                    
		% we have to find out whether we continue a proper set of prwaitbar calls, or whether
		% we have to recover from some error or possible interrupt and have to re-initialise 
		% the admin
	
    if (isempty(err.message)) | (~ishandle(WHANDLE)) | ...
            (~strcmp(err.message,'##') &  ...
            isempty(strfind(err.identifier,'prtools:')))
        % yes, restart after interrupt
      if REP
        prwaitbar report;
      else
        prwaitbar;
      end
      DEPTH = 1;
%     elseif strcmp(err.message,'')
%       err.message = '##';
%       lasterror(err);
    end
			                                   % find the name of the calling routine
		[ss,ii] = dbstack;
		if length(ss) == 1
			name = '';                         % call from command line!!?? debugging?
		else
			[path,name] = fileparts(ss(2).name);
			if strcmp(name,mfilename) | strncmp(name,'prwait',6) % search deeper in case of internal call
				[path,name] = fileparts(ss(3).name);
			end
		end
		FNAME{DEPTH} = name;                 % store it, and display it in waitbar
		if ~REP
			set(WHANDLE,'name',['                         PRWAITBAR:   ' name]);
		end
		
		N(DEPTH) = n;                        % complete the admin
		M(DEPTH) = 0;
		MESS{DEPTH} = m;                     % note: m contains text
    if DEPTH == 1 
      OLDPROG = 0;
    end
		
	elseif DOIT                            % ---- call: prwaitbar(n,m) ----
		DEPTH = max(DEPTH,1);                % store admin
		N(DEPTH) = n;
		M(DEPTH) = m-1;                      % store m-1, assuming call is in the start of loop
	else                                   % skip and return if progress should not be reported
		return
	end
	
elseif nargin == 3
	
	if isstr(m)                            % ---- call: prwaitbat(n,text,flag) ----
		if ~DOIT                             % there is a flag, but we already skip progress               
			DEPTH = DEPTH + 1;                 % just keep track of loop nesting.
		else                                 % we are in reporting mode  
			DOIT = (text > 0);                 % check flag 
			if ~DOIT                           % stop progress reporting!!!
				DEPTH = DEPTH + 1;               % keep track of loop nesting
				NSKIP = DEPTH;                   % store level
			else                               % flag does not apply,
				prwaitbar(n,m);                  % just print waitbar
			end
		end
		return
		
	else                                   % ---- call: prwaitbar(n,m,text) ----
		
		if ~DOIT                             % skip progress report if needed
			return
		end
		
		if isempty(WHANDLE) & ~REP           % switch waitbar on if needed
			prwaitbar;                         % this only happens if the loop is not properly opened
    end
    if DEPTH > 0
      N(DEPTH) = n;                      % update admin
      M(DEPTH) = m-1;
      MESS{DEPTH} = text;
    end
	end
end

if ~STAT | ~DOIT                         % dont display, (just to be sure)
	return
end

% waitbar creation and display

if isempty(WHANDLE) & ~REP               % if no waitbar, make one
	
	s = sprintf(' \n \n \n' );             % room for 4 levels of text
	
	if ~isempty(WPOS)
		h = waitbar(0,'','position',WPOS);
		h = waitbar(0,h,s);
	else
		h = waitbar(0,s);                    % create waitbar
		WPOS = get(h,'position');
	end
	OLDPROG1 = 0;
	OLDPROG2 = 0;
	if DEPTH == 0, fnam = [];
	else, fnam = FNAME{DEPTH}; end         % retrieve name of calling routine
 	set(h,'name',['                         PRWAITBAR:   ' fnam])
	%set(h,'HandleVisibility','on');
	WHANDLE = h;
	
	err.message = '##';                    % store for proper continuation in next call
	lasterror(err);
	
end
	
progress = 0;                            % progress
if DEPTH == 0                            % still on level zero
	s = sprintf(' \n \n \n' );             % no message
else
	for j=DEPTH:-1:1                       % run over all levels for progress
		progress = (progress + M(j))/(N(j)); % compute total progress
	end
	s = wtext(MESS,N,DEPTH,REP);           % construct message
  progress = max(progress,OLDPROG);      % avoid progress to decrease
end

significant1 = (progress - OLDPROG1 >= 0.005 | OLDPROG1 == 0 | decrease_depth);
significant2 = (progress - OLDPROG2 >= 0.0005);
if ishandle(WHANDLE)
	if significant1                          % show significant updates only                
		waitbar(progress,WHANDLE,s,'position',WPOS);            % update waitbar
		OLDPROG1 = progress;
		OLDPROG2 = progress;
	elseif significant2                       % update text more often
		set(get(get(WHANDLE,'children'),'title'),'string',s); drawnow
		OLDPROG2 = progress;
	end
else                                      % handle was not proper (anymore)
	if REP 
    if progress - OLDPROG3 >= 0.001
      fprintf(['\n%6.3f ' s],progress)
      OLDPROG3 = progress;
    end
  else % reconstruct waitbar
    if isempty(WPOS)
      WHANDLE = waitbar(progress,s);
    else
      WHANDLE = waitbar(progress,s,'position',WPOS);    
    end
	end
	OLDPROG1 = progress;
	OLDPROG2 = progress;
end
OLDPROG = progress;

if decrease_depth                         % update admin
	DEPTH = DEPTH - 1;
	if DEPTH == 0                           % take care that new waitbar starts fresh
		if ~REP
			WPOS = get(WHANDLE,'position');       % but in old position
      %delete(WHANDLE);
      %we don't delete but preserve the waitbar for future use
      %so the following cleanup is needed
      p=get(get(WHANDLE,'children'),'children');
      set(p(2),'erasemode','normal');
      waitbar(0,WHANDLE,sprintf(' \n \n \n' ));
      set(WHANDLE,'name',['                         PRWAITBAR:      ']);
      set(p(2),'erasemode','none');
      prwaitbar;  
      lasterror('reset');
      %now we are ready to leave
    end
	else
		if significant1
	 		set(WHANDLE,'name',['                         PRWAITBAR:   ' FNAME{DEPTH}]);
		end
	end
end

if nargout > 0
	if REP
		varargout = 0;
	else
		varargout = {WHANDLE};
	end
end

return

function s = wtext(MESS,N,DEPTH,REP)

	t = cell(1,DEPTH);
	n = 0;
	for j=1:DEPTH
		if N(j) > 1
			n = n+1;
			t{n} = MESS{j};
		end
	end

	if REP
		if n > 0
			s = sprintf('%s\n',t{1:n});
		else
			s = '';
		end
	elseif n == 0
		s = sprintf([' \n \n \n \n ']);
	elseif n == 1
		s = sprintf([' \n \n' t{1} '\n ']);
	elseif n == 2
		s = sprintf([' \n' t{1} '\n' t{2} '\n ']);
	elseif n == 3
		s = sprintf([' \n' t{1} '\n' t{2} '\n' t{3}]);
	else
		s = sprintf([t{n-3} '\n' t{n-2} '\n' t{n-1} '\n' t{n}]);
	end
	
return
	
