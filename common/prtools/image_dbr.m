function varargout = image_dbr(varargin)
%IMAGE_DBR M-file for image_dbr.fig
%      IMAGE_DBR, by itself, creates a new IMAGE_DBR or raises the existing
%      singleton*.
%
%      H = IMAGE_DBR returns the handle to a new IMAGE_DBR or the handle to
%      the existing singleton*.
%
%      IMAGE_DBR('Property','Value',...) creates a new IMAGE_DBR using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to image_dbr_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      IMAGE_DBR('CALLBACK') and IMAGE_DBR('CALLBACK',hObject,...) call the
%      local function named CALLBACK in IMAGE_DBR.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help image_dbr

% Last Modified by GUIDE v2.5 16-Aug-2008 23:50:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @image_dbr_OpeningFcn, ...
                   'gui_OutputFcn',  @image_dbr_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



function image_dbr_OpeningFcn(hObject, eventdata, handles, varargin)
% Main routine, this starts after creation of the guide

                              % Handling of input parameters
	handles.output = hObject;
	h = my_handles(handles);
	if length(varargin) < 4
		cclasf = meanc;           % default combiner
	else 
		cclasf = varargin{4};
	end
	                              
	if length(varargin) < 3
		clasf = knnc([],1);       % default classifier
	else
		clasf = varargin{3};
	end
	
	featsets = varargin{2};      % arrange feature sets in cells
	if ~iscell(featsets), featsets = {featsets}; end
	if ~iscell(clasf)
		clasf = repmat({clasf},1,length(featsets));
	else
		if length(clasf) ~= length(featsets)
			error('Number of classifiers and number of feature sets do not match')
		end
	end
	
	database = varargin{1};       % dataset or datafile, just used for image display

	m = size(database,1);
	labels = zeros(1,m);          % user give labels: 0 (no), 1 (target, 2 (outlier)
	S = randperm(m);              % initial set of images to be shown

	for j=1:4                     % feature weight slider bars default invisible
		set(h.feat_weight(j),'visible','off');
		set(handles.(['text' num2str(j+2)]),'visible','off');
	end
	                              % optimize location sliders in case of 2 or 3 sliders 
	if length(featsets) == 1, feath = 2; end
	if length(featsets) == 2, feath = [2 3]; end
	if length(featsets) == 3, feath = [1 2 3]; end
	if length(featsets) == 4, feath = [1 2 3 4]; end
	                        
	for j=1:length(featsets)      % initialize feature sets and related sliders
		if size(featsets{j},1) ~= m
			error('Feature set(s) have wrong size')
		end
		set(handles.(['text' num2str(feath(j)+2)]),'visible','on');
		set(handles.(['text' num2str(feath(j)+2)]),'string',getname(featsets{j}));
		% re-arrange datasets with features, use 'target' and 'outlier' as labels
		featsets{j} = prdataset(featsets{j});
		featsets{j} = featsets{j}*scalem(featsets{j},'variance');
		featsets{j} = prdataset(+featsets{j});
		featsets{j} = setlablist(featsets{j},char('target','outlier'));
		featsets{j} = setnlab(featsets{j},labels);
		set(h.feat_weight(feath(j)),'value',1);       % initial weights: 1
		set(h.feat_weight(feath(j)),'visible','on');
	end
	
	if length(featsets) == 1 % no feature slider bar in case of one feature set
		set(h.feat_weight(feath(1)),'visible','off');
		set(handles.(['text' num2str(feath(j)+2)]),'visible','off');
	end
	
	% store useful data in figure user data field
	guidata(hObject,{h,database,featsets,clasf,cclasf}); 
	
	set(h.classify,'userdata','reset'); % initial procedure: reset
	set(h.all,'value',0);               % initial ranking:
	set(h.unlabeled,'value',2);         % unlabeled objects
	labels_im = zeros(1,10);            % labels of shown images
	weights   = zeros(1,4);             % weights of feature sets
	Targets = [];
	Outliers= [];
	show_targ_out(database,Targets,1,h,'Targets');  % initialise target /
	show_targ_out(database,Outliers,1,h,'Outliers');% outlier images
	S = randperm(m);
	waitstate = prwaitbar;
	prwaitbar off
	                            % main loop
	while(1)                               % S, images to be shown
		show_images(database,S,1,h);         % show images on axes	
							                           % wait for the user 
		                                     % ..........
		uiwait                               % uiresume is activated by one of the 
		                                     % CLASSIFY / LABEL / RESET / QUIT buttonS.
		                                     % get labels and weights
		Targets = get(h.targets,'userdata'); % get Tagerts
		Outliers = get(h.outliers,'userdata');%abd Outliers
		if isempty(guidata(hObject))         % we need to stop (QUIT button is pressed)
			guidata(hObject,{S,Targets,Outliers});% 
			prwaitbar(waitstate);
			return                             % return ranking S, Targets and Outliers.
		end
		
		proch = proc(h);                   % classify / label / reset
		
		if ~strcmp(proch,'reset');
		
			for j=1:10                         % labels are 1-2, convert from 0-1 check boxes
				labels_im(j) = 2 - get(h.handlab(j),'value');
			end

			for j=1:length(featsets)           % get feature weights from sliders
				weights(j) = get(h.feat_weight(feath(j)),'value');
			end

			T = [1:m];                         % potential trainingset
	%		Targets = get(h.targets,'userdata');
	%		Outliers = get(h.outliers,'userdata');
			t = get(h.next,'userdata'); n = t{2};
       
			SS = [];                           % training set indices in case of 'label' procedure
			for j=1:10                         % update labeling
				SAct = S(j-1+n);                 % object number j in S
				if labels_im(j) == 1
					if isempty(find(Targets==SAct)), Targets = [Targets SAct]; end
				else
					if isempty(find(Outliers==SAct)), Outliers = [Outliers SAct]; end
				end
				SS = [SS SAct];
			end

			set(h.targets,'userdata',Targets);             % store indices of Targets
			set(h.outliers,'userdata',Outliers);           % and outliers
			labels = zeros(1,m);                           % construct labels training set
			labels(Targets)  = ones(1,length(Targets));
			labels(Outliers) = 2*ones(1,length(Outliers));

			if strcmp(proch,'label')                   % Label: learn from present set only
				if all(labels_im == labels_im(1))        % present set is uniformly labeled
					if labels_im(1) == 1 % no outliers, use remaining unlabeled dataset as outlier
						labels_train = ones(size(labels))*2; % Label entire set as outlier
						labels_train(SS) = ones(size(SS));   % Label present set as target
					else                 % no targets, use remaining dataset as target
						labels_train = ones(size(labels));   % Label entire set as target
						labels_train(SS) = 2*ones(size(SS)); % Label present set as outlier
					end
				else
					T = SS;
					labels_train = labels_im;
				end
				[S,W] = train_classify(featsets,T,labels_train,clasf,cclasf,weights);
			else                                       % Classify procedure: 
				T = find(labels > 0);                    % learn from all Targets and Outliers
				labels_train = labels(T);
				[S,W] = train_classify(featsets,T,labels_train,clasf,cclasf,weights);
			end
			U = find(labels(S) == 0);
			set(h.all,'userdata',{S,U});
			if (get(h.all,'value') == 0)
				S = S(U);
			end
		else                 % reset
			S = randperm(m);
		end
		% S is the ranked set of indices, either for all objects, or for just
		% the unlabeled objects
	end
return
		

function h = my_handles(handles)
% construct meaningfull handle names / numbers
	T = [1 2 3 4 7 8 9 10 11 12]; % correct checkbox numbering
	h.image   = zeros(1,10);
	h.handlab = zeros(1,10);
	h.feat_weight = zeros(1,4);
	for j=1:10
		h.image(j) = handles.(['axes' num2str(j)]);
		h.imnum(j) = handles.(['text' num2str(j+8)]);
		h.obnum(j) = handles.(['text' num2str(j+18)]);
		h.handlab(j) = handles.(['checkbox' num2str(T(j))]);
	end
	for j=1:4
		h.feat_weight(j) = handles.(['slider' num2str(j)]);
	end
	h.all  = handles.radiobutton4;
	h.unlabeled  = handles.radiobutton5;
	h.classify  = handles.pushbutton1;
	h.targets = handles.axes11;
	h.outliers = handles.axes12;
	h.target_slider = handles.slider5;
	h.outliert_slider = handles.slider6;
	h.target_title = handles.text1;
	h.outlier_title = handles.text2;
	h.target_obnum = handles.text29;
	h.outlier_obnum = handles.text30;
	h.target_delete = handles.pushbutton2;
	h.target_move = handles.pushbutton3;
	h.outlier_delete = handles.pushbutton5;
	h.outlier_move = handles.pushbutton4;
	h.next = handles.pushbutton8;
	h.previous = handles.pushbutton9;
return

function varargout = image_dbr_OutputFcn(hObject, eventdata, handles)
% takes care of final return of main return
	varargout = guidata(hObject);
	delete(hObject);
return

function pushbutton1_Callback(hObject, eventdata, handles)
% CLASSIFY button, this reactivates the UIWAIT function in the main loop
	s = guidata(hObject);
	h = s{1};
	set(h.classify,'userdata','classify');
	uiresume
return

function show_images(database,S,n,h)
% show images on axes
	Targets = get(h.targets,'userdata');
	Outliers = get(h.outliers,'userdata');
	im = data2im(database,S(n:n+9));
	               % the new images
	for j=1:10
		show_im(im,j,h.image(j));
		set(h.handlab(j),'value',1); % set target default
		set(h.imnum(j),'string',num2str(n+j-1));
		set(h.obnum(j),'string',num2str(S(n+j-1)));
	end
	set(h.next,'userdata',{S,n});
	              % set visibility of next / previous buttons
	if n > 1
		set(h.previous,'visible','on');
	else
		set(h.previous,'visible','off');
	end
	if n < (length(S)-9)
		set(h.next,'visible','on');
	else
		set(h.next,'visible','off');
	end
	               % the targets, show most recent one
	if ~isempty(Targets)
		show_targ_out(database,Targets,length(Targets),h,'Targets'); 
		set(h.target_obnum,'string',num2str(Targets(end)));
	end
	               % the outliers, show most recent one
	if ~isempty(Outliers)
		show_targ_out(database,Outliers,length(Outliers),h,'Outliers');
		set(h.outlier_obnum,'string',num2str(Outliers(end)));
	end
return

function show_targ_out(database,Pointers,n,h,name)
% show targets or outliers and update properties

	if strcmp(name,'Targets') % create single handle variable for targets and outliers
		h_axes = h.targets;
		h_slider = h.target_slider;
		h_title = h.target_title;
		h_delete = h.target_delete;
		h_move = h.target_move;
		h_obnum = h.target_obnum;
	else % Outliers
		h_axes = h.outliers;
		h_slider = h.outliert_slider;
		h_title = h.outlier_title;
		h_delete = h.outlier_delete;
		h_move = h.outlier_move;
		h_obnum = h.outlier_obnum;
	end
	
	if isempty(Pointers)  % all Targets / Outliers deleted, make all related objects invisible
		set(h_title,'String',['No ' name ' defined'])
		set(get(h_axes,'children'),'visible','off')
		set(h_axes,'visible','off');
		set(h_slider,'visible','off');
		set(h_delete,'visible','off');
		set(h_move,'visible','off');
		set(h_obnum,'visible','off');
		set(h_axes,'userdata',[]);
	else
		show_im(data2im(database,Pointers(n)),1,h_axes); % show image
		set(h_axes,'userdata',Pointers); 
		num_pointers = length(Pointers); % number of images
		set(h_title,'String',[num2str(num_pointers) ' ' name]); % set title
		set(h_obnum,'String',num2str(Pointers(n)));
		set(h_slider,'Min',0.99999/num_pointers);               % min (about 1/n)
		set(h_slider,'Max',1);                                  % max
		if num_pointers > 1                                     % step size
			set(h_slider,'SliderStep',[1/(num_pointers-1) 1/(num_pointers-1)]);
			set(h_slider,'value',n/num_pointers);
			set(h_slider,'visible','on');
		else                              % special case: one target/outlier image
			set(h_slider,'SliderStep',[1 1]);
			set(h_slider,'value',n/num_pointers);
			set(h_slider,'visible','off');
		end
		                                  % make image and buttons visible
		set(get(h_axes,'children'),'visible','on')
		set(h_delete,'visible','on');
		set(h_move,'visible','on');
		set(h_obnum,'visible','on');
		
	end
return

function proch = proc(h)
% get desired procedure
	proch = get(h.classify,'userdata'); % classify / label / reset
return

function show_im(im,n,h_axes)
% low level display routine
	if iscell(im), y = im{n};
	else, y = squeeze(im(:,:,:,n));
	end
	if size(y,3) == 1    % reset gray images to color
		y = cat(3,y,y,y);
	end
	axes(h_axes);        % activate the right axes
	image(y);            % display
	                     % prepare for callback
	set(get(gca,'children'),'ButtonDownFcn',{@resetlab,n})
	axis off             % make axes invisible
	axis equal
return

function resetlab(hObject, eventdata,n)
% Reset the image labels (target/outlier) by clicking in the image
	s = guidata(hObject);
	h = s{1};
	get(h.handlab(n),'value');
	set(h.handlab(n),'value',1 - get(h.handlab(n),'value'));
	get(h.handlab(n),'value');
return

function [S,W] = train_classify(featsets,T,labels_train,clasf,cclasf,weights);
% train, combine, classify and rank
% S will be the ranked object indices of all or unlabeled objects
	d = [];
	W = [];
	for j=1:length(featsets)
		b = featsets{j};
		trainset = setnlab(b(T,:),labels_train);
		%trainset = setprior(trainset,getprior(trainset,0));
		trainset = setprior(trainset,0);
		if ~isvaldset(trainset,2)
			v = trainset*knnc([],1);
		else
			v = trainset*clasf{j}*classc;
		end
		d = [d featsets{j}*v*weights(j)]; 
		W = [W;v];
	end
	d = d*cclasf;
	d = +d(:,'target');
	[dd,S] = sort(-d);
%	W = v*affine(weights)*cclasf; % to be corrected
	W = [];
return
	

function slider5_Callback(hObject, eventdata, handles)
% target slider
	s = guidata(hObject);
	h = s{1};
	database = s{2};
	Targets = get(h.targets,'userdata');
	n = round(get(h.target_slider,'value')*length(Targets));
	show_im(data2im(database,Targets(n)),1,h.targets);
	set(h.targets,'userdata',Targets); % image() destroys userdata, restore it
	set(h.target_title,'String',[num2str(length(Targets)) ' Targets']); % needed??? 
	set(h.target_obnum,'String',num2str(Targets(n)));
return

function slider6_Callback(hObject, eventdata, handles)
% outlier slider
	s = guidata(hObject);
	h = s{1};
	database = s{2};
	Outliers = get(h.outliers,'userdata');
	n = round(get(h.outliert_slider,'value')*length(Outliers));
	show_im(data2im(database,Outliers(n)),1,h.outliers);
	set(h.outliers,'userdata',Outliers); % image() destroys userdata, restore it
	set(h.outlier_title,'String',[num2str(length(Outliers)) ' Outliers']); % needed???
	set(h.outlier_obnum,'String',num2str(Outliers(n)));
return

function pushbutton2_Callback(hObject, eventdata, handles)
% delete target 
	s = guidata(hObject);
	h = s{1};
	database = s{2};
	Targets = get(h.targets,'userdata');
	n = round(get(h.target_slider,'value')*length(Targets));
	Targets(n) = [];
	n = max(n-1,1);
	show_targ_out(database,Targets,n,h,'Targets');
return

function pushbutton3_Callback(hObject, eventdata, handles)
% move target to outlier
	s = guidata(hObject);
	h = s{1};
	database = s{2};
	Targets = get(h.targets,'userdata');
	Outliers = get(h.outliers,'userdata');
	n = round(get(h.target_slider,'value')*length(Targets));
	Outliers = [Outliers Targets(n)];
	Targets(n) = [];
	n = max(n-1,1);
	show_targ_out(database,Outliers,length(Outliers),h,'Outliers');
	show_targ_out(database,Targets,n,h,'Targets');
return

function pushbutton4_Callback(hObject, eventdata, handles)
% move outlier to target
	s = guidata(hObject);
	h = s{1};
	database = s{2};
	Outliers = get(h.outliers,'userdata');
	Targets = get(h.targets,'userdata');
	n = round(get(h.outliert_slider,'value')*length(Outliers));
	Targets = [Targets Outliers(n)];
	Outliers(n) = [];
	n = max(n-1,1);
	show_targ_out(database,Targets,length(Targets),h,'Targets');
	show_targ_out(database,Outliers,n,h,'Outliers');
return

function pushbutton5_Callback(hObject, eventdata, handles)
% delete outlier 
	s = guidata(hObject);
	h = s{1};
	database = s{2};
	Outliers = get(h.outliers,'userdata');
	n = round(get(h.outliert_slider,'value')*length(Outliers));
	Outliers(n) = [];
	n = max(n-1,1);
	show_targ_out(database,Outliers,n,h,'Outliers');
return

function pushbutton6_Callback(hObject, eventdata, handles)
% Reset
	s = guidata(hObject);
	h = s{1};
	database = s{2};
	set(h.targets,'userdata',[]);
	set(h.outliers,'userdata',[]);
	show_targ_out(database,[],0,h,'Targets');
	show_targ_out(database,[],0,h,'Outliers');
	set(h.classify,'userdata','reset');
	uiresume
return



function pushbutton7_Callback(hObject, eventdata, handles)
% Quit
	guidata(hObject,[]);
	uiresume
return

function radiobutton4_Callback(hObject, eventdata, handles)
% All
	s = guidata(hObject);
	h = s{1};
	database = s{2};
	if (get(h.all,'value') == 0)
		set(h.unlabeled,'value',2);
	else
		set(h.unlabeled,'value',0);
	end
	t = get(h.all,'userdata');
	S = t{1};
	show_images(database,S,1,h)
return

function radiobutton5_Callback(hObject, eventdata, handles)
% unlabeled
	s = guidata(hObject);
	h = s{1};
	database = s{2};
	if (get(h.unlabeled,'value') == 0)
		set(h.all,'value',2);
	else
		set(h.all,'value',0);
	end
	t = get(h.all,'userdata');
	S = t{1};
	U = t{2};
	S = S(U);
	show_images(database,S,1,h)
return

function pushbutton8_Callback(hObject, eventdata, handles)
% Next
	s = guidata(hObject);
	h = s{1};
	database = s{2};
	t = get(h.next,'userdata');
	S = t{1};
	n = t{2};
	n = min(length(S)-9,n+10);
	show_images(database,S,n,h);
return
	

function pushbutton9_Callback(hObject, eventdata, handles)
% Previous
	s = guidata(hObject);
	h = s{1};
	database = s{2};
	t = get(h.next,'userdata');
	S = t{1};
	n = t{2};
	n = max(1,n-10);
	show_images(database,S,n,h);
return

function pushbutton10_Callback(hObject, eventdata, handles)
% Label
	s = guidata(hObject);
	h = s{1};
	set(h.classify,'userdata','label');
	uiresume
return

function pushbutton11_Callback(hObject, eventdata, handles)
% All Target
	s = guidata(hObject);
	h = s{1};
	for j=1:10
		set(h.handlab(j),'value',1);
	end
return

function pushbutton12_Callback(hObject, eventdata, handles)
% All Outlier
	s = guidata(hObject);
	h = s{1};
	for j=1:10
		set(h.handlab(j),'value',0);
	end
return

