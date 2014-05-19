
function w = dbnc(dsTrain, varargin)

	prtrace(mfilename);

	% Settings for the different training algorithms.
    if(~exist('nnet','dir'))
		error('Neural network toolbox not found')
    end

    cInput_types = {'gaussian', 'logsigm'};
    
    %argument definition
    p = inputParser;   % Create instance of inputParser class.
    p.addRequired('dsTrain', @(x)( isdataset(x)));
    
    p.addParamValue('input_type', 'logsigm', @(x)( any(strcmpi(x,cInput_types))) );
    p.addParamValue('num_hidden', 3, @(x)( isnumeric(x) && x > 0 ));
    p.addParamValue('size_hidden', [], @(x)( (all(isnumeric(x)) && all(x > 0)) || isempty(x) ) );
    p.addParamValue('Initialization_dataset', [], @(x)( isdataset(x) )  );
    p.addParamValue('Validation_dataset', [], @(x)( isdataset(x) )  );

	% Check arguments
    try
        p.parse( dsTrain, varargin{:} );
    catch MyError
        rethrow(MyError);    
    end

    dsTrain = p.Results.dsTrain;
    dsInitialization = p.Results.Initialization_dataset;
    dsValidation = p.Results.Validation_dataset;
    num_hidden = p.Results.num_hidden;
    size_hidden = p.Results.size_hidden;
    input_type = p.Results.input_type;

    % Dont know why this variable uses a lot of bytes to store at disk.
    clear p

    %Check the datasets
    islabtype(dsTrain,'crisp');
    isvaldfile(dsTrain,1,2); 							% At least 1 object per class, 2 classes
    dsTrain = testdatasize(dsTrain);
    %force domain to be in [0-1]
    w_scale_domain = scalem(dsTrain, 'domain');
    dsTrain = dsTrain * w_scale_domain;
    
    if( isempty(dsValidation) )
        dsValidation = dsTrain;
    else
        dsValidation = testdatasize(dsValidation);
        iscomdset(dsTrain, dsValidation);   							% Check whether training and tuning set match
        dsValidation = dsValidation * w_scale_domain;
    end
    
    if( isempty(dsInitialization) )
        dsInitialization = dsTrain;
    else
        dsInitialization = dsInitialization * w_scale_domain;
    end
    
    [train_m, train_k, train_c] = getsize(dsTrain); 
    [ train_labels, train_lablist] = getnlab(dsTrain); 
    
    if( isempty(size_hidden) )
        % num_hidden layers each of train_k size by default
        size_hidden = repmat(train_k, num_hidden, 1);
    end
    
    %train the first RBM on data
    top = +dsInitialization;
    %train all other RBM's on top of each other
    for ii=1:num_hidden
        [model{ii} , top ] = rbmBB(top, size_hidden(ii));
    end
    clear top

    model{num_hidden}.Wc = 0.1*randn(train_c, size(model{num_hidden}.W,2) );
    model{num_hidden}.cc = 0.1*randn(1, train_c);
    model{num_hidden}.labels = 1:train_c;
    
    %fine-tuning of the DBN
    targets = zeros(train_m, train_c);
    for ii = 1:train_c
        targets(train_labels == ii, ii) = 1;
    end
    
    for epoch = 1:20

        %%%%%%%%%%%%%%% PERFORM CONJUGATE GRADIENT WITH 3 LINESEARCHES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        max_iter=3;

        if epoch<6  % First update top-level weights holding other weights fixed. 

            top = +dsTrain;
            for ii=1:num_hidden
                top = rbmVtoH(model{ii}, top);
            end

            VV = colvec([ model{num_hidden}.Wc'; model{num_hidden}.cc ]);
            Dim = [ size_hidden(end); train_c];
            XX = minimize(VV,'CG_CLASSIFY_INIT',max_iter,Dim,top,targets);
            aux = reshape(XX, size_hidden(end)+1, train_c);
            model{num_hidden}.Wc = aux(1:end-1,:)';
            model{num_hidden}.cc = aux(end,:);

        else

            VV = [];
            for ii=1:num_hidden
                VV = [VV; colvec([model{ii}.W; model{ii}.b ]) ];
            end
            VV = [VV; colvec([model{ii}.Wc'; model{ii}.cc ])];

            Dim = [ train_k; colvec(size_hidden); train_c];

            XX = minimize(VV,'CG_CLASSIFY',max_iter,Dim, +dsTrain,targets);

            xxx = 0;
            for ii=1:num_hidden
                aux = reshape(XX(xxx+1:xxx+(Dim(ii)+1)*Dim(ii+1)),Dim(ii)+1,Dim(ii+1));
                model{ii}.W = aux(1:end-1,:);
                model{ii}.b = aux(end,:);
                xxx = xxx+(Dim(ii)+1)*Dim(ii+1);
            end
            aux = reshape(XX(xxx+1:xxx+(Dim(ii+1)+1)*Dim(ii+2)),Dim(ii+1)+1,Dim(ii+2));
            model{ii}.Wc = aux(1:end-1,:)';
            model{ii}.cc = aux(end,:);

        end    

    end

    % Create mapping.

    w = w_scale_domain*mapping('dbnPredict', 'trained', model, train_lablist, train_k, train_c);
    w = setname(w, 'DBN classifier');
%     w = setcost(w,a);

return
