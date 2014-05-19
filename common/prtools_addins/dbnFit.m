function model= dbnFit(X, numhid, y, varargin)
%fit a DBN to bianry data in X

%INPUTS: 
%X              ... data. should be binary, or in [0,1] interpreted as
%               ... probabilities
%numhid         ... list of numbers of hidden units
%y              ... List of discrete labels

%OUTPUTS:
%model          ... A cell array containing models from all RBM's

%varargin may contain options for the RBM's of this DBN, in row one by one
%for example:
%dbnFit(X, [500,400], opt1, opt2) uses opt1 for 500 and opt2 for 400
%dbnFit(X, [500,400], opt1) uses opt1 only for 500, and defaults for 400

numopts=length(varargin);
H=length(numhid);
model=cell(H,1);
cant_clases = length(unique(y));

if H>=2
    
    %train the first RBM on data
    if(numopts>=1)
        model{1}= rbmBB(X, numhid(1),varargin{1});
    else
        model{1}= rbmBB(X, numhid(1));
    end
    
%     %extra-debug
%     visualize(model{1}.W);
%     drawnow;
    
    %train all other RBM's on top of each other
    for ii=2:H
        if(numopts>=ii)
            model{ii}=rbmBB(model{ii-1}.top, numhid(ii), varargin{ii});
        else
            model{ii}=rbmBB(model{ii-1}.top, numhid(ii));
        end
        
%         %extra-debug
%         visualize(model{ii}.W);
%         drawnow;
        
    end

%     %the last RBM has access to labels too
%     if(numopts>=H)
%         model{H}= rbmFit(model{H-1}.top, numhid(end), y, varargin{H});
%     else
%         model{H}= rbmFit(model{H-1}.top, numhid(end), y);
%     end
    model{H}.Wc = 0.1*randn(cant_clases, size(model{H}.W,2) );
    model{H}.cc = 0.1*randn(1, cant_clases);
    model{H}.labels = 1:cant_clases;

    
    %fine-tuning of the DBN
    
    targets = zeros(length(y), cant_clases);
    for ii = 1:cant_clases
        targets(y == ii, ii) = 1;
    end
    
    for epoch = 1:200
    
        %%%%%%%%%%%%%%% PERFORM CONJUGATE GRADIENT WITH 3 LINESEARCHES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        max_iter=3;

        if epoch<6  % First update top-level weights holding other weights fixed. 
%             N = size(data,1);
%             XX = [data ones(N,1)];
%             w1probs = 1./(1 + exp(-XX*w1)); w1probs = [w1probs  ones(N,1)];
%             w2probs = 1./(1 + exp(-w1probs*w2)); w2probs = [w2probs ones(N,1)];
%             w3probs = 1./(1 + exp(-w2probs*w3)); %w3probs = [w3probs  ones(N,1)];

            testdata = X;
            for ii=1:length(model)
                testdata= rbmVtoH(model{ii}, testdata);
            end
            
            VV = colvec([ model{H}.Wc'; model{H}.cc ]);
            Dim = [ numhid(end); cant_clases];
            XX = minimize(VV,'CG_CLASSIFY_INIT',max_iter,Dim,testdata,targets);
            aux = reshape(XX, numhid(end)+1, cant_clases);
            model{H}.Wc = aux(1:end-1,:)';
            model{H}.cc = aux(end,:);

        else
            
%             VV = [w1(:)' w2(:)' w3(:)' w_class(:)']';
%             Dim = [l1; l2; l3; l4; l5];
            VV = [];
            for ii=1:H
                VV = [VV; colvec([model{ii}.W; model{ii}.b ]) ];
            end
            VV = [VV; colvec([model{ii}.Wc'; model{ii}.cc ])];
            
            Dim = [ size(X,2); colvec(numhid); cant_clases];

            XX = minimize(VV,'CG_CLASSIFY',max_iter,Dim,X,targets);

            xxx = 0;
            for ii=1:length(model)
                aux = reshape(XX(xxx+1:xxx+(Dim(ii)+1)*Dim(ii+1)),Dim(ii)+1,Dim(ii+1));
                model{ii}.W = aux(1:end-1,:);
                model{ii}.b = aux(end,:);
                xxx = xxx+(Dim(ii)+1)*Dim(ii+1);
            end
            aux = reshape(XX(xxx+1:xxx+(Dim(ii+1)+1)*Dim(ii+2)),Dim(ii+1)+1,Dim(ii+2));
            model{ii}.Wc = aux(1:end-1,:)';
            model{ii}.cc = aux(end,:);
            
%             w1 = reshape(X(1:(l1+1)*l2),l1+1,l2);
%             xxx = (l1+1)*l2;
%             w2 = reshape(X(xxx+1:xxx+(l2+1)*l3),l2+1,l3);
%             xxx = xxx+(l2+1)*l3;
%             w3 = reshape(X(xxx+1:xxx+(l3+1)*l4),l3+1,l4);
%             xxx = xxx+(l3+1)*l4;
%             w_class = reshape(X(xxx+1:xxx+(l4+1)*l5),l4+1,l5);

        end
        %%%%%%%%%%%%%%% END OF CONJUGATE GRADIENT WITH 3 LINESEARCHES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        prediction = dbnPredict(model, X);
        errors = sum(prediction~=y);
        fprintf(1, 'Errores %d\n', errors);
        
    end
    
    
%     %extra-debug
%     visualize(model{H}.W);
%     drawnow;
    
else
    
    %numhid is only a single layer... but we should work anyway
    if (numopts>=1)
        model{1}= rbmFit(X, numhid(1), y, varargin{1});
    else
        model{1}= rbmFit(X, numhid(1), y);
    end
end    

