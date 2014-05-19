%ROC Receiver-Operator Curve
% 
%   E = ROC(B,C,N)
%
% INPUT
%   A  Dataset
%   W  Trained classifier, or
%   B  Classification result, B = A*W*CLASSC
%   C  Index of desired class (default: C = 1)
%   N  Number of points on the Receiver-Operator Curve (default: 100)
%
% OUTPUT
%   E  Structure containing the error of the two classes
%
% DESCRIPTION
% Computes N points on the receiver-operator curve of the classifier W for
% class C in the labeled dataset B, which is typically the result of
% B = A*W; or for the dataset A labelled by applying the (cell array of)
% trained classifiers W.
%
% Note that a Receiver-Operator Curve is related to a specific class (class C)
% for which the errors are plotted horizontally. The total error on all other
% classes is plotted vertically. The class index C refers to its position in
% the label list of the dataset (A or B). It can be found by GETCLASSI.
%
% The curve is computed for N thresholds of the posteriori probabilities
% stored in B. The resulting error frequencies for the two classes are
% stored in the structure E. E.XVALUES contains the errors in the first
% class, E.ERROR contains the errors in the second class. In multi-class
% problems these are the mean values in a single class, respectively the
% mean values in all other classes. This may not be very useful, but not
% much more can be done as for multi-class cases the ROC is equivalent to a
% multi-dimensional surface.
%
% Use PLOTE(E) for plotting the result. In the plot the two types of error
% are annotated as 'Error I' (error of the first kind) and 'Error II' (error
% of the second kind). All error estimates are weighted according the class
% prior probabilities. Remove the priors in A or B (by setprior(A,[])) to
% produce a vanilla ROC.
%
% EXAMPLES
%	Train set A and test set T:
%	  B = T*NMC(A); E = ROC(T,50); PLOTE(E); % Plots a single curve
%	  E = ROC(T,A*{NMC,UDC,QDC});  PLOTE(E); % Plots 3 curves
%

function [roc, thr] = roc_detection(a,n)

	% Depending on the call, CLAS may the third or second argument.
	% and N the third or the fourth.
	
	if nargin < 2 || isempty(n), n = 100; end
	
	datname = getname(a);
	lablist = getlablist(a,'string');
    
    required_labs = {'FP' 'TP'};
    aux_lablist = intersect(cellstr(lablist), required_labs );

    cant_required_labs = length(required_labs);

    if( length(aux_lablist) ~= cant_required_labs )
        fprintf(2, ['Esta funcion esta pensada para usarse en datasets de deteccion con labels:\n' colvec([char(required_labs) repmat('\n', cant_required_labs, 1)]')' ] );
        error();
    end
    
    clas = find( strcmpi(cellstr(lablist), 'TP') );
	clasname = lablist(clas,:);
	%DXD: also check the class sizes:
	cs = classsizes(a);
	if any(cs == 0)
		error('Ambas clases deben contener ejemplos');
	end

	% Set up the ROC structure.


    % If a cell array of classifiers was given, apply each one here.

    a = a*normm; % make sure we have a normalised classification matrix

    [m,c] = size(a); 
    nlab = getnlab(a); 
    d = sort(a(:));

    % Attempt to compute a good threshold range on the first classified
    % dataset.
%     thr = [max(0, min(d)-eps) rowvec(d(round(linspace(2,length(d)-1, n-1)))) 1];
    thr = linspace(0,1,n+1);

    % NLAB_OUT will be one where B is larger than THR.
    I = matchlablist(getlablist(a),getfeatlab(a)); % Correct possible changes class orders
    nlab_out = (repmat(+a(:,I(clas)),1,n+1) >= repmat(thr,m,1));

    % aciertos will be 1 where the numeric label is unequal to NLAB_OUT
    % (i.e., where errors occur).
    bTP = nlab == clas;
    aciertos = (repmat(bTP,1,n+1) == nlab_out);

    sensitivity = mean(aciertos(bTP,:),1); % S
    pospred = sum(aciertos(bTP,:)) ./ sum(nlab_out); % +P

    roc = [ colvec(sensitivity) colvec(pospred) ];
  
    mod = sqrt(sum(roc.^2,2));
    [max_mod max_mod_idx] = max(mod);
    
    figure(100);
    h = subplot(1,2,1);
    plot(h(1), sensitivity(:), pospred(:), 'bo-' )
    hold(h(1), 'on');
    plot(h(1), sensitivity(max_mod_idx), pospred(max_mod_idx), 'rx:', 'MarkerSize',11)
    plot(h(1), sensitivity(max_mod_idx), pospred(max_mod_idx), 'mo:', 'MarkerSize',11)
    hold(h(1), 'off');
    axis(h(1), 'square');
    box(h(1), 'off')
    xlabel(h(1), 'Sensitivity')
    ylabel(h(1), 'Positive predictivity')
%     x_lim = xlim();
%     y_lim = ylim();
%     xlim([min([x_lim ylim]) 1.1]);
%     xlim([min([x_lim ylim]) 1.1]);
    
    h(2) = subplot(1,2,2);
    plot(h(2), thr, mod, 'bo-' )
    hold(h(2), 'on');
    plot(h(2), thr(max_mod_idx), mod(max_mod_idx), 'rx', 'MarkerSize',11)
    plot(h(2), thr(max_mod_idx), mod(max_mod_idx), 'mo', 'MarkerSize',11)
    plot(h(2), [0.5 0.5], ylim(), 'k--')
    hold(h(2), 'off');
    box(h(2), 'off')
    ylabel(h(2), 'mod(S,P+)')
    xlabel(h(2), 'Operating point')

    ConfusionMat = [ sum(aciertos(~bTP,max_mod_idx))  sum(~aciertos(~bTP,max_mod_idx)); sum(~aciertos(bTP,max_mod_idx)) sum(aciertos(bTP,max_mod_idx)) ];
    
    DisplayResults('dsResult', ConfusionMat, 'SupportDataset', a);
    