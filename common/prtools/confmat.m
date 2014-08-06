%CONFMAT Construct confusion matrix
% 
%  [C,NE,LABLIST1,LABLIST2] = CONFMAT(LAB1,LAB2,METHOD,FID)
%
% INPUT
%  LAB1        Set of labels
%  LAB2        Set of labels
%  METHOD      'count' (default) to count number of co-occurences in
%	             LAB1 and LAB2, 'disagreement' to count relative
%		           non-co-occurrence.
%  FID         Write text result to file
%
% OUTPUT
%  C           Confusion matrix
%  NE          Total number of errors (empty labels are neglected)
%  LABLIST1    Label list for LAB1
%  LABLIST2    Label list for LAB2
%
% DESCRIPTION
% Constructs a confusion matrix C between two sets of labels LAB1 
% (corresponding to the rows in C) and LAB2 (the columns in C). The order of 
% the rows and columns is returned in LABLIST. NE is the total number of 
% errors (sum of non-diagonal elements in C).
%
% If METHOD = 'count' (default), co-occurences in LAB1 and LAB2 are counted 
% and returned in C. 
% For METHOD = 'disagreement', the relative disagreement is returned in NE,
% and is split over all combinations of labels in C (such that the rows sum
% to 1). (The total disagreement for a class equals one minus the 
% sensitivity for that class as computed by TESTC).
% For METHOD = 'ids' a cell array C is returned in which C(i,j) contains
% the indices of the objects for which LAB1 equals LABLIST1(i,:) and LAB2
% equals LABLIST2(j,:).
%
%   [C,NE,LABLIST] = CONFMAT(D,METHOD)
%
% If D is a classification result D = A*W, the labels LAB1 and LAB2 are 
% internally retrieved by CONFMAT before computing the confusion matrix.
%
%   C = CONFMAT(D)
%
% This call also applies in case in D = A*W the dataset A has soft labels 
% W is trained by a soft labeld classifier. 
%
% When no output argument is specified, or when FID is given, the
% confusion matrix is displayed or written a to a text file. It is assumed
% that LAB1 contains true labels and LAB2 stores estimated labels.
%
% EXAMPLE
% Typical use of CONFMAT is the comparison of true and and estimated labels
% of a testset A by application to a trained classifier W: 
% LAB1 = GETLABELS(A); LAB2 = A*W*LABELD.
% More examples can be found in PREX_CONFMAT, PREX_MATCHLAB.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, GETLABELS, LABELD

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: confmat.m,v 1.8 2010/06/01 08:45:31 duin Exp $

function [CC,ne,lablist1,lablist2] = confmat (arg1,arg2,arg3,fid)

		% Check arguments.
  if nargin < 4, fid = 1; end
	if nargin < 3 | isempty(arg3)
		if isdataset(arg1)
			if islabtype(arg1,'crisp')
				lablist1 = getlablist(arg1);
				nlab1 = getnlab(arg1);
				lab2 = arg1*labeld;
				nlab2 = renumlab(lab2,lablist1); % try to fit on original lablist
				if any(nlab2==0)                 % doesn't fit:  contruct its own one 
					[nlab2,lablist2] = renumlab(lab2);
				else
					lablist2 = lablist1;
				end
				if nargin < 2| isempty(arg2)
					method = 'count';
					prwarning(4,'no method supplied, assuming count');
				else
					method = arg2;
				end
			elseif islabtype(arg1,'soft')
				true_labs = gettargets(arg1);
				if any(true_labs > 1 | true_labs < 0)
					error('True soft labels should be between 0 and 1')
				end
				est_labs = +arg1;
				if any(est_labs > 1 | est_labs < 0)
					error('Estimated soft labels should be between 0 and 1')
				end
				C = true_labs'*est_labs;
				output(C,getlablist(arg1),getfeatlab(arg1),'soft',fid);
				return
			else
				error('Confusion matrix can only be computed for crisp and soft labeld datasets')
			end
		else
			method = 'count';
			prwarning(4,'no method supplied, assuming count');
			if (nargin < 2 | isempty(arg2))
				error('Second label list not supplied')
			end
			[nlab1,lablist1] = renumlab(arg1);
			[nlab2,lablist2] = renumlab(arg2);
		end
	else
		[nlab1,lablist1] = renumlab(arg1);
		[nlab2,lablist2] = renumlab(arg2);
		%lab1 = arg1;
		%lab2 = arg2;
		method = arg3;
	end
	if nargin < 2
		if ~isdataset(arg1)
			error('two labellists or one dataset should be supplied')
		end
	end
		
	% Renumber LAB1 and LAB2 and find number of unique labels.

	m = size(nlab1,1);
	if (m~=size(nlab2,1))
		error('LAB1 and LAB2 have to have the same lengths.');
	end
	n1 = size(lablist1,1);
	n2 = size(lablist2,1);
	n = max(n1,n2); 
	
	% Construct matrix of co-occurences (confusion matrix).

  if strcmp(method,'ids')
    C = cell(n1+1,n2+1); % we need to store object id's
  else
    C = zeros(n1+1,n2+1);
  end
	for i = 0:n1
		K = find(nlab1==i);
		if (isempty(K))
      if strcmp(method,'ids')
        C(i+1,:) = repmat({[]},1,n2+1);
      else
        C(i+1,:) = zeros(1,n2+1);
      end
		else
			for j = 0:n2
        if strcmp(method,'ids')
          C(i+1,j+1) = {K(nlab2(K)==j)};
        else
          C(i+1,j+1) = sum(nlab2(K)==j);
        end
			end
		end
	end

	% position rejects and unlabeled object at the end of the matrix
	
	D = C;
	D(1:end-1,1:end-1) = C(2:end,2:end);
	D(end,:) = [C(1,2:end) C(1,1)];
	D(1:end-1,end) = C(2:end,1);
	C = D;
	D = D(1:end-1,1:end-1);
	DD = D(1:min(n1,n2),1:min(n1,n2));
	
	% Calculate number of errors ('count') or disagreement ('disagreement').
	% Neglect rejects
  
	switch (method)
		case {'count','ids'}
			J = find(nlab1~=0 & nlab2~=0);
			ne = nlabcmp(lablist1(nlab1(J),:),lablist2(nlab2(J),:));																															
		case 'disagreement'
			ne = (sum(sum(D)) - sum(diag(DD)))/m;  % Relative sum of off-diagonal    
			ne = ne/m;														 % entries.
			E = repmat(sum(D,2),1,n2);             % Disagreement = 1 - 
			D = ones(n1,n2)-D./E;                  % relative co-occurence.
			D = D / (n-1);
		otherwise
			error('unknown method');
  end

  if ~strcmp(method,'ids')
    %Distinguish 'rejects / no_labels' from 'non_rejects / fully labeled'
    if (any(C(:,end) ~= 0) | any(C(end,:)~=0)) & strcmp(method,'count')
      n1 = n1+1; n2 = n2+1;
      labch = char(strlab(lablist2),'reject');
      labcv = char(strlab(lablist1),'No');
    else
      labch = strlab(lablist2);
      labcv = strlab(lablist1);
      %labcv = labch;
      C = D;
    end
  end

    % If no output argument is specified, pretty-print C.

	if ((nargout == 0) || nargin == 4) && (~strcmp(method,'ids'))

		output(C,labcv,labch,method,fid)
	end
	
	if nargout > 0 || strcmp(method,'ids')
		CC = C;
	end
	
return

function output(C,labcv,labch,method,fid)
	
	n1 = size(labcv,1);
	n2 = size(labch,1);
		
		% Make sure labels are stored in LABC as matrix of characters, 
    % max. 6 per label.

		if (size(labch,2) > 6)
			labch = labch(:,1:6); 
			%labcv = labcv(:,1:6); 
		end
		if (size(labch,2) < 5)
			labch = [labch repmat(' ',n2,ceil((5-size(labch,2))/2))]; 
			labcv = [labcv repmat(' ',n1,ceil((5-size(labcv,2))/2))]; 
		end

		%C = round(1000*C./repmat(sum(C,2),1,size(C,2)));
		
		nspace = max(size(labcv,2)-7,0);
		cspace = repmat(' ',1,nspace);
		%fprintf(fid,['\n' cspace '        | Estimated Labels']);
		fprintf(fid,['\n  True   ' cspace '| Estimated Labels']);
		fprintf(fid,['\n  Labels ' cspace '|']);
		for j = 1:n2, fprintf(fid,'%7s',labch(j,:)); end
		fprintf(fid,'|');
		fprintf(fid,' Totals');
		fprintf(fid,'\n ');
		fprintf(fid,repmat('-',1,8+nspace));
		fprintf(fid,'|%s',repmat('-',1,7*n2));
		fprintf(fid,'|-------');
		fprintf(fid,'\n ');
	
		for j = 1:n1
			fprintf(fid,' %-7s|',labcv(j,:));
			switch (method)
				case 'count'
					fprintf(fid,'%5i  ',C(j,:)');
					fprintf(fid,'|');
					fprintf(fid,'%5i',sum(C(j,:)));
				case 'disagreement'
					fprintf(fid,' %5.3f ',C(j,:)');
					fprintf(fid,'|');
					fprintf(fid,' %5.3f ',sum(C(j,:)));
				case 'soft'
					fprintf(fid,' %5.2f ',C(j,:)');
					fprintf(fid,'|');
					fprintf(fid,' %5.2f ',sum(C(j,:)));
			end
			fprintf(fid,'\n ');
		end

		fprintf(fid,repmat('-',1,8+nspace));
		fprintf(fid,'|%s',repmat('-',1,7*n2));
		fprintf(fid,'|-------');
		fprintf(fid,['\n  Totals ' cspace '|']);

		switch (method)
			case 'count'
				fprintf(fid,'%5i  ',sum(C));
				fprintf(fid,'|');
				fprintf(fid,'%5i',sum(C(:)));
			case 'disagreement'
				fprintf(fid,' %5.3f ',sum(C));
				fprintf(fid,'|');
			%	fprintf(fid,' %5.3f ',sum(C(:)));
			case 'soft'
				fprintf(fid,' %5.2f ',C(j,:)');
				fprintf(fid,'|');
				fprintf(fid,' %5.2f ',sum(C(:)));
		end
		fprintf(fid,'\n\n');
	return
