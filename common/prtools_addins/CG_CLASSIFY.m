% Version 1.000
%
% Code provided by Ruslan Salakhutdinov and Geoff Hinton
%
% Permission is granted for anyone to copy, use, modify, or distribute this
% program and accompanying programs and documents for any purpose, provided
% this copyright notice is retained and prominently displayed, along with
% a note saying that the original programs are available from our
% web page.
% The programs and documents are distributed without any warranty, express or
% implied.  As the programs were written for research purposes only, they have
% not been tested to the degree that would be advisable in any important
% application.  All use of these programs is entirely at the user's own risk.


function [f, df] = CG_CLASSIFY(VV,Dim,XX,target)

% l1 = Dim(1);
% l2 = Dim(2);
% l3= Dim(3);
% l4= Dim(4);
% l5= Dim(5);

N = size(XX,1);
hid_num = length(Dim)-2;

% Do decomversion.

xxx = 0;
for ii=1:hid_num
    model{ii}.W = reshape(VV(xxx+1:xxx+(Dim(ii)+1)*Dim(ii+1)),Dim(ii)+1,Dim(ii+1));
    xxx = xxx+(Dim(ii)+1)*Dim(ii+1);
end
w_class = reshape(VV(xxx+1:xxx+(Dim(ii+1)+1)*Dim(ii+2)),Dim(ii+1)+1,Dim(ii+2));

%  w1 = reshape(VV(1:(l1+1)*l2),l1+1,l2);
%  xxx = (l1+1)*l2;
%  w2 = reshape(VV(xxx+1:xxx+(l2+1)*l3),l2+1,l3);
%  xxx = xxx+(l2+1)*l3;
%  w3 = reshape(VV(xxx+1:xxx+(l3+1)*l4),l3+1,l4);
%  xxx = xxx+(l3+1)*l4;
%  w_class = reshape(VV(xxx+1:xxx+(l4+1)*l5),l4+1,l5);


%   XX = [XX ones(N,1)];
%   w1probs = 1./(1 + exp(-XX*w1)); 
%   w1probs = [w1probs  ones(N,1)];
%   w2probs = 1./(1 + exp(-w1probs*w2)); 
%   w2probs = [w2probs ones(N,1)];
%   w3probs = 1./(1 + exp(-w2probs*w3)); 
%   w3probs = [w3probs  ones(N,1)];

XX = [XX ones(N,1)];
wxprobs = cell(hid_num,1);
wxprobs{1} = 1./(1 + exp(-XX * model{1}.W)); 
wxprobs{1} = [wxprobs{1} ones(N,1)];
for ii=2:hid_num
    wxprobs{ii} = 1./(1 + exp(-wxprobs{ii-1} * model{ii}.W)); 
    wxprobs{ii} = [wxprobs{ii}  ones(N,1)];
end

% w1probs = wxprobs{1};
% w2probs = wxprobs{2};
% w3probs = wxprobs{3};
% w2 = model{2}.W;
% w3 = model{3}.W;

% targetout = exp(w3probs*w_class);
targetout = exp(wxprobs{hid_num}*w_class);
targetout = targetout./repmat(sum(targetout,2),1,size(target,2));
f = -sum(sum( target .* log(targetout))) ;

%     [~, J]=max(targetout,[],2);
%     [~, J1]=max(target,[],2);
%     errors = sum(J~=J1);
%     fprintf(1, 'Error rate %3.2f\n', errors/N);
  
% IO = (targetout-target);
% Ix_class=IO; 
% dw_class =  w3probs'*Ix_class; 
% 
% Ix3 = (Ix_class*w_class').*w3probs.*(1-w3probs);
% Ix3 = Ix3(:,1:end-1);
% dw3 =  w2probs'*Ix3;
% 
% Ix2 = (Ix3*w3').*w2probs.*(1-w2probs); 
% Ix2 = Ix2(:,1:end-1);
% dw2 =  w1probs'*Ix2;
% 
% Ix1 = (Ix2*w2').*w1probs.*(1-w1probs); 
% Ix1 = Ix1(:,1:end-1);
% dw1 =  XX'*Ix1;
% 
% df = [dw1(:)' dw2(:)' dw3(:)' dw_class(:)']'; 

Ix = (targetout-target);
df = colvec(wxprobs{hid_num}' * Ix);

Ix = (Ix*w_class') .* wxprobs{hid_num} .* (1-wxprobs{hid_num});    
Ix = Ix(:,1:end-1);
df = [ colvec( wxprobs{hid_num-1}'*Ix ) ; df];

if( hid_num > 2 )
    for ii=hid_num-1:-1:2

        Ix = (Ix * model{ii+1}.W') .* wxprobs{ii} .* (1-wxprobs{ii});    
        Ix = Ix(:,1:end-1);
        df = [ colvec( wxprobs{ii-1}'*Ix ) ; df];

    end
end

Ix = (Ix * model{2}.W') .* wxprobs{1} .* (1-wxprobs{1});    
Ix = Ix(:,1:end-1);
df = [ colvec( XX'*Ix ) ; df];
