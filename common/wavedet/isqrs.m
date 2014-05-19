function [ann]=isqrs(ann,heasig,t)
% isqrs.m is a function for getting only QRS marks from annotator
% Syntaxis:
%       function [ann]=isqrs(ann,heasig,t)
% Input parameters:
%  	ann: annotation structure to analyze
%	heasig: signal header structure
%	t=[t0 t1]: time interval vector
% Output parameter:
%	ann: new annotation structure

% Salvador Olmos
% e-mail: olmos@posta.unizar.es


% Removing non-QRS annotations

aux=find(ann.anntyp~='N' & ann.anntyp~='A' & ann.anntyp~='V' & ...
         ann.anntyp~='L' & ann.anntyp~='R' & ann.anntyp~='J' & ...
         ann.anntyp~='F' & ann.anntyp~='S' & ann.anntyp~='j' & ...
         ann.anntyp~='e' & ann.anntyp~='a' & ann.anntyp~='r' & ...
         ann.anntyp~='/' & ann.anntyp~='E' & ann.anntyp~='f' & ...
         ann.anntyp~='Q' & ann.anntyp~='?' & ann.anntyp~='B' & ...
	 ann.anntyp~='n' & ann.anntyp~='P');

ann.anntyp(aux)=[];
ann.time(aux)=[];
ann.num(aux)=[];
ann.subtyp(aux)=[];
ann.chan(aux)=[];
ann.aux(aux,:)=[];

if nargin > 1
    % Removing QRS annotations outside analysis time interval t

    aux=find(ann.time<t(1) | ann.time>(t(2)+2*heasig.freq) );
    if ~isempty(t)
        ann.anntyp(aux)=[];
        ann.time(aux)=[];
        ann.num(aux)=[];
        ann.subtyp(aux)=[];
        ann.chan(aux)=[];
        ann.aux(aux,:)=[];
    end
end
