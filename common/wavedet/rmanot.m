% function ann=rmanot(ann,vec) removes annotations in vec from annotation structure
% Input parameters:
%	ann: input annotations structure
%	vec: vector with records to remove
% Output parameters:
% 	ann: output annotation structure

% Salvador Olmos: olmos@posta.unizar.es  Last revised 7/3/98

function ann=rmanot(ann,vec)

ann.anntyp(vec)=[];
ann.time(vec)=[];
ann.subtyp(vec)=[];
ann.num(vec)=[];
ann.chan(vec)=[];
ann.aux(vec,:)=[];
