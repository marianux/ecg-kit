function [ann]=isqrs(ann)
% isqrs.m is a function for removing non-WAVEFORM LIMITS marks from annotator
% function [ann]=iswfl(ann)
% Input parameter:
%	ann: annotation structure
% Output parameter:
%	ann: annotation structure with only QRS marks.


% ------------------------------------------
% Salvador Olmos:
% e-mail: olmos@posta.unizar.es
% Last modified: 20/Nov/1996
% ------------------------------------------



kk=find(ann.anntyp~='N' & ann.anntyp~='A' & ann.anntyp~='V' & ann.anntyp~='L' & ...
        ann.anntyp~='R' & ann.anntyp~='J' & ann.anntyp~='F' & ann.anntyp~='S' & ...
        ann.anntyp~='j' & ann.anntyp~='e' & ann.anntyp~='a' & ann.anntyp~='r' & ...
        ann.anntyp~='?' & ann.anntyp~='B' & ann.anntyp~='n' & ann.anntyp~='P' & ...
        ann.anntyp~='/' & ann.anntyp~='E' & ann.anntyp~='f' & ann.anntyp~='Q' & ...
        ann.anntyp~='(' & ann.anntyp~=')' & ann.anntyp~='p' & ann.anntyp~='t');


ann.anntyp(kk)=[];
ann.time(kk)=[];
ann.num(kk)=[];
ann.subtyp(kk)=[];
ann.chan(kk)=[];
ann.aux(kk,:)=[];

