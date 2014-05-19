function POS=posmat(ann)
% function POS=posmat(ann) generates a Waveform Limit MATRIX with annotation ann
% Input parameters:
%	ann: annotation structure
% Output parameters:
%	POS: matrix with waveform limits

% Salvador Olmos. olmos@posta.unizar.es


ann=iswfl(ann);
q=isqrs(ann);

beats=length(q.time);

%%%%%%%%%
% Removing annotations before the first complete beat
%%%%%%%%%
ai=find(ann.anntyp=='('  & ann.num=='0' & ann.time<q.time(1)); % previous P onset
ann=rmanot(ann,ai(1:end-1));
ai=find(ann.anntyp=='p'  & ann.num=='0' & ann.time<q.time(1)); % Previous P
ann=rmanot(ann,ai(1:end-1));
ai=find(ann.anntyp==')'  & ann.num=='0' & ann.time<q.time(1)); % Previous P offset
ann=rmanot(ann,ai(1:end-1));

%%%%%%%%%
% FALTA ELIMINAR ANOTACIONES POSTERIORES AL ULTIMO LATIDO COMPLETO
%%%%%%%%%

POS=nan(beats,10);

POS(:,5)=q.time; 	% QRS fiducial point

ai=find(ann.anntyp=='(' & ann.num=='1');	% QRS onset

while (max(ann.time(ai))>max(POS(:,5))) 
    ai(end)=[]; 
end

if length(ai)==beats
   POS(:,4)=ann.time(ai);	
else
   for i=1:length(ai)
     p=min(find((ann.time(ai(i))-POS(:,5))<0));
     POS(p,4)=ann.time(ai(i));
   end
end

ai=find(ann.anntyp==')' & ann.num=='1');	% QRS offset
%while (max(ann.time(ai))>max(POS(:,5)+60)) ai(end)=[]; end
if length(ai)==beats
   POS(:,6)=ann.time(ai);	
else
   for i=1:length(ai)
     p=max(find((ann.time(ai(i))-POS(:,5))>0));
     POS(p,6)=ann.time(ai(i));
   end
end

ai=find(ann.anntyp=='(' & ann.num=='0');
while (max(ann.time(ai))>max(POS(:,4))) 
    ai(end)=[]; 
end
if length(ai)==beats
   POS(:,1)=ann.time(ai);	% P wave onset
else
   for i=1:length(ai)
     p=min(find((ann.time(ai(i))-POS(:,5))<0));
     POS(p,1)=ann.time(ai(i));
   end
end

ai=find(ann.anntyp=='p');
while (max(ann.time(ai))>max(POS(:,4)))
    ai(end)=[]; 
end
if length(ai)==beats
   POS(:,2)=ann.time(ai);	% P wave 
else
   for i=1:length(ai)
     p=min(find((ann.time(ai(i))-POS(:,5))<0));
     POS(p,2)=ann.time(ai(i));
   end
end

ai=find(ann.anntyp==')' & ann.num=='0');
while (max(ann.time(ai))>max(POS(:,4)))
    ai(end)=[]; 
end
if length(ai)==beats
   POS(:,3)=ann.time(ai);	% P wave offset
else
   for i=1:length(ai)
     p=min(find((ann.time(ai(i))-POS(:,5))<0));
     POS(p,3)=ann.time(ai(i));
   end
end

ai=find(ann.anntyp=='(' & ann.num=='2');
%while (max(ann.time(ai))<max(POS(:,4))) ai(end)=[]; end
if length(ai)==beats
   POS(:,7)=ann.time(ai);	% T wave onset
else
   for i=1:length(ai)
     p=max(find((POS(:,5)-ann.time(ai(i)))<0));
     POS(p,7)=ann.time(ai(i));
   end
end

ai=find(ann.anntyp==')' & ann.num=='2');
if length(ai)==beats
   POS(:,10)=ann.time(ai);	% T wave offset
else
   for i=1:length(ai)
     p=max(find((POS(:,5)-ann.time(ai(i)))<0));
     POS(p,10)=ann.time(ai(i));
   end
end

% T peak 
ai=find(ann.anntyp=='t');
for i=1:length(ai)
  p=max(find((POS(:,5)-ann.time(ai(i)))<0));
  if POS(p,8)==0 
      POS(p,8)=ann.time(ai(i));
  else
      POS(p,9)=ann.time(ai(i));
  end
end
