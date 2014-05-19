function ann = pos2ann(pos);

% function to convert the position structure to a MIT annotation structure
% Note: Q,R,S and R' individual peaks are lost in this conversion
% Changed on 31/07/02 Rute

types = '(p)(N)(tt)';
nbeats = length(pos.qrs);

matriz = zeros(3,nbeats*10);

matriz(1,:) = [pos.Pon pos.P pos.Poff pos.QRSon pos.qrs pos.QRSoff ...
              pos.Ton pos.T pos.Tprima pos.Toff];  % Time

matriz(2,:) = kron([1:10],ones(1,nbeats)) ;         % anntyp

matriz(3,:) = kron([0 0 0 1 0 1 2 0 0 2],ones(1,nbeats));
matriz(3,7*nbeats+1:9*nbeats) = [pos.Ttipo pos.Ttipo];

ind = find(isnan(matriz(1,:)));
matriz(:,ind)=[];

matriz = matriz';

matriz = sortrows(matriz,1);

ind = find(diff(matriz(:,1))<=0);

if ~isempty(ind),
   fprintf(1,'Warning: There is a repeated annotation at sample: %f tipo %d %d \n',matriz(ind,1:3)')

%%%%%Rute 31/07/02%%%%%%%% restrictions on the elimination of repeted elements %%%%%%%%%%%%%%
aux=(matriz(ind,2)~=1) & (matriz(ind,2)~=4) & (matriz(ind,2)~=7);
%aux=1;
%in these cases eliminate the element ind+1 otherwise eliminate the element ind %%%%%%%%%%%%%
matriz(ind+aux,:)=[]; %% Rute 31/07/02
end

ann.time = matriz(:,1);
ann.anntyp = types(matriz(:,2))';
ann.subtyp = blanks(size(matriz,1))';
ann.chan = char('0'*ones(size(matriz,1),1));
ann.num = char((matriz(:,3))+48);
ann.aux = char(zeros(size(matriz,1),1));
ann.aux (:,1)=[];
