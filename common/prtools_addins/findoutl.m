%FINDOUTL Find outliers from a dataset
%
%  B = FINDOUTL(A,T,P)
%  B = A*FINDOUTL([],T,P)
%
% INPUT
%  A  Dataset
%  P  [0 to 1] Maximum fraction of outliers to be detected. In case more outliers than
%     P would be detected, then only the worst P fraction according to wi
%     score are considered as outliers.
%
% OUTPUT
%  B  Dataset
%
% DESCRIPTION
% Outliers in A are removed, other objects are copied to B. Class by class
% a distance matrix is constructed and objects are removed that have a fraction
% P of their distances larger than the average distance in the class + T times
% the standard deviation of the within-class distances. This routine works
% also on unlabeled datasets. In partially labeled datasets the unlabeled
% objects are neglected.

function R = findoutl(a,p);
if nargin < 2, p = []; end
if nargin < 1 | isempty(a)
	a = mapping(mfilename,'fixed',{t,p});
	a = setname(a,'rem_outliers');
	return
end

cs = classsizes(a);
c = sum(cs > 0);

if c == 0
    [E,V] = eig(covm(a)); 
    aa = (+a)*E*sqrt(inv(V));
	d = sqrt(distm(aa));
	J = findoutd(d,t,p);
	a(J,:) = [];
else
    classes_idx = unique(getnlab(a));
    classes_idx = (classes_idx(:))';
	R = [];
	for j = classes_idx
		
        L = findnlab(a,j);

% %         Primer forma en que lo hice basado en los algoritmos del PRtools.
% %         Tiene la desventaja que la presencia de outliers perjudica mucho
% %         las estimaciones de media y covarianza y puede no encontrar bien
% %         los outliers.
% 
%         %Lo primero que hacemos es asegurarnos de "esferizar" las
%         %distribuciones para compatibilizar las varianzas de cada
%         %caracterísitica.
%         [E,V] = eig(covm(a(L,:))); 
%         aa = (+a(L,:))*E*sqrt(inv(V));
% 		
%         %Luego una vez esferizada es válida la distancia euclidea como
%         %medida de distancia.
%         d = sqrt(distm(aa));
% 		J = findoutd(d,t,p);
% 		R = [R;L(J)];

%         Segunda alternativa basada en “Robust covariance matrix
%         estimation and multivariate outlier detection” (con F. J. Prieto). Artículo 
%         En este articulo se propone una remoción de outliers proyectando
%         en las direcciones de maxima (provocada por outliers) y minima
%         (provocada por la multimodalida) curtosis y estandarizando
%         robustamente (con medianas y dispercion, ver Capítulo 4: GRAPHICAL
%         ANALYSIS AND OUTLIERS de Daniel Peña) las variables para la
%         remocion.

%         [idx,dm,meanData,CovData] = kur_rce(data,-1);

        %El parámetro -1 indica que solo iterará una vez sobre la
        %información y en la dirección de máxima curtosis, que es la que
        %suele indicar la presencia de outliers.

%         try
            
            %En caso de tener colinearidad, descarto esa columna para el
            %análisis de outliers.
            matrix_slice = +a(L,:);

%             linear_indep_cols = rank(matrix_slice);
% 
%             %Busco la variable con alta correlacion
%             if( linear_indep_cols < k )
% 
%                 cols_subset = 1:k;
%                 CCmatrix_slice = corrcoef(matrix_slice);
%                 %elimino la diagonal.
%                 CCmatrix_slice(1:k+1:k*k) = 0;
%                 [row, col] = find(CCmatrix_slice>=1-eps);
%                 subset_size = k;
% 
%                 while( ~isempty(row) )
% 
%                     colinear_cols = unique(row);
%                     cols_subset(cols_subset == colinear_cols(end) ) = [];
%                     subset_size = subset_size - 1;
%                     matrix_slice = matrix_slice(:,cols_subset);
% 
%                     CCmatrix_slice = corrcoef(matrix_slice);
%                     %elimino la diagonal.
%                     CCmatrix_slice(1:subset_size+1:subset_size*subset_size) = 0;
%                     [row, col] = find(CCmatrix_slice>=1-eps);
% 
%                 end
% 
%             end
% 
%             idx = kur_rce( matrix_slice, -1);
            
            [idx wi] = pcout(matrix_slice);            
            
            lmatrix_slice = size(matrix_slice,1);
            
            if( ~isempty(p) && sum(~idx)/lmatrix_slice > p)
                
                [dummy ind] = sort(wi);
                
                idx = idx | 1;
                
                idx(ind(1:floor(p*lmatrix_slice))) = false;
                
            end
            
            idx = ~idx;
            
%         catch
%             %Si hubo algun problema con la busqueda de outliers lo atajo, y
%             %aborto la detección para esta clase.
%             warning('Problemas con kur_rce! No se quitan outliers.');
% 
%             error = lasterror;
%             disp(error.message);
%             for i = 1:length(error.stack)
%                 disp(error.stack(i).file);
%                 disp(error.stack(i).name);
%                 disp(error.stack(i).line);
%                 disp('-------------------')
%             end
%             
%             idx = false( length(L),1 );
%         end
        
        R = [R;L(find(idx))];
        
	end
end

%FINDOUTD Detect outliers in distance matrix
%
% J = FINDOUTD(D,T,P)
%
% Find the indices of the objects in the dissimilarity representation D
% that have a fraction P (default P = 0.10) of their distances larger than
% mean(D(:)) + T * std(D), default T = 3.

function J = findoutd(d,t,p);

if nargin < 3 | isempty(p), p = 0.10; end
if nargin < 2 | isempty(t), t = 3; end
x = +d;
x = x(:);
L = find(x~=0);
x = x(L);
s = mean(x) + t * std(x);
J = find(sum(+d > s,2) > p*size(d,2));


