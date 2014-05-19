function clust_labels = cluster_data_em(data2clust, CantClusters, iter_times)

if nargin < 3 || isempty(iter_times)
    iter_times = 1;
end
if nargin < 2 || isempty(CantClusters)
    CantClusters = 5;
end

cant_ejemplos = getsize(data2clust,1);

std_train = std(+data2clust);
data2clust = setdata(data2clust, bsxfun(@rdivide, +data2clust, std_train ));

clustered_labels = repmat('0', cant_ejemplos, iter_times);
jj = 1;
%         for CantClusters = 2:5
for jj = 1:iter_times

    %EMclust
    ii = 0;
    bContinuar = true;
    while(bContinuar && ii < CantClusters)
        try
            dbclear if caught error            
            clust_labels = em_loop(data2clust, CantClusters-ii );
            dbstop if caught error            
            bContinuar = false;
        catch ME
            if(strcmpi(ME.message, 'Not possible to find desired number of components'))
                ii = ii + 1;
            else
                rethrow(ME)
            end
        end
    end
    clustered_labels(:,jj) = char(97+clust_labels);
%             jj = jj + 1;
end

bClustered = true;
while(bClustered)
    %analizo las coincidencias en las distintas iteraciones respecto data2clust
    %los etiquetados.
    [~, sort_idx] = sort( cellstr(clustered_labels) );
    [all_clusters, aux_location] = unique(clustered_labels(sort_idx,:), 'rows', 'first');

    aux_location = [colvec(aux_location); cant_ejemplos+1];
    cluster_sizes = diff(aux_location);
    [~, clust_sort_idx] = sort(cluster_sizes, 'descend');

    %agrupo todos los subclusters data2clust una dist maxima max_distance
    max_distance = round(0.2*iter_times); %distancia maxima para considerarse un cluster
    cant_clusters = size(all_clusters,1);
    cant_iter = size(all_clusters,2);
    aux_1 = repmat( all_clusters, cant_clusters, 1);
    aux_idx = colvec(repmat(1:cant_clusters,cant_clusters,1));
    aux_2 = all_clusters(aux_idx,:);
    distances = reshape( sum(aux_1 ~= aux_2,2), cant_clusters, cant_clusters );
    remaining_clusters = 1:cant_clusters;

    %clusterizo igual todo los clusters iguales hasta max_distance
    for ii = rowvec(clust_sort_idx)
        bClustered = false;
        cluster2fusion_idx = find(distances(:,ii) <= max_distance);
        cluster2fusion_idx = cluster2fusion_idx(cluster2fusion_idx < ii | cluster2fusion_idx > ii );
        for jj = 1:length(cluster2fusion_idx)
            aux_idx = find(strcmpi(cellstr(clustered_labels), cellstr(all_clusters(cluster2fusion_idx(jj),:))));
            if( ~isempty(aux_idx) )
                bClustered = true;
                remaining_clusters( remaining_clusters == cluster2fusion_idx(jj)) = [];
                clustered_labels( aux_idx ,:) = repmat(all_clusters(ii,:), length(aux_idx), 1 );
                all_clusters(cluster2fusion_idx(jj),:) = all_clusters(ii,:);
            end
        end
        if(bClustered)
            remaining_clusters( remaining_clusters == ii ) = [];
            %fuerzo el recalclo de distancias, para que no se junte todo y solo
            %los clusters grandes se coman data2clust los mas chicos.
            break
        end
    end
end

%REanalizo las coincidencias en las distintas iteraciones respecto data2clust
%los etiquetados.
[~, sort_idx] = sort( cellstr(clustered_labels) );
[all_clusters, aux_location] = unique(clustered_labels(sort_idx,:), 'rows', 'first');

aux_location = [colvec(aux_location); cant_ejemplos+1];
cluster_sizes = diff(aux_location);

cant_clusters = size(all_clusters,1);

clust_labels = nan(cant_ejemplos,1);
for ii = 1:cant_clusters
    clust_labels(sort_idx(aux_location(ii):(aux_location(ii+1)-1))) = ii;
end

if(any(isnan(clust_labels)))
   error() 
end



function new_lab = em_loop(data2clust, CantClusters )

n_ini		= 500;			% Maximum size of subset to use for initialisation.

cant_ejemplos = size(data2clust,1);

not_found = true;
iter = 0;
while(not_found)
    
    % try to find an initialisation with all class sizes > 1

    if (cant_ejemplos > n_ini)       						% ... data2clust random subset of data2clust.
        %sampleo uniformemente el feature-space
        a_ini = data2clust;
        a_ini = sum(abs(a_ini),2);
        [~, sort_idx ] = sort(a_ini);
        a_ini = data2clust(sort_idx(randsample(cant_ejemplos, n_ini)),:);
    else
        a_ini = data2clust;								% ... the entire set data2clust.
    end

    iter = iter + 1;
    if iter > 3
        error('Not possible to find desired number of components')
    end
    % 50 trials
    assign = zeros(1,CantClusters-1);
    std_a_ini = std(a_ini);
    while( length(unique(assign)) ~= CantClusters )
        % add some noise to data to avoid problems
        a_ini = a_ini.*(ones(size(a_ini))+ bsxfun(@times, max( [repmat(0.001,size(std_a_ini)); 0.01*std_a_ini] ), randn(size(a_ini))));
        %construct the distance matrix
        % The order of operations below is good for the accuracy.
        dist_mat = ones(n_ini,1)*sum(a_ini'.^2,1);
        dist_mat = dist_mat + sum(a_ini.^2,2)*ones(1,n_ini);
        dist_mat = dist_mat - 2 .* (a_ini)*(a_ini)';

        % Check for a numerical inaccuracy. 
        dist_mat(dist_mat<0) = zeros(size(J));          % dist_mat should be nonnegative.

        % take care of symmetric distance matrix
        dist_mat(1:n_ini+1:n_ini*n_ini) = zeros(1,n_ini);
        
        try
            assign  = kcentres(dist_mat,CantClusters,50);
        catch ME
            error('Not possible to find desired number of components')
        end
    end
    % Train initial classifier on labels generated by KCENTRES and find
    % initial hard labels. Use NMC instead of W_CLUST to make sure that we 
    % always have enough data to estimate the parameters.
    a_ini = dataset(a_ini,assign); 
    a_ini = setprior(a_ini,getprior(a_ini,0));
    dist_mat = data2clust*(a_ini*nmc);
    
    new_lab = dist_mat*labeld;
    cs_new_lab = classsizes(dataset(dist_mat,new_lab));
    if length(cs_new_lab) == CantClusters && all( cs_new_lab > 1)
        not_found = 0;
    end
end

lablist_org = [];

% Ready for the work.
iter = 0;
change = 1;
lab = ones(cant_ejemplos,1);

while (change ~= 0 && any(lab ~= new_lab)) % EM loop, run until labeling is stable.

    prwaitbar(100,100-100*exp(-iter/10));
    data2clust = setlabels(data2clust,new_lab); 		% 0. Set labels and store old labels.
    data2clust = setprior(data2clust,getprior(data2clust,0));%    Set priors to class frequencies
    lab = new_lab;								% 
    data2clust = remclass(data2clust,1);            %    demand class sizes > 2 objects
    iter = 0;

    while getsize(data2clust,3) < CantClusters        %    increase number of classes if necessary

        iter = iter + 1;

        if iter > 5
            error('Not possible to find desired number of components')
        end
        laba = getlablist(data2clust);
        labmax = max(laba);
        CantClusters = classsizes(data2clust);
        [Nmax,cmax] = max(CantClusters);        % find largest class
        aa = seldat(data2clust,cmax);         % select just that one

        std_aa = std(+aa);
        bContinue = true;
        while( bContinue )
            try
                new_lab_aa = kmeans(aa .* (ones(size(+aa))+ bsxfun(@times, max( [repmat(0.001,size(+aa)); 0.01*std_aa] ), randn(size(+aa)))) ,2);   % split it by kmeans
                bContinue = false;
            catch ME
            %                         if( strcmpi(ME.message, 'kcentres fails as some objects are identical: add some noise'))
            %                             std_aa = std_aa * 2;
            %                         else
            %                             rethrow(ME)
            %                         end
                error('Not possible to find desired number of components')
            end
        end

        N1 = sum(new_lab_aa == 1);   
        N2 = sum(new_lab_aa == 2);
        if (N1 > 1 & N2 > 1) % use it if both classes have more than one sample
            J = findlabels(data2clust,laba(cmax,:));
            data2clust = setlabels(data2clust,new_lab_aa + labmax,J);
        end
    end

    std_a = std(+data2clust);
    w_em = (data2clust .* (ones(size(+data2clust))+ bsxfun(@times, max( [repmat(0.001,size(+data2clust)); 0.01*std_a] ), randn(size(+data2clust))))) * w_clust;             % 1. Compute classifier, crisp output.
    b = data2clust*w_em;               		% 2. Classify training samples.
    new_lab = labeld(b);      		% 3. Insert classification as new labels.
    iter = iter+1;             %DXD Added also the iter for the crisp labels
    
    if iter > 50
        change = 0;
    end

    %para ver la velocidad de convergencia.
    %             disp(sum(lab ~= new_lab))

end


if ~isempty(lablist_org) % substitute original labels if desired
    new_lab = lablist_org(new_lab);
    wlab = getlabels(w_em);
    wlab = lablist_org(wlab);
    w_em = setlabels(w_em,wlab);
end


function labels = kcentres( dist_mat, CantClusters, itern)

dmax = max(max(dist_mat));
dopt = inf;

for tri = 1:itern

    M = randperm(cant_ejemplos); M = M(1:CantClusters);		% Random initializations

    J = zeros(1,CantClusters); 					    	  % Center objects to be found.

    % Iterate until J == M. See below.
    while true
        
        [dm,I] = min(dist_mat(M,:));

        % Find CantClusters centers. 
        for ii = 1:CantClusters                   
            JJ = find(I==ii); 
            if (isempty(JJ))
                %JJ can be empty if two or more objects are in the same position of 
                % feature space in dataset
                J(ii) = 0;									
            else
                % Find objects closest to the object M(ii)
                [dummy,j,dm] = kcentres(dist_mat(JJ,JJ),1,1);
                J(ii) = JJ(j);
            end
        end

        J(find(J==0)) = [];
        CantClusters = length(J);
        if CantClusters == 0
            error('kcentres fails as some objects are identical: add some noise')
        end
        if (length(M) == CantClusters) & (all(M == J))
            % CantClusters centers are found.
            [dmin,labs] = min(dist_mat(J,:));
            dmin = max(dmin);
            break;
        end
        M = J;
    end

    % Store the optimal centers in JOPT.

    if (dmin <= dopt)
        dopt = dmin;
        labels = labs';
        Jopt = J;
    end

end

prwaitbar(0)

% Determine the best centers over the N trials.
dm = zeros(1,CantClusters);   
for ii=1:CantClusters
    L = find(labels==ii);
    dm(ii) = max(dist_mat(Jopt(ii),L));
end


