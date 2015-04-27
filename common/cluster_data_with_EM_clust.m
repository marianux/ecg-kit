%% (Internal) Cluster data with expectation-maximization algorithm
%   
%   [clust_labels w_Trained_Classifier] = cluster_data_with_EM_clust(dsTrain, w_mapp, CantClusters, iter_times)
% 
% Arguments:
% 
%      + dsTrain: PRdataset with data
% 
%      + w_mapp: mapping to perform clustering.
% 
%      + CantClusters: clusters to discover.
% 
%      + iter_times: iterations to perform clustering.
% 
% Output:
% 
%      + clust_labels: Cluster pertenence labels
% 
%      + w_Trained_Classifier: Mapping that perform the clustering found.
% 
% Example:
% 
%             warning off all; prwarning(0);
%             Clust_Labels = cluster_data_with_EM_clust(featMat_clust(pending_hb_idx,:), qdc_new([],1e-6,1e-6, []), CantClusters, iter_times);
%             warning on all; prwarning(1);
% 
% See also a2hbc_main
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function [clust_labels w_Trained_Classifier] = cluster_data_with_EM_clust(dsTrain, w_mapp, CantClusters, iter_times)

if nargin < 4 || isempty(iter_times)
    iter_times = 1;
end
if nargin < 3 || isempty(CantClusters)
    CantClusters = 5;
end

m = getsize(dsTrain,1);

std_train = std(+dsTrain);
dsTrain = setdata(dsTrain, bsxfun(@rdivide, +dsTrain, std_train ));

clustered_labels = repmat('0', m, iter_times);
jj = 1;
%         for CantClusters = 2:5
for jj = 1:iter_times

    %EMclust
    ii = 0;
    bContinuar = true;
    while(bContinuar && ii < CantClusters)
        try
%             dbclear if caught error            
%             [clust_labels w_Trained_Classifier ] = emclust_new(dsTrain, w_mapp, CantClusters-ii );
            [clust_labels w_Trained_Classifier ] = emclust(dsTrain, w_mapp, CantClusters-ii );
%             dbstop if caught error            
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
    %analizo las coincidencias en las distintas iteraciones respecto a
    %los etiquetados.
    [~, sort_idx] = sort( cellstr(clustered_labels) );
    [all_clusters, aux_location] = unique(clustered_labels(sort_idx,:), 'rows', 'first');

    aux_location = [colvec(aux_location); m+1];
    cluster_sizes = diff(aux_location);
    [~, clust_sort_idx] = sort(cluster_sizes, 'descend');

    %agrupo todos los subclusters a una dist maxima max_distance
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
            %los clusters grandes se coman a los mas chicos.
            break
        end
    end
end

% %luego clusterizo todo lo que fue quedando a distancias mayores, al estilo
% %clustering jerarquico.
% kk = max_distance+1;
% while( kk <= iter_times)
%     for ii = remaining_clusters
%         bClustered = false;
%         cluster2fusion_idx = find(distances(:,ii) ==  kk );
%         cluster2fusion_idx = cluster2fusion_idx(cluster2fusion_idx > ii );
%         cluster2fusion_idx = intersect(cluster2fusion_idx, remaining_clusters);
%         
%         for jj = 1:length(cluster2fusion_idx)
%             aux_idx = find(strcmpi(cellstr(clustered_labels), cellstr(all_clusters(cluster2fusion_idx(jj),:))));
%             if( ~isempty(aux_idx) )
%                 bClustered = true;
%                 remaining_clusters( remaining_clusters == cluster2fusion_idx(jj)) = [];
%                 clustered_labels( aux_idx ,:) = repmat(all_clusters(ii,:), length(aux_idx), 1 );
%             end
%         end
%         if(bClustered)
%             remaining_clusters( remaining_clusters == ii ) = [];
%         end
%     end
%     kk = kk + 1;
% end

%REanalizo las coincidencias en las distintas iteraciones respecto a
%los etiquetados.
[~, sort_idx] = sort( cellstr(clustered_labels) );
[all_clusters, aux_location] = unique(clustered_labels(sort_idx,:), 'rows', 'first');

aux_location = [colvec(aux_location); m+1];
cluster_sizes = diff(aux_location);

% esto es por si interesaria filtrar clusters chicos, no parece buena idea
% porque agrupa juntas las clases que aparecen ocacionalmente, haciendolas
% indetectables.
% 
% big_clusters_start_idx = find(cluster_sizes > 50);
% big_clusters_range_idx = [ colvec(aux_location( big_clusters_start_idx )) colvec(aux_location( big_clusters_start_idx+1 )-1)];
% 
% big_clusters_idx = [];
% group_labels = [];
% cant_big_clusters = length(big_clusters_start_idx);
% for ii = 1:cant_big_clusters
%     big_clusters_idx = [ big_clusters_idx; colvec(big_clusters_range_idx(ii,1):big_clusters_range_idx(ii,2))];
%     group_labels = [group_labels ; repmat(ii, cluster_sizes(big_clusters_start_idx(ii)) ,1)];
% end
% 
% clust_labels = repmat(cant_big_clusters+1,m,1);
% cl    ust_labels(sort_idx(big_clusters_idx)) = group_labels;

cant_clusters = size(all_clusters,1);

clust_labels = nan(m,1);
for ii = 1:cant_clusters
    clust_labels(sort_idx(aux_location(ii):(aux_location(ii+1)-1))) = ii;
end

if(any(isnan(clust_labels)))
   error() 
end
