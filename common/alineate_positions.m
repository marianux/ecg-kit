%% (Internal) Alineate wave positions (onset, peak and offsets commonly) 
%
%     [unique_elements idx ] = alineate_positions( {onset, peak, offset}, max_tol )
% 
% 
% Arguments:
% 
%   + aux_mat: data elements, 
% 
%   + win_size: tolerance to consider aux_mat(i) == aux_mat(i+1).
% 
% Output:
% 
%   + sft_intersect1: the elements in the soft intersection
% 
%   + idx1, idx2: the indexes of aux_mat(idx1) which are in the soft
%   intersection.
% 
% Example:
% 
% 
% See also soft_set_intersect
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate: 14/05/2016
% Last update: 14/05/2016
% Copyright 2008-2016
% 
function varargout = alineate_positions( val1, max_tol )

if( nargin < 2 || isempty(max_tol) )
    max_tol = realmax;
end

sizes = cellfun( @(a)(length(a)), val1 );
can_val1 = length(sizes);
max_size = max(sizes);

nnargout = nargout;
if nnargout ~= can_val1
    error('alineate_positions:OutputError', 'You must expect the same amount of arguments you provide. %d arguments provided, %d expected.\n', can_val1, nnargout);
end

% indexes of can_val1 pairwise distances 
aux_idx = colvec(1:sizes(end));
jj = 1;
for ii = can_val1-1:-1:1
    aux_idx = [ colvec(repmat(colvec(1:max_size),1,max_size^jj)') repmat(aux_idx,max_size,1) ];
    jj = jj + 1;
end

% all pairwise distances indexes
dis_mat_idx = nchoosek(1:can_val1,2);
can_pw = size(dis_mat_idx,1);

% index matrix to address the distance matrix corresponding to pairwise
% comparison
dis_mat2_idx = nan(can_pw);
for ii = 1:can_pw
    dis_mat2_idx(dis_mat_idx(ii,1), dis_mat_idx(ii,2)) = ii;
    dis_mat2_idx(dis_mat_idx(ii,2), dis_mat_idx(ii,1)) = ii;
end

aux_mat = cell2mat(arrayfun( @(a)( [val1{a}; nan(max_size-length(val1{a}),1)]) , 1:can_val1, 'UniformOutput', false));
dis_mat = cell2mat(arrayfun( @(a,b)( abs( repmat( aux_mat(:,a) , 1, max_size) - repmat(aux_mat(:,b), 1, max_size)' ) ), reshape(dis_mat_idx(:,1),1,1,size(dis_mat_idx,1) ), reshape(dis_mat_idx(:,2),1,1,size(dis_mat_idx,1) ), 'UniformOutput', false));

total_comp = max_size^can_val1;
dist_array = nan(max_size^can_val1,1);

for ii = 1:total_comp
    dist_array(ii) = dis_mat(aux_idx(ii,1), aux_idx(ii,end), dis_mat2_idx(1, can_val1) );
    for jj = 2:can_pw
        dist_array(ii) = dist_array(ii) + dis_mat(aux_idx(ii,jj-1), aux_idx(ii,jj), dis_mat2_idx(jj-1, jj) );
    end
end

[~, s_idx] = sort(dist_array);

aligned_pos = nan(max_size,can_val1);
for ii = 1:max_size
%     [] = ind2sub([max_size,max_size,can_val1], s_idx(ii));
%     aligned_pos(ii) = ;
    aligned_pos(ii,:) = arrayfun( @(a)(aux_mat(aux_idx(s_idx(ii),a),a)), 1:can_val1 );
end

[~, s_idx] = sort(aligned_pos(:,1));
aligned_pos = aligned_pos(s_idx,:);

for ii = 1:nnargout
   varargout{ii} = aligned_pos(:,ii);
end
