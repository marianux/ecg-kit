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
function [s, varargout] = alineate_positions( val1, max_tol )

if( nargin < 2 || isempty(max_tol) )
    max_tol = realmax;
end

sizes = cellfun( @(a)(length(a)), val1 );
can_val1 = length(sizes);
max_size = max(sizes);

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

dis_mat = cell2mat(arrayfun( @(a,b)( abs( repmat( [val1{a}; nan(max_size-length(val1{a}),1)] , 1, max_size) - repmat([val1{b}; nan(max_size-length(val1{b}),1)], 1, max_size)' ) ), reshape(dis_mat_idx(:,1),1,1,size(dis_mat_idx,1) ), reshape(dis_mat_idx(:,2),1,1,size(dis_mat_idx,1) ), 'UniformOutput', false));


total_comp = max_size^can_val1;
dist_array = zeros(total_comp,1);

for ii = 1:total_comp
    dist_array(ii) = 
    for jj = 1:can_pw-1
        dist_array(ii) = dist_array(ii) + dis_mat( aux_idx(ii,jj), aux_idx(ii,jj+1), dis_mat2_idx(aux_idx(jj), aux_idx(jj+1)) )
    end
end

[~, s_idx] = sort(dis_mat);

aligned_pos = nan(max_size,can_val1);
for ii = 1:max_size
%     [] = ind2sub([max_size,max_size,can_val1], s_idx(ii));
%     aligned_pos(ii) = ;
    
end

sub2ind()

aux_dist = sum(aux_dist,3);

[~, aux_idx] = sort(aux_dist(:));

ind2sub([cant_anns, length(val1)], aux_idx(1:cant_anns) );

nout = max(nargout,1) - 1;
s = size(x);
for k = 1:nout
   varargout{k} = s(k);
end
