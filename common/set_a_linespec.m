function set_a_linespec(line_hdl, all_proerties)

if(nargin < 1 || any(isempty(line_hdl)) || any(~ishandle(line_hdl)) )
    return;
end

mrk = all_proerties{1};
linestyle = all_proerties{2};
markerSize = all_proerties{3};

colspec_mf = all_proerties{4};
colspec_me = all_proerties{5};
colspec_line = all_proerties{6};

markers_idx = all_proerties{7};
linestyles_idx = all_proerties{8};
markerSize_idx = all_proerties{9};

set(line_hdl,   {'Marker'}, mrk(markers_idx), ... 
                {'MarkerSize'}, markerSize(markerSize_idx), ...
                {'MarkerFaceColor'}, num2cell(colspec_mf,2), ...
                {'MarkerEdgeColor'}, num2cell(colspec_me,2), ...
                {'LineStyle'},linestyle(linestyles_idx) , ...
                {'Color'}, num2cell(colspec_line,2));

