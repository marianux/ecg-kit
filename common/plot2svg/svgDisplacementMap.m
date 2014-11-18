function svgDisplacementMap(s, source1, source2, scale, xChannel, yChannel, result)
% Adds a feDisplacementMap SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
%
% svgDisplacementMap(s, source1, source2, scale, xChannel, yChannel, result)
% Parameters:
%   s : Array of plot object handles
%   source1 : Any previous defined filter result string, 'SourceGraphic',
%             or 'SourceAlpha'.
%   source2 : Any previous defined filter result string, 'SourceGraphic',
%             or 'SourceAlpha'.
%   scale : 
%   xChannel:
%   yChannel:
%   result : String that identifies the filter result for following filter
%            stages.   
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feDisplacementMap';
    userdata.svg.Filter(next).Subfilter.Source1 = source1;
    userdata.svg.Filter(next).Subfilter.Source2 = source2;
    userdata.svg.Filter(next).Subfilter.Scale = scale;
    userdata.svg.Filter(next).Subfilter.xChannel = xChannel;
    userdata.svg.Filter(next).Subfilter.yChannel = yChannel;    
    userdata.svg.Filter(next).Subfilter.Result = result;
    set(s(i),'UserData', userdata);
end