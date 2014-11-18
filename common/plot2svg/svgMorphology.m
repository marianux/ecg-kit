function svgMorphology(s, source, operator, radius, result)
% Adds a feMorphology SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
%
% svgMorphology(s, source, operator, radius, result)
% Parameters:
%   s : Array of plot object handles
%   source : Any previous defined filter result string, 'SourceGraphic',
%            or 'SourceAlpha'.
%   operator : Type of morphology [erode,dilate]
%   radius : Radius
%   result : String that identifies the filter result for following filter
%            stages.   
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feMorphology';
    userdata.svg.Filter(next).Subfilter.Source = source;    
    userdata.svg.Filter(next).Subfilter.Operator = operator;
    userdata.svg.Filter(next).Subfilter.Radius = radius;
    userdata.svg.Filter(next).Subfilter.Result = result;
    set(s(i),'UserData', userdata);
end
