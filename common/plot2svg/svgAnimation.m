function svgAnimation(s, type, key, value, duration)
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Animation')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Animation(next).SubAnimation.Type = type;
    userdata.svg.Animation(next).SubAnimation.Key = key;
    userdata.svg.Animation(next).SubAnimation.Value = value;
    userdata.svg.Animation(next).SubAnimation.Duration = duration;
    set(s(i),'UserData', userdata);
end
