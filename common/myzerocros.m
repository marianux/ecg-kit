function index = myzerocros(x)
% Returns the index of the input vector in which the first zero crossing is located.

sx = sign(x);
diffsx = diff(sx);
diffsx = [diffsx(1); diffsx];
index = find( diffsx ~= 0);

