%PRTOOLS_MOD Read PRTools Message Of the Day

function mod = prtools_mod

url = 'http://37steps.com/files/prtools_mod.text';

if ~usejava('jvm') & isunix 
	[stat,s] = unix(['wget -q -O - ' url]);
	status = (stat == 0);
else
  [s,status] = urlread(url);
end

if ~status
  mod = [];
else
  mod = s;
end