%PRTOOLS_MOD Read PRTools Message Of the Day

function mod = prtools_mod

url = 'http://37steps.com/files/prtools_mod.text';

if ~usejava('jvm') & isunix 
	[stat,s] = unix(['wget -q -O - ' url]);
	status = (stat == 0);
else
  if verLessThan('matlab','8.0')
    status = 0; % Cannot read if server is down
  else
    [s,status] = urlread(url,'TimeOut',5);
  end
  %[s,status] = urlread(url);
  % needs updatating. s may contain a firewall message while status = 1
end

if ~status
  mod = [];
else
  mod = s;
end