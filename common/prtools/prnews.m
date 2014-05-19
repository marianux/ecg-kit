%PRNEWS Opens browser with PRTools news

url = 'http://www.37steps.com/prtools/updates/';
[s,status] = urlread(url);

if ~status
  error('Website could not be reached. Please retry later')
else
  web(url,'-browser')
end