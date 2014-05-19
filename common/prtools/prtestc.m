%PRTESTC Test routine for the PRTOOLS classifier
%
% This script tests a given, untrained classifier w, defined in the
% workspace, e.g. w = my_classifier. The goal is to find out whether 
% w fulfills all the requirements of a PRTools classifier. 
% 

if exist('w') ~= 1 | ~ismapping(w) | ~isuntrained(w)
	error('No untrained classifier w found')
end

m = 50;
a = gendath([m,m]);
a = setprior(a,0);
[b,c] = gendat(a,0.5);
name = getname(w);
if (isempty(name))
	error('No name found for the untrained classifier.')
end

disp('')
disp(['     Testing Classifier ' name])
disp(['     ------------------ '])
v = b*w;
newfig(1,3);
scatterd(b);
plotc(v);
if (isempty(getname(v)))
	disp('No name found for the classifier.')
end
disp('Classification error for the Higleyman data: ')
disp(c*v*testc)

disp('Direct output: ')
disp(c(1:5,:)*v)

disp('Classifier output: ')
disp(c(1:5,:)*v*classc)

disp('Soft labels, classification error: ')
a = gendath([m,m],'soft');
a = setprior(a,0);
[b,c] = gendat(a,0.5);
disp(c*(b*w)*testc);

newfig(2,3);
learnsizes = [3,5,7,10,15];
e = cleval(w,b,learnsizes,2,c);
plote(e)
a = gendatm(repmat(20,1,8));
a = setprior(a,0);
[b,c] = gendat(a,0.5);
newfig(3,3)
scatterd(b)
u = b*w;
plotc(u)
disp('Classification error for a multi-class problem: ')
disp(c*u*testc)

if (isaffine(v))
	load nist16_38
	a = prdataset(a);
	v = a*w;
	newfig(6,3);
	show(v);
end

