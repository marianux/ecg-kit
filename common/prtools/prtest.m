%PRTEST test routine for PRTools commands
%
% This functions runs countless basic Prtools operations on basic
% Prtools objects. Set variable 'item' for specific choices.
% Possible 'item' values are: dataset, labtype, mapping, fixed, classifiers

if ~exist('item'), item = 'prdataset'; end

echo on

switch item
case 'prdataset'

	disp('dataset definition')
	r = rand(10,2);
	lab = genlab([3,7])
	
	a = prdataset(r,lab)
	lab = genlab([3,7],['apple1';'apple2'])
	a = prdataset(r,lab)
	lab = genlab([3,7],{'pear';'apple'})
	a = prdataset(r,lab)
	struct(a)
	
	getdata(a)
	a.data
	get(a,'data')
	
	getlabels(a)
	a.labels
	get(a,'labels')
	
	getnlab(a)
	a.nlab
	get(a,'nlab')
	
	getfeatlab(a)
	a.featlab
	get(a,'featlab')
	
	a = set(a,'featlab',{'feat1';'f2'})
	get(a,'featlab')
	
	getprior(a)
	a.prior
	get(a,'prior')
	
	getlablist(a)
	a.lablist
	get(a,'lablist')
	
	getlabtype(a)
	a.labtype
	get(a,'labtype')
	
	size(a)
	getobjsize(a)
	a.objsize
	get(a,'objsize')
	
	getfeatsize(a)
	a.featsize
	get(a,'featsize')
	
	getident(a)
	a.ident
	get(a,'ident')
	
	getversion(a)
	a.version
	a.version{1}
	a.version{2}
	get(a,'version')
	
	a = setname(a,'test of dataset definition')
	getname(a)
	a.name
	get(a,'name')
	
	a = setuser(a,'some rubbish')
	getuser(a)
	a.user
	get(a,'user')

case 'labtype'
	a = gendath(5)
	a = setlabels(a,genlab([3 2 3 2],[4 3 1 2]'))
	getlablist(a)'
	getlabels(a)'
	gettargets(a)'
	
	b = setlabtype(a,'soft')
	getlablist(b)'
	getlabels(b)'
	gettargets(b)'
	
	c = setlabtype(a,'targets')
	getlablist(c)'
	getlabels(c)'
	gettargets(c)'
	
	a = setlabtype(b,'crisp')
	getlablist(a)'
	getlabels(a)'
	gettargets(a)
	
	a = setlabtype(b,'targets')
	getlablist(a)'
	getlabels(a)'
	gettargets(a)
	
	a = setlabtype(c,'crisp')
	getlablist(a)'
	getlabels(a)'
	gettargets(a)
	
	a = setlabtype(c,'soft')
	getlablist(a)'
	getlabels(a)'
	gettargets(a)
	
	

case 'prmapping'

	a = prdataset(rand(10,2),genlab([5,5]))
	w = affine(rand(2,2),[1 0])
	
	struct(w)
	
	getmapping_file(w)
	w.mapping_file
	get(w,'mapping_file')
	
	getmapping_type(w)
	w.mapping_type
	get(w,'mapping_type')
	
	getdata(w)
	w.data
	get(w,'data')
	
	getlabels(w)
	w.labels
	get(w,'labels')
	
	getsize_in(w)
	w.size_in
	get(w,'size_in')
	
	getsize_out(w)
	w.size_out
	get(w,'size_out')
	
	getsize(w)
	w.size
	get(w,'size')
	
	getscale(w)
	w.scale
	get(w,'scale')
	
	getout_conv(w)
	w.out_conv
	get(w,'out_conv')
	
	getname(w)
	w.name
	get(w,'name')
	
	getuser(w)
	w.user
	get(w,'user')
	
	v = [w w w]
	
	v = [w; w; w]
	
case 'fixed'

	a = prdataset(rand(5,2),[1 1 2 2 2]');
	disp('Testing fixed mappings')
	disp('Dataset')
	+a
	disp('feature selection')
	b = [a a a];
	+b
	w = cmapm(size(b,2),[1 3 4 6])
	+(b*w)
	disp('sigmoid')
	w = sigm([],2)
	+(a*w)
	disp('exp')
	w = cmapm(2,'exp')
	+(a*w)
	disp('nexp')
	w = cmapm(2,'nexp')
	+(a*w)
	disp('log')
	w = cmapm(2,'log')
	+(a*w)
	disp('randrot')
	w = cmapm(2,'randrot')
	+(a*w)
	disp('rot')
	w = cmapm([1 1;1 -1],'rot')
	figure
	scatterd(a);
	figure
	scatterd(a*w)
	disp('shift')
	w = cmapm([100 50],'shift');
	disp([mean(a);mean(a*w)])
	w = cmapm([100,0.01],'scale');
	disp([std(a);std(a*w)])

case 'classifiers'

	W = {fisherc, quadrc, qdc, ldc, udc, nmc, nmsc, knnc, knnc([],1),parzenc, ...
	     klldc, pcldc, svc, loglc, treec, neurc, lmnc, bpxnc}
	echo off
	for j=1:length(W)
		w = W{j};
		prtestc
	end

case 'mappings'
	W = {cmapm, scalem, klm, pca, klms, reducm, fisherm, nlfisherm, kernelm, ...
	      sigm, invsigm, distm, proxm, lmnm, normm, svm}

case 'combining'
	echo off
	W = {maxc, minc, meanc, prodc, medianc, votec}
	a = gendath;
	[b,c] = gendat(a,0.2);
	w = [fisherc, nmc, nmsc, qdc, knnc([],3)];
	for j=1:length(W)
		v1 = b*w*W{j};
		v2 = b*(w*W{j});
		v3 = b*w*classc*W{j};
		v4 = b*(w*classc*W{j});
		e1 = c*{v1 v2 v3 v4}*testc;
		e2 = c*(b*w*W{j})*testc;
		e3 = c*(b*w*classc*W{j})*testc;
		e4 = c*(b*w*W{j}*classc)*testc;
		disp(getname(W{j}))
		disp([e1{:} e2 e3 e4])
	end
		  
end

echo off

