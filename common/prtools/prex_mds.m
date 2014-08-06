%PREX_MDS   PRTools example on multi-dimensional scaling
%
% Show the training, generalisation of some non-linear mappings for
% visualisation.
%
help prex_mds
delfigs

echo on

  a = satellite;         % 36D dataset, 6 classes, 6435 objects
  [x,y] = gendat(a,0.1); % split in train and test set
  
  % TSNEM
  wt = x*tsnem;          % train TSNEM
  figure; 
  scattern(x*wt);        % show 2d result for trainset
  title('tSNEM trainset')
  figure; 
  scattern(y*wt);        % show 2d result for testset
  title('tSNEM testset')
  showfigs
  
  % SAMMONM
  ws = x*sammonm;         % train SAMMONM     
  figure; 
  scattern(x*ws);         % show 2d result for trainset
  title('Sammon trainset')
  figure; 
  scattern(y*ws);         % show 2d result for testset
  title('Sammon testset')
  showfigs
  
  % MDS
  dxx = sqrt(distm(x,x)); % dissimilarity matrix of trainset
  wm = mds(dxx);          % train MDS           
  figure; 
  scattern(dxx*wm);       % show 2d result for trainset
  title('MDS trainset')
  figure; 
  dyx = sqrt(distm(y,x)); % dissimilarity between testset and trainset
  scattern(dyx*wm);       % show 2d result for testset
  title('MDS testset')
  showfigs
  
echo off
