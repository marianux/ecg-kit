%PREX_DATAFILE  PRTools example of the datafile usage

echo on

prdatafiles
            % makes sure that datafiles are available
						% will download them from the PRTools website when needed

a = highway 
            % load datafile
            % 100 observations of 5 images and a label image
						% R,G,B 
						% and two pixel features (comparisons with previous frames)
						
b = gendat(a,0.06)
            % random selection of 6 observations
						
figure; show(b); drawnow

c = selectim(b,[1 2 3])
            % select RGB
						
figure; show(c); drawnow
showfigs

x = b(2,:)
           % select one observation
figure; show(x,3); drawnow
showfigs

y = data2im(x);
           % select features by retrieving images

data = squeeze(y(:,:,1:5)); % data
lab  = squeeze(y(:,:,6));   % labels
datasize = size(data);

z = im2feat(data);  % store images as feature, pixels become objects
z = setlabels(z,lab(:));
z = setobjsize(z,datasize(1:2));
            % stores pixels as objects with 5 features

trainset = gendat(z,[1000,1000]);
w = qdc(trainset) % train classifier on trainset only
d = z*w*classc;   % classify entire image
figure; show(z*w*classc);
figure; imagesc(d*classim);
figure; imagesc(lab);
showfigs

echo off
