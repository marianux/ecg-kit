%PREX_LOGDENS PRTools example on density based classifier improvement
%
% This example shows the use and results of LOGDENS for improving
% the classification in the tail of the distributions
%
% Note that the use of CLASSC now includes the use of logaritmic densities
% in the tails of distributions. There is no need anymore to call LOGDENS
% explicitely.

	help prex_logdens
	delfigs
	figure
	echo on
				% Generate a small two-class problem
	randreset(7);
	a = gendatb([20 20]);
				% Compute two classifiers: Mixture of Gaussians and Parzen
	w_mogc = mogc(a)*classc;    w_mogc = setname(w_mogc,'MoG');
	w_parz = parzenc(a)*classc; w_parz = setname(w_parz,'Parzen');
				% Scatterplot with MoG classifier
	subplot(3,2,1);
	scatterd(a);
	plotc(w_mogc); xlabel(''); ylabel(''); 
	set(gca,'xtick',[],'ytick',[])
	title('MoG density classifier','fontsize',12)
	drawnow
				% Scatterplot with Parzen classifier
	subplot(3,2,2);
	scatterd(a);
	plotc(w_parz); xlabel(''); ylabel(''); 
	set(gca,'xtick',[],'ytick',[])
	title('Parzen density classifier','fontsize',12)
	drawnow
				% Scatterplot from a distance : 
				% far away points are inaccurately classified
	subplot(3,2,3);
	%scatterd([a; [150 100]; [-150 -100]]);
	scatterd([a; [75 50]; [-75 -50]]);
	plotc(w_mogc); xlabel(''); ylabel(''); 
	set(gca,'xtick',[],'ytick',[])
	title('MoG: bad for remote points','fontsize',12)
	drawnow
				% Scatterplot from a distance : 
				% far away points are inaccurately classified
	subplot(3,2,4);
	scatterd([a; [20 12]; [-20 -12]]); 
	plotc(w_parz); xlabel(''); ylabel(''); 
	set(gca,'xtick',[],'ytick',[])
	title('Parzen: bad for remote points','fontsize',12)
	drawnow
				% Improvement of MOGC by LOGDENS
	subplot(3,2,5);
	%scatterd([a; [150 100]; [-150 -100]]);
	scatterd([a; [75 50]; [-75 -50]]);
	plotc({w_mogc,logdens(w_mogc)},['k--';'r- ']); legend off
	xlabel(''); ylabel(''); set(gca,'xtick',[],'ytick',[])
	title('MoG improved by Log-densities','fontsize',12)
	drawnow
				% Improvement of PARZEN by LOGDENS
	subplot(3,2,6);
	scatterd([a; [20 12]; [-20 -12]]);
	plotc({w_parz,logdens(w_parz)},['k--';'r- ']); legend off
	xlabel(''); ylabel(''); set(gca,'xtick',[],'ytick',[])
	title('Parzen improved by Log-densities','fontsize',12)
	
	echo off
	disp(' ')
	disp('    This example shows the use of the logdens() routine. It')
	disp('    improves the classification in the tails of the distribution,')
	disp('    which is especially important in high-dimensional spaces.')
	disp('    To this end it is combined with normalisation, generating')
	disp('    posterior probabilities. Logdens() can only be applied to')
	disp('    classifiers based on normal densities and Parzen estimates.')
	disp(' ')
	showfigs
