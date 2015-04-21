%% (Internal) Add the path only if not was already added.
%   
%   paths_added = addpath_if_not_added(str_path, position)
% 
% Arguments:
% 
%      + str_path: paths to add
% 
%      + position: where to add. See addpath doc.
% 
% Output:
% 
%      + paths_added: Only the paths added
% 
% Example:
% 
% See also CalcRRserieQuality
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function [corr_gain noise_power] = calc_correlation_gain(signal, heasig, references, limits, noise_power, bRobust)

corr_gain = nan;
[nsamp, nsig] = size(signal);
pattern_size = sum(limits);
lreferences = length(references);

sig_pack = pack_signal(signal, references(references > limits(1) & references < (nsamp-limits(2) )), limits, true);    

if( isempty(sig_pack) )
    noise_power = nan;
    return
end

if( nargin < 5 || isempty(bRobust) )
    pattern_avg = squeeze(median(sig_pack,3));
else
    pattern_avg = squeeze(mean(sig_pack,3));
end

sig_pack_sub_avg = bsxfun(@minus, sig_pack, pattern_avg);

corr_gain = 10*log10( mean(pattern_avg.^2) ./ 10.^squeeze(mean(log10(mean(sig_pack_sub_avg.^2,1)),3)) );
 
% % Obsolete: estimate noise from a random sampling simulating guessing the QRS locations. 
% 
% if( nargin < 4 || isempty(noise_power) )
%     noise_power = nan(30,nsig);
%     noise_disp = nan(30,nsig);
%     for ii = 1:30
%         noise_pack = pack_signal(signal, sort(randsample(limits(1):(nsamp-limits(2)), lreferences)), limits);    
% %         noise_avg = flipud(reshape(squeeze(median(noise_pack,3)), pattern_size, nsig));
%         noise_limits = prctile(noise_pack, [5 50 95], 3);
%         
%         noise_down = squeeze(noise_limits(:,:,1));
%         noise_avg = squeeze(noise_limits(:,:,2));
%         noise_up = squeeze(noise_limits(:,:,3));
%         
%         noise_power(ii,:) = mean(noise_avg.^2);
%         noise_disp(ii,:) = mean(noise_up - noise_down);
%     end
%     noise_power = nanmedian(noise_power);
%     noise_disp = nanmedian(noise_disp);
% end
% 
% corr_gain = 10*log10(mean(pattern_avg.^2) ./ noise_power );


% pattern_limits = prctile(sig_pack, [5 50 95], 3);
% 
% pattern_down = squeeze(pattern_limits(:,:,1));
% pattern_up = squeeze(pattern_limits(:,:,3));
% mean_pattern_disp = mean(pattern_up - pattern_down);

% corr_gain_2 = 20*log10( noise_disp ./ mean_pattern_disp );


% debug
% aux_val = [];
% for ii  = 1:100
%     noise_pack = pack_signal(signal, sort(randsample(limits(1):(nsamp-limits(2)), length(references))), limits);    
%     noise_avg = flipud(reshape(squeeze(median(noise_pack)), pattern_size, nsig));
%     aux_val = [aux_val; mean(noise_avg.^2)];
% end
% 
% figure
% hist( aux_val, 100)
% legend()
