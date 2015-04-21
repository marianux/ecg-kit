function [rotation_matrix, rotated_signal, autovals, reconstruction_matrix] = PiCA(x, peaks, peaks2, win_size, sig2, cant_delays)
%
% [y,W,A] = PiCA(x,peaks1,peaks2)
% Pseudo-Periodic Component Analysis
%
% input:
% x: input data array (channels x samples)
% peaks1: indexes of first signal peaks (R-wave locations for ECG signals)
% peaks2 (optional): indexes of second signal peaks (R-wave locations for ECG signals)
%
% outputs:
% y: pseudo-periodic components ranked in order of resemblance with the
% first to second (if available) desired signals.
% W: the decomposition matrix
% A: the mixing matrix
%
% comment: if the signal peaks2 is not available the components are ranked
% from the most similar to the least similar to the signal with
% corresponding to peaks1.
% 
% Reference:
% R. Sameni, C. Jutten, and M. B. Shamsollahi. Multichannel electrocardiogram
%  decomposition using periodic component analysis. IEEE Transactions on Biomedical Engineering,
%  55(8):1935-1940, Aug. 2008.
%
% Open Source ECG Toolbox, version 2.0, October 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com
% 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% See also autovec_calculation
% 
% Adapted by Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar 
% to ecg-kit toolbox for Matlab.
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
if( nargin < 6 || isempty(cant_delays) )
    cant_delays = 1;
end

if( nargin < 4 )
    win_size = 100;
end

if( nargin < 3 )
    peaks2 = [];
end

[nsamp, nsig] = size(x);

if( nsamp < nsig )
    x = x';
    [nsamp, nsig] = size(x);
end

if( cant_delays < 2)
    
    % traditional PiCA
    
    A = DelayedCovMat(x, peaks, win_size, nsamp, 1);
    B = DelayedCovMat(x, peaks, win_size, nsamp, 0);

    if( isempty(peaks2) )
        [V,D] = eig(A,B,'chol');
    else
        % PM time calculation for the second peaks sequence
        if( nargin < 5 || isempty(sig2) )
            sig2 = x;
        end

    %     figure(99)

        AA = DelayedCovMat(sig2, peaks, win_size, nsamp, 1);

        [V,D] = eig(A-AA,B,'chol');
    end

    d = diag(D);
    
else
    
    % joint diagonalization of several cov matrices for cant_delays delays
    
%     CovMats = nan(4,4,cant_delays+1);
%     for ii = 0:cant_delays
%         CovMats(:,:,ii+1) = DelayedCovMat(x, peaks, win_size, nsamp, ii);
%     end    
%     
%     [V,D]=acdc(CovMats,[],[],[],[],max(max(max(CovMats)))*1e-3);
%     d = mean(D,2);
    
    CovMats = [];
    for ii = 0:cant_delays
        CovMats = [CovMats DelayedCovMat(x, peaks, win_size, nsamp, ii)];
    end    
    
    [ V ,  D ] = joint_diag_r(CovMats,1e-3);    
    
    d = diag(mean(reshape(D,4,4,cant_delays+1), 3));
    
    
end

[autovals, aux_idx] = sort(d, 'descend');
rotation_matrix = V(:,aux_idx);

if nargin > 1
    rotated_signal = x*rotation_matrix;
end

if nargin > 3
    reconstruction_matrix = pinv(rotation_matrix);
end


% function [T0, T1] = CalcSyncIndexes( peaks, win_size, nsamp )
% 
%     peaks = peaks( (peaks - win_size) > 0  & (peaks + win_size ) <= nsamp );
%     aux_seq_idx = 1:length(peaks);
%     T0 = cell2mat(arrayfun(@(a)( peaks(a) - win_size:peaks(a) + win_size ), aux_seq_idx(1:end-1), 'UniformOutput', false));
%     T1 = cell2mat(arrayfun(@(a)( peaks(a) - win_size:peaks(a) + win_size ), aux_seq_idx(2:end), 'UniformOutput', false));
% 
% %     cla
% %     hold on
% %     aux_idx = arrayfun(@(a)( peaks(a) - win_size:peaks(a) + win_size ), aux_seq_idx(1:end-1), 'UniformOutput', false);
% %     cellfun( @(a)( plot(x(a,:)) ), aux_idx )
% %     hold off



    

