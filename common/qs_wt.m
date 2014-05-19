% 
% Prototype:
% ----------
% wt = qs_wt(x, scales, fs, q_filters)
% 
% Description: 
% ------------
% Calculates the Wavelet transform on signal "x", for the scales defined in
% "scales" at sampling frequency "fs". As a result the "wt" signal is
% returned.
% Optionally the "q_filters" object can be defined externally to perform
% the design of the filters only one time.
% 
% Arguments :
% ----------
%       x: k Signal/s of length n. A matrix of n by k
%       scales: A vector indicating which scales to calculate. Default 1:4
%       fs: Sampling frequency. Default 250 Hz
%       q_filters: Optional filter objects designed with
%                  "qs_filter_design".
% 
% Examples :
% ----------
% 
% % Wavelet transform filter bank design
% scales = 3:5;
% sampling_rate = 512;
% CantScales = length(scales);
% MaxScales = max(scales);
% scale_idx = nan(MaxScales,1);
% 
% % Design the filters externally only once.
% filters_cache_filename = ['wt_filters_' num2str(MaxScales) ' scales_' num2str(sampling_rate) ' Hz.mat' ];
% if( exist(filters_cache_filename, 'file') )
%     load( filters_cache_filename );
% else
%     q_filters = qs_filter_design(MaxScales, sampling_rate);
% end
% 
% % Generate a mapping variable to identify the wavelet scales in case not
% % all scales are needed
% for ii = 1:CantScales
%     scale_idx(scales(ii)) = ii;
% end
% 
% wtECG = qs_wt(ECG, scales, sampling_rate, q_filters);
% 
% % All signals at scale 4 are accesed like this.
% wtECG(:,:,scale_idx(4))
% 
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar || unizar.es})
% Birthdate: 17/2/11
% Last update: 20/2/13


function wt = qs_wt(x, scales, fs, q_filters)

    if( nargin < 3 || isempty(fs) )
        fs = 250; 
        warning('Assuming fs=250Hz.');
    end

    if( nargin < 2 || isempty(scales) )
        scales = 1:4; 
        warning('Assuming 4 scales of wavelet decomposition.');
    end

    cant_scales = length(scales);
    max_scale = max(scales);
    
    if( nargin < 4 || isempty(q_filters) )
        filters_cache_filename = ['wt_filters_' num2str(scales) ' scales_' num2str(fs) ' Hz.mat' ];
        if( exist(filters_cache_filename, 'file') )
            load( filters_cache_filename );
        else
            q_filters = qs_filter_design(scales, fs);
            if( usejava('desktop') )
                disp('Please check the filters high-pass frequency response. Press Run (F5) to save these filters for future use.')
                fvtool(q_filters, 'Fs', fs)
                keyboard
            end
            this_path = mfilename('fullpath');
            this_path = fileparts(this_path);
            save([this_path filesep filters_cache_filename], 'q_filters');
        end
    end


    [CantSamp, CantSig] = size(x);
    
    wt = zeros(CantSamp, CantSig, cant_scales);
%     scale_idx = nan(max_scale,1);
    
    for ii = 1:cant_scales
        %indice para saber que escala corresponde a cada columna de MyWavelets
%         scale_idx(scales(ii)) = ii;
        
        aux_delay = grpdelay(q_filters(ii),1);
        cant_casc = length(aux_delay);
        aux_delay = round( sum(aux_delay) );
        aux_filter_length = 2*aux_delay+cant_casc;
        aux = filter(q_filters(ii), x);
        wt(:,:,ii) = [ zeros(aux_delay,CantSig); aux(aux_filter_length:end, :); zeros(aux_delay+1,CantSig) ]; 
    end
