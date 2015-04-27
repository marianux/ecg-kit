%% (Internal) Design the wavelet decomposition filters for wavedet algorithm
%
% Prototype:
% ----------
% q_filters = qs_filter_design(scales, fs, N);
% 
% Description: 
% ------------
% Mimics the transfer function of the filters used for ECG delineation in
% [Martinez et al. 2004] for an arbitrary sampling frequency and filter order N.
% 
% Martinez et al. "A Wavelet-Based ECG Delineator: Evaluation on Standard
% Databases" IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. 51, NO. 4,
% APRIL 2004.
% 
% WARNING :
% ---------
% As this routines iterates through several configurations in order to
% converge, the user should check the transfer functions of the filters. My
% suggestion is once you obtain a desired filter bank for a given Fs or
% configuration, save or cache it in a .mat file in order to use it during
% operation. For example, if you usually work with signals sampled at 360
% Hz, a good choice is to have a cached version of the filters for this Fs
% in a .mat file called "wt_filters_6 scales_360 Hz.mat". You can use this
% function on-line with your algorithm at your own risk.
% 
% Examples :
% ----------
% 
% q_filters = qs_filter_design(4, 250);
% 
% q_filters = qs_filter_design(5, 360);
% 
% q_filters = qs_filter_design(6, 1000);
% 
% you can check filter characteristics using:
% 
% fss = [250 360 500 1000]; %Hz
% for fs  = fss
%     fvtool(q_filters, 'fs', fs )
% end
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate: 17/2/11
% Last update: 22/02/13
% Copyright 2008-2015
% 
function q_filters = qs_filter_design(scales, fs, N)

design_iter_attemps = 10;

if( nargin < 2 || isempty(fs) )
    fs = 250; %Hz
end

if( nargin < 3 || isempty(N) )
    %valor empírico obtenido de varios diseños.
    N = max(10, round(fs*4/30 + 16+2/3)); 
end

%Pruebo que N no sea demasiado distinto a lo recomendable.
recommended_N = max(10, round(fs*4/30 + 16+2/3));
if( abs(N - recommended_N) > 0.1*N )
    warning(['Check the transfer functions of the differentiator filters designed. Recommended order N = ' num2str(recommended_N) ]);
end

%frecuencia a la que está diseñado el delineador, y que se toma para
%referencia para que las escalas signifiquen lo mismo a cualquier Fs.
f_ref = 250; %Hz
f_ratio = f_ref/fs;


% Funciones de transferencia correctas para el diseño de los filtros utilizadas 
% por el wavedet a 250 Hz.
empirical_tf = { ...
    [2;-2] ... 
    [0.250000000000000;0.750000000000000;0.500000000000000;-0.500000000000000;-0.750000000000000;-0.250000000000000;] ... 
    [0.0312500000000000;0.0937500000000000;0.187500000000000;0.312500000000000;0.343750000000000;0.281250000000000;0.125000000000000;-0.125000000000000;-0.281250000000000;-0.343750000000000;-0.312500000000000;-0.187500000000000;-0.0937500000000000;-0.0312500000000000;] ... 
    [0.00390625000000000;0.0117187500000000;0.0234375000000000;0.0390625000000000;0.0585937500000000;0.0820312500000000;0.109375000000000;0.140625000000000;0.160156250000000;0.167968750000000;0.164062500000000;0.148437500000000;0.121093750000000;0.0820312500000000;0.0312500000000000;-0.0312500000000000;-0.0820312500000000;-0.121093750000000;-0.148437500000000;-0.164062500000000;-0.167968750000000;-0.160156250000000;-0.140625000000000;-0.109375000000000;-0.0820312500000000;-0.0585937500000000;-0.0390625000000000;-0.0234375000000000;-0.0117187500000000;-0.00390625000000000;] ... 
    [0.000488281250000000;0.00146484375000000;0.00292968750000000;0.00488281250000000;0.00732421875000000;0.0102539062500000;0.0136718750000000;0.0175781250000000;0.0219726562500000;0.0268554687500000;0.0322265625000000;0.0380859375000000;0.0444335937500000;0.0512695312500000;0.0585937500000000;0.0664062500000000;0.0727539062500000;0.0776367187500000;0.0810546875000000;0.0830078125000000;0.0834960937500000;0.0825195312500000;0.0800781250000000;0.0761718750000000;0.0708007812500000;0.0639648437500000;0.0556640625000000;0.0458984375000000;0.0346679687500000;0.0219726562500000;0.00781250000000000;-0.00781250000000000;-0.0219726562500000;-0.0346679687500000;-0.0458984375000000;-0.0556640625000000;-0.0639648437500000;-0.0708007812500000;-0.0761718750000000;-0.0800781250000000;-0.0825195312500000;-0.0834960937500000;-0.0830078125000000;-0.0810546875000000;-0.0776367187500000;-0.0727539062500000;-0.0664062500000000;-0.0585937500000000;-0.0512695312500000;-0.0444335937500000;-0.0380859375000000;-0.0322265625000000;-0.0268554687500000;-0.0219726562500000;-0.0175781250000000;-0.0136718750000000;-0.0102539062500000;-0.00732421875000000;-0.00488281250000000;-0.00292968750000000;-0.00146484375000000;-0.000488281250000000;] ... 
    [6.10351562500000e-05;0.000183105468750000;0.000366210937500000;0.000610351562500000;0.000915527343750000;0.00128173828125000;0.00170898437500000;0.00219726562500000;0.00274658203125000;0.00335693359375000;0.00402832031250000;0.00476074218750000;0.00555419921875000;0.00640869140625000;0.00732421875000000;0.00830078125000000;0.00933837890625000;0.0104370117187500;0.0115966796875000;0.0128173828125000;0.0140991210937500;0.0154418945312500;0.0168457031250000;0.0183105468750000;0.0198364257812500;0.0214233398437500;0.0230712890625000;0.0247802734375000;0.0265502929687500;0.0283813476562500;0.0302734375000000;0.0322265625000000;0.0339965820312500;0.0355834960937500;0.0369873046875000;0.0382080078125000;0.0392456054687500;0.0401000976562500;0.0407714843750000;0.0412597656250000;0.0415649414062500;0.0416870117187500;0.0416259765625000;0.0413818359375000;0.0409545898437500;0.0403442382812500;0.0395507812500000;0.0385742187500000;0.0374145507812500;0.0360717773437500;0.0345458984375000;0.0328369140625000;0.0309448242187500;0.0288696289062500;0.0266113281250000;0.0241699218750000;0.0215454101562500;0.0187377929687500;0.0157470703125000;0.0125732421875000;0.00921630859375000;0.00567626953125000;0.00195312500000000;-0.00195312500000000;-0.00567626953125000;-0.00921630859375000;-0.0125732421875000;-0.0157470703125000;-0.0187377929687500;-0.0215454101562500;-0.0241699218750000;-0.0266113281250000;-0.0288696289062500;-0.0309448242187500;-0.0328369140625000;-0.0345458984375000;-0.0360717773437500;-0.0374145507812500;-0.0385742187500000;-0.0395507812500000;-0.0403442382812500;-0.0409545898437500;-0.0413818359375000;-0.0416259765625000;-0.0416870117187500;-0.0415649414062500;-0.0412597656250000;-0.0407714843750000;-0.0401000976562500;-0.0392456054687500;-0.0382080078125000;-0.0369873046875000;-0.0355834960937500;-0.0339965820312500;-0.0322265625000000;-0.0302734375000000;-0.0283813476562500;-0.0265502929687500;-0.0247802734375000;-0.0230712890625000;-0.0214233398437500;-0.0198364257812500;-0.0183105468750000;-0.0168457031250000;-0.0154418945312500;-0.0140991210937500;-0.0128173828125000;-0.0115966796875000;-0.0104370117187500;-0.00933837890625000;-0.00830078125000000;-0.00732421875000000;-0.00640869140625000;-0.00555419921875000;-0.00476074218750000;-0.00402832031250000;-0.00335693359375000;-0.00274658203125000;-0.00219726562500000;-0.00170898437500000;-0.00128173828125000;-0.000915527343750000;-0.000610351562500000;-0.000366210937500000;-0.000183105468750000;-6.10351562500000e-05;] ... 
                };

Grid_size = 1024;

% Creo una grilla de muestreo en frecuencia logaritmica para que la
% zona en que derivan los filtros sea una recta de pendiente constante.
F_log = logspace(-3, min(0, log10(fs/f_ref)), Grid_size) * pi;

filter_count = 1;

for ii = rowvec(scales)

    [g_amp F] = freqz(empirical_tf{ii},1, F_log);
    F = F * f_ratio / pi;
    g_amp = abs(g_amp);

    %averiguo hasta qué muestra se comporta como un derivador, ya que será
    %un parámetro de diseño.
    slope_aux = diff( log10(g_amp) );
    end_diff_idx = find(slope_aux < 0.95*slope_aux(1), 1, 'first');
    

    %diseño un derivador hasta dicha frecuencia, con una banda de
    %transición dada por la expresion , cuando sea posible.
    ftrans = min(70, 251 * exp(-0.63*ii));
    diff_order = round(N*1.8.^((ii*0.4 -0.6))); 
    %Los N tienen que ser necesariamente pares para el diseño del derivador.
    if( rem(diff_order,2) ~= 0 )
        diff_order = diff_order + 1;
    end
    
    [msgstr, msgid] = lastwarn;
    %itero hasta que se diseña correctamente.
    jj = 0;
    effective_order = diff_order;
    while( jj < design_iter_attemps )
        d = fdesign.differentiator('n,fp,fst', effective_order, F(end_diff_idx), min( 0.95, F(end_diff_idx)+(ftrans*2/fs) ) );
        try
            Hd = design(d,'equiripple');  
            [~, msgid] = lastwarn('','');
        catch MException
            %on error, force iteration with different config.
            msgid = 'signal:firpm:DidNotConverge';
        end
        if(strcmpi(msgid,'signal:firpm:DidNotConverge'))
            %fallo en la convergencia, iteramos de nuevo.
            %recorro linealmente por el rango +10:-50% 
            effective_order = round(diff_order * (jj*(0.5-1.1)/design_iter_attemps+1.1));
            %Los effective_order tienen que ser necesariamente pares para el diseño del derivador.
            if( rem(effective_order,2) ~= 0 )
                effective_order = effective_order + 1;
            end            
            jj = jj + 1;
        else
            jj = design_iter_attemps;
        end
    end    
    
    if(strcmpi(msgid,'signal:firpm:DidNotConverge'))
    %fallo en la convergencia, reportamos el error.
        error('qs_filter_design:Impossible2Design', 'Impossible to design the differentiator filter. Please try another N value, or check the filters transfer functions manually.')
    end
    
    %Vuelvo a ver la transferencia real del derivador diseñado y del filtro
    %a emular para que tengan una respuesta similar en la zona en que se
    %comporta como derivador. Averiguo el factor de escala entre ambas
    %transferencias.
    g_amp_emp = freqz(empirical_tf{ii},1, F_log);
    g_amp_emp = abs(g_amp_emp);
    
    g_amp_real = freqz(Hd, F * pi);
    g_amp_real = abs(g_amp_real);
    
    aux_idx = round(end_diff_idx/2);
    aux_scale = g_amp_emp(aux_idx)/g_amp_real(aux_idx);
    
    %Escalo la zona del derivador.
    Hd.Numerator = Hd.Numerator * aux_scale;
    g_amp_real = g_amp_real * aux_scale;
    

    %ahora diseño un filtro de compensación para que la zona en que no se
    %deriva sea similar al filtro de referencia. Esta irá desde el máximo
    %de la transferencia pasobanda hasta donde haya una amplitud mayor a
    %0.5 veces.
    [~, max_idx] = max(g_amp);

    aux_3db_unique = find( g_amp_real > 0.5 & g_amp_emp > 0.5, 1, 'last');
        
    if( max_idx >= aux_3db_unique)
        max_idx = end_diff_idx;
    end
    aux_idx = max_idx:aux_3db_unique; 

    %Creo una grilla de frecuencia - magnitud arbitraria para el diseño del
    %filtro. Esta transferencia no deberá afectar la zona derivadora (H(w)
    %= 1) y compensar la zona indicada por aux_idx, emulando la
    %transferencia de referencia y teniendo un valor de atenuación
    %considerable en Nyquist (H(w) = 1e-3)
    F_comp = [0 max(0.9*F(end_diff_idx), (F(max_idx)-F(end_diff_idx))/2) F(aux_idx) 1];
    g_amp_comp = [1 1 g_amp(aux_idx)./g_amp_real(aux_idx) 1e-3];
    W = ones(1, length(F_comp));
    W(3:end-1) = 10;

    %Se diseña el filtro 
    d = fdesign.arbmag('N,F,A', diff_order, F_comp, g_amp_comp);
    Hd_comp = design(d,'firls', 'weights', W, 'FilterStructure', 'dfsymfir');     

    % y se cascadea al original.
    Hd = dfilt.cascade(Hd, Hd_comp);

    q_filters(filter_count) = Hd;
    
    filter_count = filter_count + 1;
    
end
