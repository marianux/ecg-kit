%% Convert adimentional sample values to target voltage units
% Convert adimentional sample values to target voltage units. The gain and
% zero-offset is included in the header structure.
% 
% Example
% 
%   [x, new_header] = ADC2units(x, header, target_units)
% 
%   where:
%     *x is a matrix with signals in the columns
%     *header is a structure (ECG_header prop. in ECGwrapper object)
%        describing the signal. Mandatory fields for this function are:
%        units, gain and zero. 
%     *target_units is a string to convert the ADC sample values to. See
%        cTypicalUnits below.
% 
% See also ADC2realunits, ECGwrapper
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Birthdate: 01/01/2012
% Last update: 18/10/2014
% Copyright 2008-2015% Version: 0.1 beta
% Birthdate: 01/01/2012
% Last update: 18/10/2014
% Copyright 2008-2015
function [x, new_header] = ADC2units(x, header, target_units)

if( nargin < 3 || isempty(target_units) )
   target_units = 'MICROVOLTIOS';
end

cTypicalUnits = { ... 
                'NV', 'NANOVOLTS' , 'NANOVOLTIOS' ;...
                'UV', 'MICROVOLTS' , 'MICROVOLTIOS'; ...
                'MV', 'MILIVOLTS' , 'MILIVOLTIOS'; ...
                'V', 'VOLTS' , 'VOLTIOS'; ...
                };

UnitFactor = colvec(10.^(-9:3:0));
            

origin_units_idx = [];
target_units_idx = [];
for ii = 1:size(cTypicalUnits,1)

    if( isempty(target_units_idx) )
        target_units_idx = find(strcmpi(cTypicalUnits(ii,:), upper(deblank(target_units))) );
        if(~isempty(target_units_idx))
            target_units_idx = ii;
        end
    end
    
    if( isempty(origin_units_idx) )
        origin_units_idx = find(strcmpi(cTypicalUnits(ii,:), upper(deblank(header.units(1,:)))) );
        if(~isempty(origin_units_idx))
            origin_units_idx = ii;
        end
    end

    if( ~isempty(origin_units_idx) && ~isempty(target_units_idx) )
        break
    end
end

if( isempty(origin_units_idx) )
%     error('ADC2units:UnknownUnits', ['Unknown units: ' header.units(1,:) ]);
%     fprintf(2, 'Unknown units: %s\n', header.units(1,:) );
    new_header = header;
    return
end

if( isempty(target_units_idx) )
%     error('ADC2units:UnknownUnits', ['Unknown units: ' target_units ]);
%     fprintf(2, ['Unknown units: %s\n', target_units]);
    new_header = header;
    return
end

x = ADC2realunits(x, header.adczero, header.gain);

x = x * UnitFactor(origin_units_idx) / UnitFactor(target_units_idx);
cant_sig = size(x,2);

new_header = header;
new_header.gain = ones(cant_sig,1);
new_header.adczero = zeros(cant_sig,1);
new_header.units = repmat(target_units, cant_sig,1);
