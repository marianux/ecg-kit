function [x new_header] = ADC2units(x, header, target_units)
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
    error('ADC2units:UnknownUnits', ['Unknown units: ' header.units(1,:) ]);
end

if( isempty(target_units_idx) )
    error('ADC2units:UnknownUnits', ['Unknown units: ' target_units ]);
end

x = ADC2realunits(x, header.adczero, header.gain);

x = x * UnitFactor(origin_units_idx) / UnitFactor(target_units_idx);
cant_sig = size(x,2);

new_header = header;
new_header.gain = ones(cant_sig,1);
new_header.adczero = zeros(cant_sig,1);
new_header.units = repmat(target_units, cant_sig,1);
