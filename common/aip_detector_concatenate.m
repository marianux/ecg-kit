function payload = aip_detector_concatenate(plA, plB)

    if( ~isempty(plA) )
        if isfield(plA,'result')
            plA = plA.result;
        end
    end

    if( ~isempty(plB) )
        if isfield(plB,'result')
            plB = plB.result;
        end
    end

    payload = ConcatenateQRSdetectionPayloads([], plA, plB);
