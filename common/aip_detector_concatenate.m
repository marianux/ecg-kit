function payload = aip_detector_concatenate(plA, plB)

    if( ~isempty(plA) )
        plA = plA.result;
    end

    if( ~isempty(plB) )
        plB = plB.result;
    end

    payload = ConcatenateQRSdetectionPayloads([], plA, plB);
