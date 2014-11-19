function bAux = isOctave()

bAux = false;

octave_ver = ver('Octave');

if( ~isempty(octave_ver) && strcmpi(octave_ver.Name, 'Octave') )
    bAux = true;
end
