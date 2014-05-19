% script which identifies P wave, as well as its onset and offset
% subtask finding extreme points and deciding if there is wave is the scale n
% Rute Almeida 23/07/02
% Last update: Rute Almeida  24/03/09
%
% Designed for MATLAB Version R12; tested with MATLAB Version R13

% Find biggest positive maximum and negative minimum
maxpos = begwin + modmax(w(begwin+1:endwin,n),2,0,+1);
minpos = begwin + modmax(w(begwin+1:endwin,n),2,0,-1);
[maxim ind] = max(w(maxpos,n));
maxpos = maxpos(ind);
[minim ind] = min(w(minpos,n));
minpos = minpos(ind);

% If no local maxima and minima, take first or last sample
if isempty(maxpos), 
    if (w(begwin,n)>=w(endwin,n)) && w(begwin,n)>0, 
        maxpos = begwin; maxim = w(maxpos,n);
    elseif (w(endwin,n)>=w(begwin,n)) && w(endwin,n)>0,
        maxpos = endwin; maxim = w(maxpos,n);
    end
end
if isempty(minpos),
    if (w(begwin,n)<=w(endwin,n)) && w(begwin,n)<0, 
        minpos = begwin; minim = w(minpos,n);
    elseif (w(endwin,n)<=w(begwin,n)) && w(endwin,n)<0,
        minpos = endwin; minim = w(minpos,n);
    end
end

absmax = abs(maxim);
absmin = abs(minim);

% Calculate Veficaz (Vrms) for the beat at scale n

if  intervalo(1)~=intervalo(2) || length(intervalo)~=2 || (length(time)==1 && intervalo(i)~=1), %#ok<IJCL>  %SL
    if i~=1 %#ok<IJCL> % First beat
        veficaz = sqrt(mean(w(time(i-1):time(i),n).^2)); %#ok<IJCL>
    else       % Other ones
        veficaz = sqrt(mean(w(1:time(i),n)).^2); %#ok<IJCL>
    end
else
    if intervalo(i)~=1 %#ok<IJCL> % First beat
        veficaz = sqrt(mean(w(time(1):time(2),n).^2));
    else       % Other ones
        veficaz = sqrt(mean(w(1:time(1),n)).^2);
    end
end

hay_onda = (abs(maxpos-minpos)<0.11*messages.setup.wavedet.freq)& ((absmax>umbraldetP*veficaz)&(absmin>umbraldetP*veficaz));
if isempty(hay_onda)
    hay_onda = 0;
end
% If maxima too small or too separated, there is no P wave

