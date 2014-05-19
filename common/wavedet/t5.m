% script which used ascale 5 in T wave delineation 
% Rute Almeida 
% Last update: Rute Almeida  19/07/2011
%
% Designed for MATLAB Version R12; tested with MATLAB Version R13

if absmax >= absmin,    % the greatest modulus maximum is the maximum
    % Now we search the two minima nearest to maxpos, one before and one after
    minapos = max(modmax(w(begwin+1:maxpos-1,scale2),2,0,-1));
    minapos = begwin + minapos; % Position of the negative minimum before the maximum
    if isempty(minapos) && (maxpos ~= begwin) && (w(begwin,scale2)<0),
        minapos = begwin;        % If no local minimum before the maximum, take the first sample
    end
    minppos = min(modmax(w(maxpos+1:endwin,scale2),2,0,-1));
    minppos = maxpos + minppos; % Position of the positive maximum after the minimum
    if isempty(minppos) && (maxpos ~= endwin) && (w(endwin,scale2)<0),
        minppos = endwin;        % If no local minimum after the maximum, take the last sample
    end
    
    mina = abs(w(minapos,scale2));     % Amplitude of minimum before maximum
    minp = abs(w(minppos,scale2));     % Amplitude of minimum after maximum
    
    if (mina < umbralsig*absmax), % If mina is not big enough
        mina =[];                  % forget it
    elseif (maxpos-minapos>Tmax_Tmin_time_min*messages.setup.wavedet.freq), %or if ther are more than 150 ms to maxpos
        mina = [];
    end
    if (minp < umbralsig*absmax), % If minp is not big enough
        minp =[];                  % forget it
    elseif (minppos-maxpos>Tmax_Tmin_time_min*messages.setup.wavedet.freq), % and also if there are more than 150 ms to maxpos
        minp =[];
    end
    
    if ~isnan(mina)&~isnan(minp),  %#ok<AND2>    %%% NUEVO JP
        if (mina >= minp)&&(minp < umbralsig*absmax*Tmax_Tmin_bifasic),
            minp = [];
        elseif (minp> mina)&&(mina < umbralsig*absmax*Tmax_Tmin_bifasic),
            mina = [];
        end
    end
    
    % Test which modulus maxima are significative and find zero crossings
    if isempty(mina),
        if isempty(minp),
            tipoT = 2;      % only upwards T wave
            if maxpos - minapos > 2,         % if not !!!!!!!!!!!
                ind = zerocros(flipud(w(minapos:maxpos,scalezerocros)));   %Scale scalezerocros !!!
                if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                    ind = zerocros(flipud(w(minapos:maxpos,scale2)));
                end %%%%%%%%%%%%%%Rute 03/09/02
                T = maxpos - ind +1;                           % Zero crossing = T wave position
                picoff = maxpos;                               %wavelet  peak to detect offset
            elseif isempty(minapos);  % If there were no minimum, there is no zero crossing				
                T = picant (w(begwin:maxpos,scale2),maxpos);       % Take the minimum at scale scale2
                if ~isempty(T) %%%%%%%%%%%% 14/06 /02 Rute 
                    picoff = maxpos;
                else
                    picoff = []; % if did not exist a peak in scale scale2 there is no T
                end %%%%%%%%%%%% 14/06/02 Rute 
            else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
              messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
            end
        else      % minp exists (is significative) but mina not
            tipoT = 0;  		%normal T wave
            if minppos -maxpos >2,       % if not!!!!!!???
                ind = zerocros(w(maxpos:minppos,scalezerocros));  %% 07/06/02 Rute
                if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                    ind = zerocros(w(maxpos:minppos,scalezerocros2));
                end %%%%%%%%%%%%%%Rute 03/09/02
                T = maxpos + ind -1;
                picon = maxpos;		% For determining onset and offset
                picoff = minppos;
            else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
            end
        end
    else
        if isempty(minp),   %  mina exists (is significative) but minp not
            tipoT = 1;  	%inverted T wave
            if maxpos -minapos >2,  % if not !!!!!!!!!
                ind = zerocros(w(minapos:maxpos,scalezerocros));  % wavelet zero crossing is T wave peak  %% 07/06/02 Rute
                if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                    ind = zerocros(w(minapos:maxpos,scalezerocros2));
                end %%%%%%%%%%%%%%Rute 03/09/02
                T = minapos + ind -1;
                picon = minapos;
                picoff = maxpos;
            else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
               messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
            end
        else    % both mina and minp are significative.  Biphasic wave.
            tipoT = 5;	% biphasic neg-pos T wave
            if maxpos - minapos > 2,    %!!!!!!!!!!!!
                ind = zerocros(flipud(w(minapos:maxpos,scalezerocros))); %% 07/06/02 Rute
                if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                    ind = zerocros(flipud(w(minapos:maxpos,scalezerocros2)));
                end %%%%%%%%%%%%%%Rute 03/09/02
                T = maxpos - ind +1;
                picon = minapos;
            else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
               messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
            end
            if minppos - maxpos > 2,    %!!!!!!!!!!!!
                ind = zerocros(flipud(w(maxpos:minppos,scalezerocros))); %% 07/06/02 Rute
                if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                    ind = zerocros(flipud(w(minapos:maxpos,scalezerocros2))); 
                end %%%%%%%%%%%%%%Rute 03/09/02
                Tprima = maxpos + ind -1;
                picoff = minppos;
            else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
               messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
            end
        end
    end
else        % If the greatest modulus maximum is the minimum
    % Search two maxima, one before and one after the minimum
    maxapos = max(modmax(w(begwin+1:minpos-1,scale2),2,0,1));
    maxapos = begwin + maxapos;
    if isempty(maxapos) && (minpos ~= begwin) && (w(begwin,scale2)>0),
        maxapos = begwin;           
    end
    maxppos = min(modmax(w(minpos+1:endwin,scale2),2,0,1));
    maxppos = minpos + maxppos;
    if isempty(maxppos) && (minpos ~= endwin) && (w(endwin,scale2)>0),
        maxppos = endwin;           
    end
    maxa = abs(w(maxapos,scale2));
    maxp = abs(w(maxppos,scale2));
    % See if they are significative
    if (maxa < umbralsig*absmin)
        maxa =[];
    elseif (minpos-maxapos>Tmax_Tmin_time_min*messages.setup.wavedet.freq),
        maxa = [];
    end
    if (maxp < umbralsig*absmin),
        maxp =[];
    elseif (maxppos-minpos>Tmax_Tmin_time_min*messages.setup.wavedet.freq),
        maxp = [];
    end
    if ~isnan(maxa)&~isnan(maxp),   %#ok<AND2>   %%% NUEVO JP
        if (maxa >= maxp)&&(maxp < umbralsig*absmin*Tmax_Tmin_bifasic),
            maxp = [];
        elseif (maxp> maxa)&&(maxa < umbralsig*absmin*Tmax_Tmin_bifasic),
            maxa = [];
        end
    end
    % Test which modulus maxima are significative and find zero crossings
    if isempty(maxa),
        if isempty(maxp),
            tipoT = 3;      % only downwards T wave
            if minpos - maxapos > 2, %!!!!!!!!
                ind = zerocros(flipud(w(maxapos:minpos,scalezerocros)));   %Scale scalezerocros
                if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                    ind = zerocros(flipud(w(maxapos:minpos,scalezerocros2))); 
                end %%%%%%%%%%%%%%Rute 03/09/02
                T = minpos - ind +1;
                picoff = minpos;
            elseif isempty(maxapos);  % If there were no maximum, there is no zero crossing.	
                T = picant (w(begwin:minpos,scale2),minpos);  % menimo en escala 4.
                if ~isempty(T) %%%%%%%%%%%% 14/06/02 Rute 
                    picoff = maxpos;
                else
                    picoff = []; % if did not exist a peak in scale 4 there is no T
                end %%%%%%%%%%%% 14/06/02 Rute 
            else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
            end
        else      % maxp is signficative, but not maxa
            tipoT = 1;  %inverted T wave
            if maxppos -minpos >2,  % !!!!!!
                ind = zerocros(w(minpos:maxppos,scalezerocros));  %% 07/06/02 Rute
                if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                    ind = zerocros(w(minpos:maxppos,scalezerocros2)); 
                end %%%%%%%%%%%%%%Rute 03/09/02
                T = minpos + ind -1;
                picon = minpos;		% For calculating onset and offset
                picoff = maxppos;
            else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
               messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
            end
        end
    else
        if isempty(maxp),   % maxa is significative, but not maxp
            tipoT = 0;       %normal T wave
            if minpos -maxapos >2,
                ind = zerocros(w(maxapos:minpos,scalezerocros)); %% 07/06/02 Rute
                if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                    ind = zerocros(w(maxapos:minpos,scalezerocros2)); 
                end %%%%%%%%%%%%%%Rute 03/09/02
                T = maxapos + ind -1;
                picon = maxapos;
                picoff = minpos;
            else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
               messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
            end
        else    % both maxa and maxp are significative.  Biphasic wave.
            tipoT = 4;	% biphasic pos-neg T wave
            if minpos - maxapos > 2,  %!!!!!!!!!!!
                % %% 07/06/02 Rute
                ind = zerocros(flipud(w(maxapos:minpos,scalezerocros))); %% 07/06/02 Rute
                if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                    ind = zerocros(flipud(w(maxapos:minpos,scalezerocros2)));
                end %%%%%%%%%%%%%%Rute 03/09/02
                T = minpos - ind +1;
                picon = maxapos;
            else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
            end
            if maxppos - minpos > 2,    %!!!!!!!!!!!!
                ind = zerocros(flipud(w(minpos:maxppos,scalezerocros))); %% 07/06/02 Rute
                if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                    ind = zerocros(flipud(w(minpos:maxpos,scalezerocros2)));
                end %%%%%%%%%%%%%%Rute 03/09/02
                Tprima = minpos + ind -1;
                picoff = maxppos;
            else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
            end
        end
    end
end

% T wave onset and offset detection
%if isempty(T), picon=[]; picoff=[]; end
if ~isempty(picon),
    Ton = searchon (picon, w(max(begwin,picon-round(T_bound_tol*messages.setup.wavedet.freq)):picon,scale2), Kton); 
end

if ~isempty(picoff),
    Toff=searchoff(picoff, w(picoff:min([size(w,1) picoff+round(T_bound_tol*messages.setup.wavedet.freq) ]) ,scale2) , Ktoff);
    if (Toff > endwin),
        Toff = endwin;
    end 
end

