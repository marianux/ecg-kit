function [q1, q2, q3, q4,q5, messages] = qspfilt5(fs,messages) 
% This function obtains the quadratic splines wavelet filterbank filters from
% scale 1 to 4 as a function of the sampling frequency, in order to use filters
% with similar analog frequency behaviour for diferent sampling frecuencies.
%
% Juan Pablo Martínez Cortés
%
% Last update: Rute Almeida 24/11/11
%
% Designed for MATLAB Version R12; tested with MATLAB Version R13
%
if nargin<2 
messages = new_empty_mesg;
end
messages.status=1;
% beutors h(n) y g(n)
h = 0.125*[1 3 3 1];   %high pass % paso baxo
g = 2*[1 -1];          %low pass  % paso alto
h2 = kron(h,[1 0]);
g2 = kron(g,[1 0]);
h4 = kron(h,[1 zeros(1,3)]);
g4 = kron(g,[1 zeros(1,3)]);
h8 = kron(h,[1 zeros(1,7)]);
g8 = kron(g,[1 zeros(1,7)]);
g16=kron(g,[1 zeros(1,15)]);  % Rute 21/05/02 %% 
q1 = g;                                                                         % Q1(w) = G(w);
q2 = conv(h,g2);                                                        % Q2(w) = H(w)G(2w);
q3 = conv(conv(h,h2),g4);                               % Q3(w) = H(w)H(2w)G(4w);
q4 = conv(conv(conv(h,h2),h4),g8);      % Q4(w) = H(w)H(2w)H(4w)G(8w);
q5 = conv(conv(conv(conv(h,h2),h4),h8),g16);      % Q5(w) = H(w)H(2w)H(4w)H(8w)G(16w); % Rute 21/05/02 %%

l1 = max(find(q1~=0)); %#ok<MXFND>
l2 = max(find(q2~=0)); %#ok<MXFND>
l3 = max(find(q3~=0)); %#ok<MXFND>
l4 = max(find(q4~=0)); %#ok<MXFND>
l5 = max(find(q5~=0));  %#ok<MXFND>  % Rute 21/05/02 %%
q1 = q1(1:l1)';
q2 = q2(1:l2)';
q3 = q3(1:l3)';
q4 = q4(1:l4)';
q5 = q5(1:l5)'; % Rute 21/05/02 %%
% nao foi alterado para adaptar q5 a outras frequencias de amostragem ainda
% Adautazión de os filtros de frec. 250 Hz ta atras frecuenzias de muestreo
% Adaptation of filters from freq. 250 Hz to other sampling frequencies
% 500 Hz
if fs~=250,
if fs == 500,
   q1 = interp([0;q1;0],2,1); q1 = q1(2:end-2);
   q2 = interp([0;q2;0],2,3); q2 = q2(2:end-2);      
   q3 = interp([0;q3;0],2,7); q3 = q3(2:end-2);
   q4 = interp([0;q4;0],2,7); q4 = q4(2:end-2);
   q5 = interp([0;q5;0],2,7); q5 = q5(2:end-2);
elseif (fs == 1000),
   q1 = interp([0;q1;0],4,1); q1 = q1(2:end-4);
   q2 = interp([0;q2;0],4,3); q2 = q2(2:end-4);      
   q3 = interp([0;q3;0],4,7); q3 = q3(2:end-4);
   q4 = interp([0;q4;0],4,7); q4 = q4(2:end-4);
   q5 = interp([0;q5;0],4,7); q5 = q5(2:end-4);
elseif fs == 200,
   %q1 = interp([0;q1;0],4,1,0.4); q1 = q1(2:end-4);
   q2 = interp([0;q2;0],4,3,0.4); q2 = q2(2:end-4);      
   q3 = interp([0;q3;0],4,7,0.4); q3 = q3(2:end-4);
   q4 = interp([0;q4;0],4,7,0.4); q4 = q4(2:end-4);
   q5 = interp([0;q5;0],4,7,0.4); q5 = q5(2:end-4);
% No cal filtrar-ie pos ya ye filtrau a 100 Hz.
   %q1 = q1(1:5:end);
   q2 = q2(4:5:end);
   q3 = q3(5:5:end);
   q4 = q4(2:5:end);
   q5 = q5(2:5:end);
   q1 = 1.1*[5/4, -5/4];
elseif fs == 360,
   % Interpolo por 6
   q1 = interp([0;q1;0],6,1); q1 = q1(2:end-6);
   q2 = interp([0;q2;0],6,3); q2 = q2(2:end-6);      
   q3 = interp([0;q3;0],6,7); q3 = q3(2:end-6);
   q4 = interp([0;q4;0],6,7); q4 = q4(2:end-6);
   q5 = interp([0;q5;0],6,7); q5 = q5(2:end-6);
   % Diezmo por 5
   q1=q1(4:5:end);
   q2=q2(1:5:end);
   q3=q3(5:5:end);
   q4=q4(3:5:end);
   q5=q5(3:5:end);
   
   % Interpolo por 6
   q1 = interp([0;q1;0],6,1); q1 = q1(2:end-6);
   q2 = interp([0;q2;0],6,3); q2 = q2(2:end-6);      
   q3 = interp([0;q3;0],6,7); q3 = q3(2:end-6);
   q4 = interp([0;q4;0],6,7); q4 = q4(2:end-6);
   q5 = interp([0;q5;0],6,7); q5 = q5(2:end-6);
   %Diezmo por 5
   q1 = q1(2:5:end);
   q2 = q2(5:5:end);
   q3 = q3(4:5:end);
   q4 = q4(4:5:end);
   q5 = q5(4:5:end);
else
    messages.errors=[messages.errors {'There are no wavelets designed for this sampling frequency.'}];
%     warning(char(messages.errors(end)))
    messages.errors_desc=[messages.errors_desc 'Default wavelet filters not defined for thia samplig frequency. Use filter design instead.'];
    messages.status=0;
end
end  
