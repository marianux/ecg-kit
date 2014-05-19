function [a,v] = optimline(pontos,weight,options)
% [a,v] = optimline(pontos)
% fit a line to a set of n points in Total least squares sense using unconstrained nonlinear optimization
%
% Input variable:
%   pontos -matrix nX3 with spatioal coordenates of points 
%
% Output variables:
%   a - sequence time occurence
%   v - vector director to the fitted line
%
% Rute Almeida 2.Dec.2004
% Last update: Rute Almeida 02Jun05
%
% Designed for MATLAB Version R12; tested with MATLAB Version R13

global loopfig
if nargin==1
    weight=ones(size(pontos,1),1);
    options=[];
end


%pontos=[TX(intreg(1):intreg(2)) TY(intreg(1):intreg(2)) TZ(intreg(1):intreg(2))];

k= -10:0.1:10;

%parameters_initial=[TX(intreg(1)) TY(intreg(1)) TZ(intreg(1)) TX(intreg(2))-TX(intreg(1)) TY(intreg(2))-TY(intreg(1)) TZ(intreg(2))-TZ(intreg(1))]%[a v]
%[x,fval,exitflag,output] = fminunc(@objfun,parameters_initial,options,pontos) 
% recta_initial=[repmat([TX(intreg(1)) TY(intreg(1)) TZ(intreg(1))],length(k),1)]+repmat([TX(intreg(2))-TX(intreg(1)) TY(intreg(2))-TY(intreg(1)) TZ(intreg(2))-TZ(intreg(1))],length(k),1).*[k' k' k'];
% recta_TLSnew=[repmat(x(1:3),length(k),1)]+[repmat(x(4:6),length(k),1)].*[k' k' k'];

%parameters_initial=[ TX(intreg(1))  TY(intreg(1))  TZ(intreg(1)) ;TX(intreg(2))-TX(intreg(1)) TY(intreg(2))-TY(intreg(1)) TZ(intreg(2))-TZ(intreg(1))]%[a v]
aux=2;
if isempty(pontos)
    a=NaN;
    v=NaN;
else
parameters_initial=[ pontos(1,:)  ; pontos(aux,:)-pontos(1,:)];%[a v]
if size(pontos,1)>2%length(pontos)>2 %RUTE 13.03.08
    while prod(parameters_initial(2,:))==0 && aux<=(length(pontos)-1)
        aux=aux+1;
        parameters_initial(2,:)=pontos(aux,:)-pontos(1,:);
    end
    if prod(parameters_initial(2,:))~=0 % admissible initial condition
        %[x,fval,exitflag,output] = fminunc(@objfun,parameters_initial,options,pontos,repmat(weight,1,size(pontos,2))); %  21.03.05
        [x,fun,funsrtq,exitflag,output] = lsqnonlin(@objfun2,parameters_initial,[],[],options,pontos,repmat(weight,1,size(pontos,2))); %#ok<ASGLU,NASGU> %  21.03.05
        % [x,fval,exitflag,output] = fminsearch(@objfun,parameters_initial,options,pontos,repmat(weight,1,size(pontos,2))); %  21.03.05
        while (exitflag<=0 || prod(parameters_initial(2,:))==0) && aux<=(length(pontos)-1)
            aux=aux+1;
            parameters_initial(2,:)=pontos(aux,:)-pontos(1,:);
            if prod(parameters_initial(2,:))==0
            parameters_initial(2,parameters_initial(2,:)==0)=0.0000001;
            end
            %[x,fval,exitflag,output] = fminunc(@objfun,parameters_initial,options,pontos,repmat(weight,1,size(pontos,2))); %  21.03.05
            [x,fun,funsrtq,exitflag,output] = lsqnonlin(@objfun2,parameters_initial,[],[],options,pontos,repmat(weight,1,size(pontos,2))); %#ok<ASGLU,NASGU> %  21.03.05
            %[x,fval,exitflag,output] = fminsearch(@objfun,parameters_initial,options,pontos,repmat(weight,1,size(pontos,2))); %  21.03.05
        end
        
    else
        exitflag=-1000;
        while sum(parameters_initial(2,:))==0 && aux<=(length(pontos)-1)
            aux=aux+1;
            parameters_initial(2,:)=pontos(aux,:)-pontos(1,:);
        end
    end
    a=parameters_initial(1,:);
    v=parameters_initial(2,:);
    if exitflag>0
        a=x(1,:);
        v=x(2,:);
        %recta_initial=[repmat([TX(intreg(2))-TX(intreg(1)) TY(intreg(2))-TY(intreg(1)) TZ(intreg(2))-TZ(intreg(1))],length(k),1).*[k' k' k']];
    end
else
    a=parameters_initial(1,:);
    v=parameters_initial(2,:);
    %recta_initial=[repmat([TX(intreg(2))-TX(intreg(1)) TY(intreg(2))-TY(intreg(1)) TZ(intreg(2))-TZ(intreg(1))],length(k),1).*[k' k' k']];
end
recta_initial=repmat(parameters_initial(2,:),length(k),1).*repmat(k',1,size(pontos,2));
recta_TLSnew=repmat(a,length(k),1)+repmat(v,length(k),1).*repmat(k',1,size(pontos,2));
if loopfig==1
    if size(pontos,2)==3
        h = figure;
        hs = axes(h);
        plot3(hs,pontos(:,1),pontos(:,2),pontos(:,3),'.g')
        hold(hs,'on');
        plot3(hs,0,0,0,'*c')
        limits=[get(hs,'XLim');get(hs,'YLim');get(hs,'ZLim')];
        plot3(hs,recta_initial(:,1),recta_initial(:,2),recta_initial(:,3),':')
        plot3(hs,recta_TLSnew(:,1),recta_TLSnew(:,2),recta_TLSnew(:,3),'r')
        legend(hs,'WT loop','zero','initial state','optimum line')
        %legend('WT4 loop','zero','search window end','Tend detected in single lead','Tand cardmean','loop in the interval','zero','initial state',['optimum line (H' num2str(H) ')'])
        set(hs,'XLim',limits(1,:))
        set(hs,'YLim',limits(2,:))
        set(hs,'ZLim',limits(3,:))
    elseif size(pontos,2)==2
        h = figure;
        hs = axes(h);
        plot(hs,pontos(:,1),pontos(:,2),'.k')
        hold(hs,'on')
        plot(hs,0,0,'*c')
        plot(hs,pontos(1,1),pontos(1,2),'s')
        plot(hs,pontos(end,1),pontos(end,2),'d')
        limits=[get(hs,'XLim');get(hs,'YLim')];
        plot(hs,recta_initial(:,1),recta_initial(:,2),':')
        plot(hs,recta_TLSnew(:,1),recta_TLSnew(:,2),'r')
        legend(hs,'WT loop','zero','interval begin','interval end','initial state','optimum line')
        %legend('WT4 loop','zero','search window end','Tend detected in single lead','Tand cardmean','loop in the interval','zero','initial state',['optimum line (H' num2str(H) ')'])
        set(hs,'XLim',limits(1,:))
        set(hs,'YLim',limits(2,:))
    end
end
end