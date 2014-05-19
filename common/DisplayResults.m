function [dsResult ConfusionMat lablist global_performances ] = DisplayResults( varargin )

    p = inputParser;   % Create instance of inputParser class.
    
    p.addParamValue('dsResult', [], @(x)( isnumeric(x) | isa(x,'prdataset') ) );
    p.addParamValue('bWithFeatMat', true, @islogical);
    p.addParamValue('priors', [], @isnumeric);
    p.addParamValue('TrueLabels', [], @(x)(isnumeric(x) || ischar(x)));
    p.addParamValue('ClassifiedLabels', [], @(x)(isnumeric(x) || ischar(x)));
    p.addParamValue('datasetName', [], @ischar);
    p.addParamValue('ClassLabels', [], @(x)( ischar(x) ));
    p.addParamValue('class_filter', [], @(x)( ischar(x) ));
    p.addParamValue('SupportDataset', [], @(x)( isa(x,'prdataset') ));
    p.addParamValue('TrainedMapping', [], @(x)( isa(x,'prmapping') && istrained(x)) );
    
    
    try
        p.parse( varargin{:});
    catch MyError
        rethrow(MyError);    
    end

    dsResult = p.Results.dsResult;
    bWithFeatMat = p.Results.bWithFeatMat;
    priors = p.Results.priors;
    TrueLabels = p.Results.TrueLabels;
    ClassifiedLabels = p.Results.ClassifiedLabels;
    datasetName = p.Results.datasetName;
    ClassLabels = p.Results.ClassLabels;
    SupportDataset = p.Results.SupportDataset;
    w_TrainedMapping = p.Results.TrainedMapping;
    class_filter = p.Results.class_filter;

    
    labels_not_learned_idx = [];
    cant_iter = 1;
    global_performances = [];
    
    if( isdataset(dsResult)  )

        if( ~isempty(w_TrainedMapping) )
            dsTrain_lablist = getlabels(w_TrainedMapping);
            dsTest_lablist = getlablist(dsResult);
            [labels_not_learned labels_not_learned_idx ] = setdiff(dsTest_lablist, dsTrain_lablist, 'rows');

            dsResult = dsResult * w_TrainedMapping;
            %Hago esto para que no batchee ...
%             dsResult = map(dsResult, w_TrainedMapping, 500e3);
            
        else
            dsTest_lablist = getlablist(dsResult);
            [labels_not_learned labels_not_learned_idx ] = setdiff(dsTest_lablist, ClassLabels, 'rows');
        end
                
        
        if(bWithFeatMat)
%             confmat(dsResult)
            confmat_new(dsResult)
        end

        if( isempty(datasetName) )
            datasetName = getname(dsResult);
        end

        [ConfusionMat, ~,lablist] = confmat_new(dsResult);
%         [ConfusionMat, ~, lablist] = confmat(dsResult);
        
        if( size(lablist,1) < size(ConfusionMat,2) )
            %hay clase reject
            bRejectClass = true;
            rejected = ConfusionMat(:,end);
            not_labeled = ConfusionMat(end,:);
            ConfusionMat = ConfusionMat(1:end-1,1:end-1);
        else
            bRejectClass = false;
        end
        
        priors = getprior(dsResult);

    elseif( isempty(dsResult) && ~isempty(TrueLabels) && ~isempty(ClassifiedLabels) )
        
        if( isempty(ClassLabels) )
            confmat_new(TrueLabels, ClassifiedLabels);
            [ConfusionMat, ~, lablist] = confmat_new(TrueLabels, ClassifiedLabels);
        else
            [labels_not_learned labels_not_learned_idx ] = setdiff( unique(TrueLabels, 'rows'), ClassLabels, 'rows');
            
            if( ischar(TrueLabels) && ischar(ClassifiedLabels) )
                confmat_new(TrueLabels, ClassifiedLabels, 'count', 1, ClassLabels );
                [ConfusionMat, ~, lablist] = confmat_new(TrueLabels, ClassifiedLabels, 'count', 1, ClassLabels );
            elseif( isnumeric(TrueLabels) && isnumeric(ClassifiedLabels) )
                confmat_new(ClassLabels(TrueLabels,:), ClassLabels(ClassifiedLabels,:));
                [ConfusionMat, ~, lablist] = confmat_new(ClassLabels(TrueLabels,:), ClassLabels(ClassifiedLabels,:) );
            else
                error('Tipo de datos no permitidos en las etiquetas.' )
            end                
        end
        
        if( size(lablist,1) < size(ConfusionMat,2) )
            %hay clase reject
            bRejectClass = true;
            rejected = ConfusionMat(:,end);
            not_labeled = ConfusionMat(end,:);
            ConfusionMat = ConfusionMat(1:end-1,1:end-1);
        else
            bRejectClass = false;
        end
        
    elseif( isnumeric(dsResult) )

        ConfusionMat = dsResult;

        iCantClases = size(ConfusionMat,2);

        bRejectClass = false;
            
        if( ~isempty(SupportDataset) )

            datasetName = getname(SupportDataset);

            priors = getprior(SupportDataset);

            lablist = getlablist(SupportDataset);                

        else

            if( isempty(priors) )
                %Equal priors
                warning('Assumnig equal priors.')
                priors = ones(iCantClases,1)/iCantClases;
            end

            if( isempty(ClassLabels) )
                lablist = [ repmat('Class ', iCantClases, 1)  num2str((1:iCantClases)')];
            else
                lablist = ClassLabels;
            end
            
        end

        cant_iter = size(ConfusionMat,3);
        
        DisplayConfusionMatrix(ConfusionMat, lablist);
        
    else
        error('Argumentos no reconocidos');
    end

    %filtro de clases
    if( isempty(class_filter) )
        %todas las clases
        class_filter_idx = 1:size(lablist,1);
        
    else
        
        [~, class_filter_idx ] = intersect(cellstr(lablist), cellstr(class_filter));

        class_filter_idx = rowvec(class_filter_idx);
        
%         ConfusionMat_orig = ConfusionMat;
%         lablist_orig = lablist;
% 
%         cant_clases_filt = length(class_filter_idx);
%         ConfusionMat = zeros(cant_clases_filt);
%         
%         for ii = 1:cant_clases_filt
%             ConfusionMat(ii,:) = ConfusionMat_orig(class_filter_idx(ii), class_filter_idx);
%             ConfusionMat(:,ii) = ConfusionMat_orig(class_filter_idx, class_filter_idx(ii));
%         end
%         lablist = lablist_orig(class_filter_idx,:);       
        

    end
    

    %Las señalo ya que requieren un tratamiento especial segun AAMI-EC57.
    Unknown_idx = find(strcmpi(cellstr(lablist),'Unknown'));
    Fusion_idx = find(strcmpi(cellstr(lablist),'Fusion'));
    Supra_idx = find(strcmpi(cellstr(lablist),'Supraventricular'));
    Ventri_idx = find(strcmpi(cellstr(lablist),'Ventricular'));
    
    
%     %Ignoramos la clase Unknown para el calculo de resultados si existe.
%     Unknown_idx = find(strcmpi(cellstr(lablist),'Unknown'));
%     
%     if( ~isempty(Unknown_idx) )
%         %Elimino esa clase.
%         lablist_orig = lablist;
%         lablist(Unknown_idx,:) = [];
%         ConfusionMat_orig = ConfusionMat;
%         ConfusionMat(Unknown_idx,:,:) = [];
%         ConfusionMat(:,Unknown_idx,:) = [];
%     end
    
    %Corrregiremos la matriz de confusion si hay clases no aprendidas.
    Cant_real_x_clase = squeeze(sum(ConfusionMat,2));        
    Cant_predecida_x_clase = squeeze(sum(ConfusionMat,1));
    ClasesPresentes_idx = find(Cant_real_x_clase > 0);
    ClasesPredichas_idx = find(Cant_predecida_x_clase > 0);

    if( ~isempty(labels_not_learned_idx) )
        %las clases no aprendidas, hago como que no estan presentes para la contabilizacion de
        %resultados.
        ClasesPresentes_idx(labels_not_learned_idx) = [];

        %y por otro lado agrupo en la diagonal las clasificaciones de las
        %clases que no fueron aprendidas.

% esto es por si quisiera contabilizarlo como acierto, creo que
% directamente hay que ignorarlo
%         for ii = find((ClasesPredichas_idx(:))')
%             ConfusionMat(ii,ii) = ConfusionMat(ii,ii) + ConfusionMat(labels_not_learned_idx,ii);
%             ConfusionMat(labels_not_learned_idx,ii) = 0;
%         end
        
        %anulo las clases no entrenadas para que no contabilicen en los
        %resultados.
        ConfusionMat(labels_not_learned_idx,:) = 0;
        Cant_real_x_clase = sum(ConfusionMat,2);        
        Cant_predecida_x_clase = sum(ConfusionMat)';
        ClasesPresentes_idx = find(Cant_real_x_clase ~= 0);
        ClasesPredichas_idx = find(Cant_predecida_x_clase ~= 0);

    end

    %Una vuelta para mostrar los resultados balanceados y otra
    %desbalanceados.
    for jj = 1:2

        performance = [];
        
        if(jj == 1)
            strAux = ['Balanced Results for ' datasetName];
        else
            strAux = ['Unbalanced Results for ' datasetName];
        end

        disp(strAux)
        disp( repmat('-',1,length(strAux)) );

        if(bRejectClass)
            strPerfMeasures = '|  Se   +P   Rej  |';
        else
            if( cant_iter == 1 )
            %     strPerfMeasures = '|  Se   +P  FPR  |';
                strPerfMeasures = '|  Se   +P  |';
            else
                strPerfMeasures = '|     Se       +P     |';
            end
        end
        
        iFieldWidth = length(strPerfMeasures) - 4;
        iLablistWidth = size(lablist,2);

        for ii = class_filter_idx
            fprintf(1, [ '| ' lablist(ii,1:min(iFieldWidth, iLablistWidth)) repmat(' ', 1, iFieldWidth-iLablistWidth) ' |'  ]);
        end
        
        if(bRejectClass)
            fprintf(1, [ '|                TOTALS                 |\n'  ]);
        else    
            if( cant_iter == 1 )
                fprintf(1, [ '|           TOTALS            |\n'  ]);
            else
                fprintf(1, [ '|                   TOTALS                   |\n'  ]);
            end
        end
        
        for ii = class_filter_idx
            fprintf(1, strPerfMeasures);
        end

        if(bRejectClass)
            fprintf(1, [ '|   Acc   |   Se    |   +P    |  Reject |\n'  ]);
        else
            if( cant_iter == 1 )
                fprintf(1, [ '|   Acc   |   Se    |   +P    |\n'  ]);
            else
                fprintf(1, [ '|     Acc      |     Se       |      +P      |\n'  ]);
            end
        end

        ConfusionMat_bal = zeros(size(ConfusionMat));

        if(jj == 1)
            %Balanceamos las clases para que las mediciones de performance no
            %tengan sesgo
            k_x_clase = repmat(max(Cant_real_x_clase, [], 1), size(Cant_real_x_clase, 1), 1) ;
            k_x_clase(ClasesPresentes_idx) =  k_x_clase(ClasesPresentes_idx) ./ Cant_real_x_clase(ClasesPresentes_idx);
            if( cant_iter == 1)
                ConfusionMat_bal = ConfusionMat .* repmat( k_x_clase, 1, size(ConfusionMat,2) );
            else
                ConfusionMat_bal = cellfun( @(a,b)(round( a .* repmat( b, 1, size(a,2) ))), mat2cell(ConfusionMat, size(ConfusionMat,1), size(ConfusionMat,2), ones(1,size(ConfusionMat,3))) , reshape(mat2cell(k_x_clase, size(k_x_clase,1), ones(1,size(k_x_clase,2)) ), [1 1 size(k_x_clase,2)] ) , 'UniformOutput', false);
                ConfusionMat_bal = cell2mat(ConfusionMat_bal);
            end
            if(bRejectClass)             
                rejected_bal = rejected(ClasesPresentes_idx) .* k_x_clase;
            end
        else
%             ConfusionMat_bal(ClasesPresentes_idx,:) = ConfusionMat(ClasesPresentes_idx,:);        
            ConfusionMat_bal = ConfusionMat;        
            if(bRejectClass)             
                rejected_bal = rejected(ClasesPresentes_idx);
            end
        end

        for ii = class_filter_idx

            if( sum(Cant_real_x_clase(ii,:)) == 0 )
                if(bRejectClass)    
                    fprintf(1,   '|  --   --   --  |' );
                    aux = nan(1,3);
                    
                else
                    if( cant_iter > 1 )
                        fprintf(1,   '|  --( --)%%  --( --)%% |' );
                        aux = nan(1,4);
                    else
                        fprintf(1,   '|  --   --  |' );
                        aux = nan(1,2);
                    end
                end
            else
                
                
                iTP = squeeze(ConfusionMat_bal(ii,ii,:));
                iFN = colvec(sum(squeeze(ConfusionMat_bal(ii,:,:)))) - iTP;
                PosPred = colvec(sum(squeeze(ConfusionMat_bal(:,ii,:))));
    
                % AAMI special treatment for P+ calculation
                if( ~isempty(Supra_idx) && ~isempty(Unknown_idx) && Supra_idx == ii )
                    %Descuento los unknown
                    PosPred = PosPred - colvec(squeeze(ConfusionMat_bal(Unknown_idx,ii,:)));
                end
                
                if( ~isempty(Ventri_idx) && Ventri_idx == ii )
                    if( ~isempty(Unknown_idx) )
                        %Descuento los unknown
                        PosPred = PosPred - colvec(squeeze(ConfusionMat_bal(Unknown_idx,ii,:)));
                    end
                    if( ~isempty(Fusion_idx) )
                        %Descuento los unknown
                        PosPred = PosPred - colvec(squeeze(ConfusionMat_bal(Fusion_idx,ii,:)));
                    end
                end
                
                if( PosPred ~= 0)
                    PosPred = iTP./PosPred*100;
                else
                    PosPred = zeros(1,cant_iter);
                end
%                 iFP = sum(ConfusionMat_bal(:,ii)) - iTP;
%                 iTN = sum(sum(ConfusionMat_bal)) - sum(ConfusionMat_bal(ii,:)) - sum(ConfusionMat_bal(:,ii)) + iTP;

                if(bRejectClass)             
                    aux = [iTP/(iTP+iFN)*100, PosPred, rejected_bal(ii)/(rejected_bal(ii)+iTP+iFN)*100 ];
                    fprintf(1, [ '| %3.0f%% %3.0f%% %3.0f%%  |'  ], aux);
                    
                else
                    if( cant_iter > 1 )
                        aux = [ median(iTP./(iTP+iFN))*100, mad(iTP./(iTP+iFN),1)*100, median(PosPred), mad(PosPred,1) ];
                        fprintf(1, [ '| %3.0f(%3.0f)%% %3.0f(%3.0f)%% |'  ], aux );
                    else
                        aux = [ iTP/(iTP+iFN)*100, PosPred ];
                        fprintf(1, [ '| %3.0f%% %3.0f%% |'  ], aux );
                    end
                end
                
            end
            
            performance = [ performance aux];
            
        end

        priors = priors(:);

        if( cant_iter == 1)
            Aciertos = diag(ConfusionMat_bal);
        else
            Aciertos = cellfun(@(a)(diag(a)), mat2cell(ConfusionMat_bal, size(ConfusionMat_bal,1), size(ConfusionMat_bal,2), ones(1,size(ConfusionMat_bal,3))), 'UniformOutput', false );
            Aciertos = squeeze(cell2mat(Aciertos));
        end
        

        Cant_real_x_clase_bal = colvec(squeeze(sum(ConfusionMat_bal,2)));        
        Cant_predecida_x_clase_bal = colvec(squeeze(sum(ConfusionMat_bal)));
        ClasesPresentes_bal_idx = find(Cant_real_x_clase_bal ~= 0);
        ClasesPredichas_bal_idx = find(Cant_predecida_x_clase_bal ~= 0);
        
        
        
        if(bRejectClass)             

            aux = [ ...
                    (sum(Aciertos(ClasesPresentes_bal_idx))/sum(sum(ConfusionMat_bal)))*100, ... 
                    mean(Aciertos(ClasesPresentes_bal_idx)./Cant_real_x_clase_bal(ClasesPresentes_bal_idx))*100,  ... 
                    mean(Aciertos(ClasesPredichas_bal_idx)./Cant_predecida_x_clase_bal(ClasesPredichas_bal_idx))*100, ...
                    sum(rejected_bal)/(sum(rejected_bal) + sum(sum(ConfusionMat_bal)))*100 ];
            
            fprintf(1, [ '|  %3.0f%%   |  %3.0f%%   |  %3.0f%%   |  %3.0f%%   |'  ], aux);
            
        else
            if( cant_iter > 1 )
                
                aux_s = nan(size(Aciertos));
                aux_s(ClasesPresentes_bal_idx) = Aciertos(ClasesPresentes_bal_idx)./Cant_real_x_clase_bal(ClasesPresentes_bal_idx);
                aux_s = nanmean(aux_s,1);
                
                aux_p = nan(size(Aciertos));
                aux_p(ClasesPresentes_bal_idx) = 0;
                aux_p(ClasesPredichas_bal_idx) = colvec(Aciertos(ClasesPredichas_bal_idx)./Cant_predecida_x_clase_bal(ClasesPredichas_bal_idx));
... %                       %funciona mal cuando no predice una clase
%                 aux_p(ClasesPresentes_bal_idx) = Aciertos(ClasesPresentes_bal_idx)./Cant_predecida_x_clase_bal(ClasesPresentes_bal_idx);
                
                aux_p = nanmean( aux_p ,1);
                
                aux = [ ...
                    median(colvec(sum(Aciertos))./squeeze(sum(sum(ConfusionMat_bal,1),2)))*100, ... 
                    mad(colvec(sum(Aciertos))./squeeze(sum(sum(ConfusionMat_bal,1),2)),1)*100, ... 
                    median(aux_s)*100,  ... 
                    mad(aux_s,1)*100,  ... 
                    median(aux_p)*100,  ... 
                    mad(aux_p,1)*100 ];
                                                        
                fprintf(1, [ '|  %3.0f(%3.0f)%%   |  %3.0f(%3.0f)%%   |  %3.0f(%3.0f)%%   |'  ], aux);
                
            else
                
                aux = [ ...
                        (sum(Aciertos(ClasesPresentes_bal_idx))/sum(sum(ConfusionMat_bal)))*100, ... 
                        mean(Aciertos(ClasesPresentes_bal_idx)./Cant_real_x_clase_bal(ClasesPresentes_bal_idx))*100,  ... 
                        mean([Aciertos(ClasesPredichas_bal_idx)./Cant_predecida_x_clase_bal(ClasesPredichas_bal_idx); zeros(abs(length(ClasesPresentes_bal_idx)-length(ClasesPredichas_bal_idx)),1) ])*100 ...
... %                       %funciona mal cuando no predice una clase
... %                         mean(Aciertos(ClasesPresentes_bal_idx)./Cant_predecida_x_clase_bal(ClasesPresentes_bal_idx))*100 ...
                        ];
                
                fprintf(1, [ '|  %3.0f%%   |  %3.0f%%   |  %3.0f%%   |'  ], aux );
            end      
            
        end
        
        performance = [ performance aux];
        
        global_performances(jj,:) = performance;
        
        fprintf(1, '\n\n');

        
    end
    
    if( bRejectClass && any( not_labeled > 0 ) )
       warning( ['Se encontraron ' num2str( sum(not_labeled) ) ' datos no etiquetados.' ] ) 
    end
    
%     if( ~isempty(Unknown_idx) )
%         %Vuelvo al estado original la CM y el lablist.
%         lablist = lablist_orig ;
%         ConfusionMat = ConfusionMat_orig ;
%     end
    
    
    %% Debug, borrarme
%     fprintf(1, 'Debug:\n');
% 
%     fprintf(1, repmat([repmat('%3.0f\t', 1, size(performance,2)) '\n'], 1,size(performance,1)), performance' );
    
    %% FIN Debug, borrarme
                                                            

        