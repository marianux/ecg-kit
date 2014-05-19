function dsOut = deNaN_dataset(dsIn, type)
    
    if( nargin < 2 || isempty(type) )
        type = 'remove';
    end
    
    
    features_with_nans = find(any(isnan(+dsIn)));
    rows_without_nans = ~any(isnan(+dsIn),2);
    
    if( any(~rows_without_nans) )

        
        rows_without_nans_idx = find(rows_without_nans);
        
        switch(type)

            case 'remove'
                % just remove nans

                dsOut = dsIn(rows_without_nans_idx,:);

            case 'change'

                % change NaNs for repeated values of the same feature, added an small variance, estimated ignoring this nan
                % values.
                dsOut = dsIn;
                

                robust_std = nanmeda(+dsIn(rows_without_nans_idx,:));

                for kk = rowvec(features_with_nans)
                    
                    nan_idx = find(isnan(+(dsIn(:,kk))));
                    not_nan_idx = find(~isnan(+(dsIn(:,kk))));
                    laux_idx = length(nan_idx);
                    
                    dsIn(nan_idx,kk) = randsample( +dsIn(not_nan_idx,kk), laux_idx, true ) + 1/100*robust_std(kk)*rand(laux_idx,1);
                    
                    dsOut(:,kk) = dsIn(:,kk);

                end

            case 'change_same_class'

                % change NaNs for repeated values of the same feature, same
                % class, added an small variance, estimated ignoring this nan
                % values.
                dsOut = dsIn;
                
                dsIn = seldat(dsIn);

                ds_aux = dsIn(rows_without_nans_idx,:);
                ds_aux = setprior(ds_aux, 0);

                w_aux = qdc(ds_aux);
                w_aux = +w_aux;
                Class_indexes = getnlab(dsIn);        

                for kk = rowvec(features_with_nans)
                    nan_idx = find(isnan(+(dsIn(:,kk))));

                    classes_involved = unique(Class_indexes(nan_idx));

                    for ll = rowvec(classes_involved)

                        aux_idx = find(Class_indexes(nan_idx) == ll);
                        this_class_not_nan_idx = find(Class_indexes == ll);
                        this_class_not_nan_idx = setdiff(this_class_not_nan_idx, nan_idx(aux_idx) );

                        laux_idx = length(aux_idx);

                        dsIn(nan_idx(aux_idx),kk) = randsample( +dsIn(this_class_not_nan_idx,kk), laux_idx, true ) + 1/100*w_aux.cov(kk,kk,ll)*rand(laux_idx,1);

                    end

                    dsOut(:,kk) = dsIn(:,kk);

                end

        end
    else
        
        dsOut = dsIn;
        
    end
    
    