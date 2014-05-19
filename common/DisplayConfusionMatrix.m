function DisplayConfusionMatrix( C, lablist )

cant_iter = size(C,3);

labcv = lablist;
n1 = size(lablist,1);
n2 = n1;


if( cant_iter > 1 )
    
    C_mean = round(mean(C,3));
    C_std = round(std(C,0,3));
    C_sum2 = sum(C,2);
    C_sum1 = sum(C,1);
    
    max_num = length(sprintf( '%i', max(max([ C_mean; round(mean(sum(C,1),3)) ] )) ));
    max_num = max_num + 2 + length(sprintf( '%i', max(max([C_std; round(std(C_sum1,0,3))])) ));

else
    
    max_num = length(sprintf( '%i', max(max([C; sum(C) ])) ));
    
end


if (size(lablist,2) > max_num)
    lablist = lablist(:,1:max_num); 
    %labcv = labcv(:,1:6); 
end
if (size(lablist,2) < 5)
    lablist = [lablist repmat(' ',n2,ceil((5-size(lablist,2))/2))]; 
    labcv = [labcv repmat(' ',n1,ceil((5-size(labcv,2))/2))]; 
end

nspace = max(size(labcv,2)-8,0);
cspace = repmat(' ',1,nspace);
%fprintf(1,['\n' cspace '        | Estimated Labels']);
fprintf(1,['\n  True   ' cspace '| Estimated Labels']);
fprintf(1,['\n  Labels ' cspace '| ']);
for j = 1:n2, fprintf(1,'%s ',lablist(j,:)); end
fprintf(1,'|');
fprintf(1,' Totals\n');
fprintf(1,' %s|%s|%s\n', repmat('-',1,8+nspace), repmat('-',1,1+(max_num+1)*n2),repmat('-',1,max_num));


for j = 1:n1
    fprintf(1,' %-7s|',labcv(j,:));
    if( cant_iter > 1 )
        c_aux = colvec([ colvec(C_mean(j,:)) colvec(C_std(j,:)) ]');
        this_row = ' ';
        for kk = 1:2:length(c_aux)
            str_aux = sprintf('%i(%i)', c_aux(kk), c_aux(kk+1));
            str_aux = [str_aux repmat(' ', 1, max(0, 1 + max_num - length(str_aux))) ];
            this_row = [this_row str_aux];
        end
        fprintf(1,'%s', this_row);
        fprintf(1,'| %i(%i)\n', round(mean(C_sum2(j,:,:),3)), round(std(C_sum2(j,:,:),0,3)) );
    else
        this_row = ' ';
        for kk = 1:n2
            str_aux = sprintf('%i', C(j,kk));
            str_aux = [str_aux repmat(' ', 1, max(0, 1 + max_num - length(str_aux))) ];
            this_row = [this_row str_aux];
        end
        fprintf(1,'%s', this_row);
        fprintf(1,'| %i\n',sum(C(j,:)));
    end
end

fprintf(1,' %s|%s|%s\n', repmat('-',1,8+nspace), repmat('-',1,1+(max_num+1)*n2),repmat('-',1,max_num));

fprintf(1,['  Totals ' cspace '|']);

if( cant_iter > 1 )
    c_aux = colvec([ colvec(round(mean(sum(C,1),3))) colvec(round(std(C_sum1,0,3))) ]');

    this_row = ' ';
    for kk = 1:2:length(c_aux)
        str_aux = sprintf('%i(%i)', c_aux(kk), c_aux(kk+1));
        str_aux = [str_aux repmat(' ', 1, max(0, 1 + max_num - length(str_aux))) ];
        this_row = [this_row str_aux];
    end
    fprintf(1,'%s', this_row);
    
    fprintf(1,'| %i(%i)\n\n', round(mean(sum(C_sum1,2),3)), round(std(sum(C_sum1,2),0,3)) );
else
    this_row = ' ';
    for kk = 1:n2
        str_aux = sprintf('%i', sum(C(:,kk)) );
        str_aux = [str_aux repmat(' ', 1, max(0, 1 + max_num - length(str_aux))) ];
        this_row = [this_row str_aux];
    end
    fprintf(1,'%s', this_row);
    fprintf(1,'| %5i\n\n',sum(C(:)));
end
