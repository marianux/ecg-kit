function colors = my_colormap( cant_colors )
colors = hsv(64);

%avoid yellowish colours
bAux = colors(:,1) > 200/255 & colors(:,2) > 200/255 & colors(:,3) == 0;
if( any(bAux) )
    colors( bAux,: ) = repmat([200 200 0]./255, sum(bAux),1 );
end

ncol = size(colors,1);

aux_idx = 1+rem(round(round(ncol/6)+linspace(0, (round(ncol/2)+1)*(ncol-1), ncol)), ncol);
if( cant_colors > 4 )
    colors = colors( aux_idx(round(linspace(1,64, cant_colors))) , :);
else
    colors = colors( aux_idx(1:cant_colors) , :);
end
