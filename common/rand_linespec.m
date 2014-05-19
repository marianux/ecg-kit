function line_spec = rand_linespec()

mrk={'+','o','*','.','x','s','d','^','v','<','>','p','h'}.';
linestyle = {'-','--',':','-.'};
colspec={'m','c','r','g','b','k'}.';

line_spec = [linestyle{randsample(1:length(linestyle),1)} colspec{randsample(1:length(colspec),1)} mrk{randsample(1:length(mrk),1)} ];
