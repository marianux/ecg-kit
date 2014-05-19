function meda = nanmeda(data)

median_data = nanmedian(data);
meda = nanmedian(abs(bsxfun(@minus, data, median_data)));
