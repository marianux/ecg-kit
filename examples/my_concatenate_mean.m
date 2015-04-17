function payload = my_concatenate_mean(plA, plB)

if( isempty(plA) )
    payload = plB;
else
    payload.the_sum = plA.the_sum + plB.the_sum;
    payload.the_size = plA.the_size + plB.the_size;
end
