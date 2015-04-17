function result_payload = my_finish_mean(payload, ECG_header)

result_payload.mean = payload.the_sum ./ payload.the_size;
