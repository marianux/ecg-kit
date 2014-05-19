
fprintf( 1, '\n%s\n', date );
fprintf( 1, 'Configuration\n' );
fprintf( 1, '-------------\n' );
fprintf( 1, '+ Recording: %s (%s)\n',  recording_name, recording_format );
fprintf( 1, '+ Mode: %s (%d clusters, %d iterations, %d%% cluster-presence)\n', op_mode, CantClusters, iter_times, cluster_presence );

