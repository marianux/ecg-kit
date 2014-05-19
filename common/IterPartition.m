function [pid_starts, pid_ends] = IterPartition( task_size, iter_size)

%no es recomendable hacerlo mas grande de cant_recs la particion.
iter_size = min(iter_size, task_size);

iters2do = floor(task_size / iter_size);
iter_size = fix(task_size / iters2do);

remainder = rem(task_size, iters2do);

cantThingsXpid = repmat(iter_size, iters2do,1);

cantThingsXpid(1:remainder) = cantThingsXpid(1:remainder) + 1;

pid_ends = cumsum(cantThingsXpid);
pid_starts = [1;pid_ends(1:end-1)+1];
