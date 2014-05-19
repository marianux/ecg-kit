function [pid_starts, pid_ends] = TaskPartition( task_size, cant_pid)

%no es recomendable hacerlo mas grande de cant_recs la particion.
cant_pid = min(cant_pid, task_size);

things2do = fix(task_size / cant_pid);

remainder = rem(task_size, cant_pid);

cantThingsXpid = repmat(things2do, cant_pid,1);

cantThingsXpid(1:remainder) = cantThingsXpid(1:remainder) + 1;

pid_ends = cumsum(cantThingsXpid);
pid_starts = [1;pid_ends(1:end-1)+1];
