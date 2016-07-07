# Todo: Newest near top
## Task power starting point
When we schedule a task, we need to know how much power the addition of a single task will cost the runtime. That means that we consider the gradient along the pareto frontier, which gives us the configuration for the quantized decrease in power.

However, these frontiers have been constructed in isolation; ie, each task is run by itself. The delta is important, but so is the absolute increase over not running the task at all (ie having zero threads running). This feeds into the energy scaling work, where we consider how much power is consumed by the device being ON, but idle.

We should run the extrapolated idle power experiments to determine this baseline, so that we know how much _additional_ power is required to run the task, on top of an idle system.

In addition, we need to take into consideration that baseline power consumption is always present. Perhaps we can subtract this value from the power cap, indicating the amount of power headroom we have to work with, assuming a powered on and running system? It definitely detracts from the race to idle camp, but from what I've heard, people say that they run at full power mode all the time anyway since they never know when things will be needed. This is true especially when you consider that, as it stands now, you can't really turn cores off. If this capability were possible, we'd have to work that into the model (ie, not just the dynamic power cost, but the incremental power cost to take the core out of a low power state).

The so each thread may be able to do more work (sqrt may not use all the core functional units) but since it's consistent, it should scale correctly (ie, the y-intercept should be correct).

We need to fix up the model data for spmv.csv, it's borked (ie, it's missing energy, power, and speedup).

### Draft data
Looks like it's about 107 W as the baseline.

##
Do I need to have a (time-based) performance model to make my scheduling decision? Would it make sense to construct some performance model, either as a best-case (ie oracle) comparison, or just as a baseline to show how much perf we give up, but how bad we blow our budget?

Is there a PAPI power limitation framework for the MIC? Nope. Need to do something more ad-hoc. Similar to power steering, we set one as a baseline? Force it as best as possible, and document how far off we are.

Are we sure that `sched_setaffinity` actually locks the threads in MKL to the correct cores? Can we check somehow?

How do we mitigate the performance overhead of busy-waiting MPI threads? If they don't have work, we want them to just chill out and not do anything until necessary. Unfortunately this might mean extra overhead during execution on the device. We might want to tune the number of worker ranks so that they're almost always in use; maybe some function of the max graph diameter? That would be conservative.

We should run the master thread on the host, which sends start messages to the device. Might not be super fast, but the master thread will never be running application tasks.

# Power logging
Keep track of average power and max power (to make sure we don't blow power budget).
Use periodic interrupts to check (over a time window), log average and max.
Report at end (in log file).

With a power cap, it's possible that some resources go unused, since they provide little performance improvement or blow the power budget. A neat graph would show the when the runtime chooses to NOT use all the available resources, and when that decision is based on little performance improvement and when it's based on possibility of blowing the power budget.

# Setting the number of threads we can use
How do we pass the number of threads to use at runtime?
Use the environment variable `NUMTHREADS`

# Scheduling with power constraints
Ultimate goal: minimize makespan (ie, minimize `max_blevel`).
We'd prefer to give more cores to the task that gives higher performance per watt.

Use POWERLIMIT environment variable to set runtime limit. Defaults to 350 watts, which I think is higher than the TDP for our part, but we can raise this if necessary.

1. Enforce power limit. Only assign cores such that we predict we remain under the power limit.
    - Dumb version: remove cores so that we go under the power limit by limiting performance the least.
    - A particular configuration is unacceptable if it goes over the power limit.
        - Going over the limit means: the predicted power for this task is greater than what we have available.
    - If we're over the budget, remove a core from whichever task increases our `max_blevel` the least.
    - Keep removing a core until we get under the power limit.
    - Need to keep track of available power (kind of like cores) so that we don't go over. This gives us a budget for this scheduling period to work from.


For each task, figure out the next least powerful number of threads, the change in power, and the change in execution time. Then select whichever thread results in the smallest change in exec time. This is the thread that loses. Do something if this results in an increase in the number of threads. Change its thread assignment and subtract the power change from `sum_p`. Then we start over.

We currently have little implicit control over the core allocation method (currently, it's very ad-hoc, with little acknowledgement to hyperthreading). We might want to allow it to configure the types of cores allowed (ie, 2-way ht, 4-way ht, no ht).
