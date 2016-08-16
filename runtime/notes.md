# Todo: Newest near top
1. Write.

What's the goal of this paper? Are we trying to show that we can schedule tasks better because we use the smart makespan-minimizing algorithm? Or are we trying to demonstrate that by using Pareto frontiers to represent the power--performance relationships of our kernels that we get better performance per watt? That sounds like a good story to me. Could the argument be made that we don't in fact use the state of the art scheduling heuristics that minimize makespan? I mean, we could go that far, but it would require figuring out the scheduling questions behind non-monotonic relationship between number of threads and performance/power as well as non-unit increases in power per point on Pareto frontier. If we can tackle that non-unit increase business, we could probably talk about thread over prescription too as a side effect.

Is this all too much to handle in the next two weeks? If we punt, that would be the absolute easiest. Might as well run those experiments to see how well they do.

Note: Pareto curves are not always monotonic in the number of threads. This makes thread-based constraining by walking down the curves impossible.

We should probably modify the advanced scheduler not to make decisions based on local execution time minimization; as indicated in the prior work paper, a list scheduler won't minimize the per-task execution time, but rather the estimated makespan of the entire graph.

Do we just perform list scheduling based on power and log the instances where we oversubscribe the cores? If that's the case, then our performance should go through the floor, but at least we'd know why. Or if that's the case, then we perform some secondary form of adjusting to fit within the number of threads available.

The goal is to MINIMIZE the MAXIMUM makespan, since that's what presumably limits the performance of the entire graph.

1. Come up with initial configuration (ie, best possible perf on frontier)
2. Estimate the makespans for each task, given the current configuration
3. Find the one with the shortest estimated makespan and walk down to next configuration
4. 

What's the difference between running the advanced scheduler with a power limit and without, but limiting the threads? I think we may have to consider using the normal baseline scheduler.


Ok no more MPI garbage. Steps:
1. Write simple sequential scheduler. While there are tasks, this scheduler makes a thread to run the task then sits on `pthread_join` until it's done. Then it goes to the next one.
2. Write a more complicated scheduler, launching as many tasks as possible and checking for completion using nonblocking joins.




Make a bunch of worker threads on the host. If a task needs to be launched, have the thread launch the task synchronously. It can then report back using all the appropriate synchronization primitives (ie, wait, checkdone, etc).

*NEED TO MEMALIGN TO 64*




Use `pragma offload` for everything. A task gets offloaded to the device. The setup phase, that allocates the heap, should have to setup the initial data. The `allocateHeap` function would be called after all buffers have been initialized with the problem data, and would marshal it all and send it across.

Could we run with just two MPI ranks, one on the host, one on the device. Then we launch tasks dynamically from the manager task on the device. How is this different than running everything natively?
Alternatively, we could launch EVERYTHING offload from the host. That requires synchronizing the data transfer.

Why not just run everything on the device? Including the scheduler.

What happens when we want to schedule a task? Do we do better to fire right away or wait until more threads are available? Whenever a task finishes and frees its threads, all pending tasks reevaluate the best configuration for the given number of available threads. Could there exist times when waiting (and not scheduling) would be better than to scheduler immediately and probably perform worse?

The baseline scheduler (and the advanced one, by extension), messes up the data. However, the profiling one and the sequential one seem to work out ok. Also it's only for some configurations its borked. I don't get it.

*For `spmv`, there is not a direct relationship between the number of non-zeros (or the histogram of non-zeros) to the performance of the kernel.*

The baseline scheduler should perform about as poorly as the power-aware one if the timing models are totally off; neither can accurately guess well. The power-aware one may exhibit different characteristics, assuming that power can be modeled more accurately than time.

*Can power be more accurately modeled than time?* Are these kernels consistent in their power usage, even if the total time is different (both avg power and max power?)?

What do we care about the distinction between avg and max power? We can only sample at 50ms granularity, so even if our peak occurs at a faster rate, who cares? It just means we won't be able to measure any events that happen faster than that.

```
    array[buffers]:
        each buf points to a location at an offset within the actual thing.
        we just need to know the offset, send that, and have it create a new buff from offset()
        do we have to do that for everyone? ie each thread? or does the voodoo happen with shmem? I think we have to send to everyone
```
Why am I doing this? To make experiments faster?

For `ddot`, even though the hand-tuned version runs an order of magnitude slower, it doesn't have as high a relative degree of variability in performance. I still think this is due to the tasks running too quickly.

So the tasks don't run with different sizes, they just run with different behaviors. This may be due to having data already loaded in the cache, or by another process. Meaning, most of the time the first iteration of a task is way slower than everything after it. NOPE: it's because `SPMV` has a variable amount of work depending on the sparseness, which is data dependent. The second time around, the matrix is different, which means that the total number of flops changes. This results in fewer total operations.

_THE GRANULARITY IS TOO SMALL FOR THESE TASKS._ It's impossible to measure individual task execution, and the variance is too high (for example, the empty start task has runtime proportional to the max number of threads, even though nothing even happens).

Maybe we can take the frontiers we have and run with it. We'll just need to come up with power numbers, which might require using our existing frontiers or by running each task a bunch of times to get avg power numbers. Although that's not as good as peak power, but I guess if we run the tasks over and over enough, that should be fine.

I don't get what's going on. I need to have tasks that run in a reasonable amount of time (~seconds), but not slower than that. Which either means massive problem sizes, or many iterations of the same kernel over and over. How would that work with the dag? I wouldn't really be able to rerun the same task over and over, that doesn't make a ton of sense (unless we do it the same number of iterations for each task? i.e., artificially increase the number of iterations by some shared ratio to make them run longer).

It's also the case that these tasks have execution time and power consumption proportional to the input size, which is not constant across executions (for example, multiple `gemm` calls in `cholesky` with different input sizes). Do we need to make a unique frontier for each input size we run into? Do we need to come up with some other method for modeling the power/performance relationships that comprise the frontiers? This should kind of be my work already I guess. 

Most of these cg tasks are still pretty short-running (milliseconds). Perhaps we should make the problem size much bigger?
Still fucked. Giving more threads to SPMV makes it slower. What?

DONE: We need to figure out a good baseline. The `BasicScheduler` currently only issues a single thread per task, which might not make a lot of sense. Perhaps we need a modified version of the `AdvancedScheduler` that does everything except for the power-aware modifications (ie, it only uses the existing thread-aware scheduler).

Also need to change how the pareto frontier is passed. Currently nothing is ordered, and we add a config for each number of threads. Perhaps this should be somehow ordered differently (ie, sorted).

DONE: Need to fix the scheduler so that if it can't find a task (due to power budget), we drop a task entirely.

## Task power starting point
When we schedule a task, we need to know how much power the addition of a single task will cost the runtime. That means that we consider the gradient along the pareto frontier, which gives us the configuration for the quantized decrease in power.

However, these frontiers have been constructed in isolation; ie, each task is run by itself. The delta is important, but so is the absolute increase over not running the task at all (ie having zero threads running). This feeds into the energy scaling work, where we consider how much power is consumed by the device being ON, but idle.

We should run the extrapolated idle power experiments to determine this baseline, so that we know how much _additional_ power is required to run the task, on top of an idle system.

In addition, we need to take into consideration that baseline power consumption is always present. Perhaps we can subtract this value from the power cap, indicating the amount of power headroom we have to work with, assuming a powered on and running system? It definitely detracts from the race to idle camp, but from what I've heard, people say that they run at full power mode all the time anyway since they never know when things will be needed. This is true especially when you consider that, as it stands now, you can't really turn cores off. If this capability were possible, we'd have to work that into the model (ie, not just the dynamic power cost, but the incremental power cost to take the core out of a low power state).

The each thread may be able to do more work (sqrt may not use all the core functional units) but since it's consistent, it should scale correctly (ie, the y-intercept should be correct).

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
