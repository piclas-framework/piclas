# Performance Analysis

This chapter describes different tools that can be utilized for measuring the computational performance of PICLas.

## Extrae and Paraver
Extra is a performance instrumentation tool that generates Paraver trace files that is distributed under LGPL-2.1 License and can be downloaded
from [https://github.com/bsc-performance-tools/extrae](https://github.com/bsc-performance-tools/extrae).
Paraver is a performance analysis GUI for visualizing the code tracing data:
[https://github.com/bsc-performance-tools/wxparaver](https://github.com/bsc-performance-tools/wxparaver)
They are part of the BSC tool set that is found under [https://tools.bsc.es/downloads](https://tools.bsc.es/downloads).

```{note}
`wxparavwer` can simply be downloaded as pre-compiled binary file as it is only used for viewing the results.
```

### Installation
See the `README` and `INSTALL` file in the git repository of the package.

(sec:code-instrumentation)=
### Code Instrumentation

```{note}
Tested with extrae version 3.8.3
```

In PICLas, the extrae code instrumentation for the very basic modules is already implemented, see the in-code statements, e.g.,

    #ifdef EXTRAE
    CALL extrae_eventandcounters(int(9000001), int8(1))
    #endif

    ! Initialization

    #ifdef EXTRAE
    CALL extrae_eventandcounters(int(9000001), int8(0))
    #endif

which spans the complete initialization phase. Other regions are the field and particle modules (pure DG, PIC, DSMC, etc.) and the
instrumentation is activated by setting the PICLas compile flag

    PICLAS_EXTRAE = ON

in the cmake settings.

Examples that are already instrumented are


```{table} Examples of instrumented code blocks
---
name: tab:ExtraeInstrumentation
---
|                 Function                | Intger Value |                            Source                            |
|             :--------------:            | :----------: |                 :--------------------------:                 |
|              Initialization             |       1      |                     `./src/piclaslib.f90`                    |
|              Load Balancing             |       2      |              `./src/loadbalance/loadbalance.f90`             |
|         Write State file to .h5         |       3      |                `./src/io_hdf5/hdf5_output.f90`               |
|    Field Solver (HDG with CG solver)    |       4      |                      `./src/hdg/hdg.f90`                     |
|             Particle Solver             |       5      | `./src/timedisc/timedisc_TimeStepPoissonByBorisLeapfrog.f90` |
|              Particle Push              |       5      |          `./src/timedisc/timedisc_TimeStep_BGK.f90`          |
|           `PerformTracking()`           |      50      |       `./src/particles/tracking/particle_tracking.f90`       |
|     `CALL UpdateNextFreePosition()`     |      51      |          `./src/timedisc/timedisc_TimeStep_BGK.f90`          |
| `CALL BGK_DSMC_main()`  or `BGK_main()` |      52      |          `./src/timedisc/timedisc_TimeStep_BGK.f90`          |
|                 Analysis                |       6      |                  `./src/analyze/analyze.f90`                 |
```

### Tracing the code

#### Load the required Modules

On the target system, the extrae software packages must be installed and loaded via, e.g.,

    module load extrae

#### Create *tracing.sh* and *extrae.xml* in the simulation directory

Create a shell script *tracing.sh* (must be executable) with the following content

    #!/bin/bash

    export EXTRAE_CONFIG_FILE=/path/to/extrae.xml
    export LD_PRELOAD=${EXTRAE_HOME}/lib/libmpitracef.so

    $*

where the path to the directory containing the *extrae.xml* file must be inserted.

```{note}
`LD_PRELOAD` might only required when no user-defined instrumentation is used. If `PICLAS_EXTRAE=ON` is used during
compilation, the line with `LD_PRELOAD` can be commented out or removed.
```

Furthermore, a configuration file *extrae.xml* is required that defines which hardware counters, which should be traced

    <?xml version='1.0'?>

    <trace enabled="yes"
     home="/opt/hlrs/non-spack/performance/extrae/3.7.1-mpt-2.23-gcc-9.2.0"
     initial-mode="detail"
     type="paraver"
    >

      <openmp enabled="no" ompt="no">
        <locks enabled="no" />
        <taskloop enabled="no" />
        <counters enabled="yes" />
      </openmp>

      <pthread enabled="no">
        <locks enabled="no" />
        <counters enabled="yes" />
      </pthread>

      <counters enabled="yes">
        <cpu enabled="yes" starting-set-distribution="1">
          <set enabled="yes" domain="all" changeat-time="0">
            PAPI_TOT_INS,PAPI_TOT_CYC
          </set>
          <set enabled="no" domain="all" changeat-time="0">
            PAPI_TOT_INS,PAPI_TOT_CYC,PAPI_VEC_SP,PAPI_SR_INS,PAPI_LD_INS,PAPI_FP_INS
            <sampling enabled="no" period="1000000000">PAPI_TOT_CYC</sampling>
          </set>
        </cpu>

        <network enabled="no" />

        <resource-usage enabled="no" />

        <memory-usage enabled="no" />
      </counters>

      <storage enabled="no">
        <trace-prefix enabled="yes">TRACE</trace-prefix>
        <size enabled="no">5</size>
        <temporal-directory enabled="yes">/scratch</temporal-directory>
        <final-directory enabled="yes">/gpfs/scratch/bsc41/bsc41273</final-directory>
      </storage>

      <buffer enabled="yes">
        <size enabled="yes">5000000</size>
        <circular enabled="no" />
      </buffer>

      <trace-control enabled="yes">
        <file enabled="no" frequency="5M">/gpfs/scratch/bsc41/bsc41273/control</file>
        <global-ops enabled="no"></global-ops>
      </trace-control>

      <others enabled="yes">
        <minimum-time enabled="no">10M</minimum-time>
        <finalize-on-signal enabled="yes"
          SIGUSR1="no" SIGUSR2="no" SIGINT="yes"
          SIGQUIT="yes" SIGTERM="yes" SIGXCPU="yes"
          SIGFPE="yes" SIGSEGV="yes" SIGABRT="yes"
        />
        <flush-sampling-buffer-at-instrumentation-point enabled="yes" />
      </others>

      <sampling enabled="no" type="virtual" period="50m" variability="10m" />

      <dynamic-memory enabled="no" />

      <input-output enabled="no" />

      <syscall enabled="no" />

      <merge enabled="no"
        synchronization="default"
        tree-fan-out="16"
        max-memory="512"
        joint-states="yes"
        keep-mpits="yes"
        sort-addresses="yes"
        overwrite="yes"
      />

    </trace>

Here, the MPI library with PAPI_TOT_INS and PAPI_TOT_CYC counters are traced. Note that the path to the extrae directory is defined
under

    home="/opt/hlrs/non-spack/performance/extrae/3.7.1-mpt-2.23-gcc-9.2.0"

### User functions

```{warning}
This section is experimental!
```

Compile the code with

    -finstrument-functions

and supply the names of the functions that are exclusively traced in a file `user-functions.dat` containing the hash and name of
each function in a comma-separated list, e.g.,

    000000000042d2e0#__mod_timedisc_MOD_timedisc

where the hash is acquired via

    nm -a bin/piclas_extrae | grep -in timedisc

or from the lib via `nm -a lib/libpiclas.so` if the function is in the shared library.
To use the `user-functions.dat` file in extrae, add the following block to the `extrae.xml` file.

    <user-functions enabled="yes" list="/absolute/path/to/user-functions.dat" exclude-automatic-functions="no">
      <counters enabled="yes" />
    </user-functions>

where the absolute path to `user-functions.dat` is supplied.

#### Run the application

Run the application and convert the output to Paraver format

```{note}
The extrae instrumented executable has a different name, which ends on `_extrae`
```

Execute `mpirun` and pass the `tracing.sh` script

    mpirun -np 32 tracing.sh piclas_extrae parameter.ini

The following command can be appended to the submit script directly after `mpirun`.

#### Convert the Extrae output for Paraver

The tracing output stored in *TRACE.mpits* is then converted to a Paraver file via

    ${EXTRAE_HOME}/bin/mpi2prv -f TRACE.mpits -o tracing.prv

e.g.,

    /opt/hlrs/non-spack/performance/extrae/3.7.1-mpt-2.23-gcc-9.2.0/bin/mpi2prv -f TRACE.mpits -o pilcas.32ranks.prv

or using mpirun

    mpirun -np 64  ${EXTRAE_HOME}/bin/mpimpi2prv -f TRACE.mpits -o tracing.prv

which will create a file containing the tracing events (.prv), list of registered events (.pcf) and cluster topology description (.row).

#### Analysing the results with Paraver

```{note}
Tested with wxparaver version 4.9.2
```

Open Paraver

    wxparaver

and load a trace file for Paraver

- Open the `.prv` file via *File* -> *Load Trace* and the possible quantities are already shown under *Workspaces*, e.g.,
  *Useful+MPI+PAPI ...*.
- To view one of these properties, go to *Hints* -> *Useful* -> *Useful Duration*, which opens a separate window displaying the data.
- It shows the MPI ranks vs. the wall time and shows the calculation time for each trace, i.e., how much of the wall time was
  actually spent for calculation (the useful part of the simulation).
- On the bottom left go to *Files \& Window Properties* and select *Window Properties*.
- Under *Properties Mode*, change the value from *Basic* to *Full* and select the drop down box
  - Have a look at the field *values* -> *...* -> `9000001` (piclas directives instrumented by hand) to see if they have been
  correctly used
  - *Filter* -> *Events* -> *Event type* and set *Function* to `=`
  - *Filter* -> *Events* -> *Event type* and set *Types* to `9000001`
  - *Filter* -> *Events* -> *Event value* and set *Function* to `=`
  - *Filter* -> *Events* -> *Event value* and set *Values* to `1`

  here, the actual tracing event number has to be used as defined in {ref}`sec:code-instrumentation`, e.g., `1` as for `int8(1)`.
  For a list of pre-defined settings, see {ref}`tab:ExtraeInstrumentation`

- Right-click into window *Useful Duration @ *.prv* -> *View* -> *Event Flags* to activate the user-instrumented events from {ref}`sec:code-instrumentation`

To synchronize the views between different windows, e.g., *Useful Duration* and *MPI call* or simply two windows *Useful Duration*
that each display a different *Event value* to show where a function instrumentation starts and ends

- Right-click into window *Useful Duration @ *.prv* -> *Synchronize* -> *[ ] 1* (select a group)
- Right-click into window *MPI call @ *.prv* -> *Synchronize* -> *[ ] 1* (select a group)




## Intel® VTune™
Intel® VTune™ is a performance analysis tool with GUI for applications running on x86 systems for Linux and Windows developed by Intel®.

### VTune Installation
Download the Intel VTune Profiler Source files for Linux and extract the files:

 - [https://software.intel.com/en-us/vtune/choose-download#standalone](https://software.intel.com/en-us/vtune/choose-download#standalone)
 - [https://software.intel.com/content/www/us/en/develop/articles/oneapi-standalone-components.html#vtune](https://software.intel.com/content/www/us/en/develop/articles/oneapi-standalone-components.html#vtune)

 A user guide can be found here:

  - [https://software.intel.com/content/www/us/en/develop/documentation/vtune-help/top/launch/getting-started.html](https://software.intel.com/content/www/us/en/develop/documentation/vtune-help/top/launch/getting-started.html)

Install VTune via the command line script (or the GUI installer)

    sudo ./install.sh

The installed environment is meant to run in a bash shell. The GUI can be started by

    bash
    source /opt/intel/vtune_profiler_2020.0.0.605129/vtune-vars.sh
    vtune-gui

Compile PICLas with "-g" to produce debugging information in the operating system's native format, which is required for detailed analysis performed by VTune.

### Batch jobs
VTune can also be run in batch mode without the GUI. For a list of available options, see

    vtune -help

and

    vtune -help collect

To run a simulation with 10 MPI threads and collect data, simply run

    mpirun -np 10 vtune -collect memory-consumption -trace-mpi -result-dir ~/intel/vtune/projects/feature_branch/r0001mc ~/piclas/particle.dev/build.vtune/bin/piclas parameter.ini

and specify where the output data of vtune should be collected.

### Usage
In the Vtune GUI, set the path to the executable, the parameters (parameter.ini DSMC.ini) and the working directory (where the executable is run).

Hit the "Start" button and wait. The piclas std-out is dumped directly into the shell where vtune-gui was launched. The output can
be redirected to a shell that is displayed in VTune by Setting: Options → General → "Product output window"

## Valgrind
Valgrind is a complete suite of tools for debugging/profiling licenced under GPL. The complete documentation can be found [here](https://www.valgrind.org/docs/manual/index.html).

### Installation of Valgrind
Valgrind is provided through the repository of all major Linux distributions. Install it with the package manager of your choice.

### Execution of Valgrind
Valgrind is composed of individual tools, each tailered to debug or profil a specific aspect. All tools need PICLas compiled with
"-g" to produce debugging information in the operating system's native format.

#### Callgrind
Callgrind tracks the call graph and duration for each function.

    valgrind --tool=callgrind ./piclas parameter.ini

The output file can be opened with KCacheGrind or converted using gprof2dot. The options `-n PERCENTAGE, --node-thres=PERCENTAGE / -e PERCENTAGE, --edge-thres=PERCENTAGE` eliminate nodes/edges below this threshold [default: 0.5].

    gprof2dot -n0 -e0 ./callgrind.out.1992 -f callgrind > out.dot
    dot -Tpng out.dot -o out.png

In both cases, make sure you have GraphViz installed.

#### Memcheck
Memcheck keeps track of every allocated memory and shows memory leaks or invalid accesses to memory/pointers. Run it with

    valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -s --suppressions=/share/openmpi/openmpi-valgrind.supp ./piclas parameter.ini

OpenMPI handles its memory independently, so memcheck will always report memory leaks due to OpenMPI. Using the provided
suppressions hides these falls flags. Combining memcheck with the GCC sanitize flag should provide full memory coverage.

#### Massif
Massif keeps track of the current memory usage as well as the overall heap usage. It helps finding code segments that hold on to
memory after they should. Run it with

    valgrind --tool=massif ./piclas parameter.ini

The output file can be opened with massif-visualizer.

#### DHAT
DHAT tracks memory allocations and inspects every memory access to a block. It is exceptionally useful to find code segments
where memory is allocated unused or rapidly re-allocated. Run it with

    valgrind --tool=dhat ./piclas parameter.ini

Afterwards, open `//usr/lib/valgrind/dh_view.html` in a web browser and load the output file.












