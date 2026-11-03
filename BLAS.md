# Choosing a BLAS for HyperCom

HyperCom's influence matrix is the inverse of a dense L × L matrix (L = cells + genes). That
inverse is computed by LAPACK, and LAPACK's speed is the speed of the BLAS library R is linked to.
R itself ships with the **reference BLAS**: correct, single-threaded, no vector instructions, about
1 GFLOP/s. An optimised BLAS (OpenBLAS, Intel MKL, Apple Accelerate, BLIS) is 10–20× faster per
core and uses all cores.

Measured on the BU SCC (R 4.5.2, Xeon E5-2660 v3):

| | reference BLAS (NETLIB) | Intel MKL |
|---|---|---|
| `solve()` of a 1,500 × 1,500 matrix, 1 thread | 8.2 s | 0.50 s |
| HyperCom influence matrix, 5k cells | 410 s | 1–2 s |
| HyperCom influence matrix, 20k cells, 8 threads | ~6–7 h | 52 s |

The LR scoring does not use BLAS (it is HyperCom's own C++ kernel), so the BLAS only matters for
`generate_infMat()`. But with the reference BLAS that one step dominates the whole run.

## 1. Check which BLAS your R uses

In R:

```r
sessionInfo()$BLAS           # path of the BLAS library, or "FlexiBLAS <BACKEND>"
sessionInfo()$LAPACK
extSoftVersion()["BLAS"]     # the same, R >= 3.5
La_library()                 # the LAPACK library
```

How to read the answer:

| What you see | What it means | Fast? |
|---|---|---|
| `.../lib/R/lib/libRblas.so`, `.../Rblas.dll`, `libblas.so.3` with no vendor name | R's bundled reference BLAS, or the system reference BLAS | **no** |
| anything containing `openblas` | OpenBLAS | yes |
| anything containing `mkl` (e.g. `libmkl_rt.so`) | Intel MKL | yes |
| `.../libRblas.vecLib.dylib` or `Accelerate.framework` (macOS) | Apple Accelerate | yes |
| `FlexiBLAS NETLIB` | FlexiBLAS wrapper, currently pointed at the reference BLAS | **no** |
| `FlexiBLAS MKLOPENMP`, `FlexiBLAS OPENBLASOPENMP`, … | FlexiBLAS pointed at an optimised backend | yes |

If the name does not tell you, time it:

```r
m <- matrix(rnorm(2000^2), 2000)
system.time(solve(m))        # reference BLAS: ~10-20 s;  optimised, 1 thread: ~1 s;  optimised, many threads: well under 1 s
```

To see how many threads the BLAS is using, `RhpcBLASctl::blas_get_num_procs()` (install the
`RhpcBLASctl` package), or watch the CPU usage of the `solve()` above.

## 2. If you only have the reference BLAS: installing an optimised one

The BLAS is chosen when R is **linked**, so the fix depends on how R was installed. Pick your case.

**Linux, R from apt (Debian / Ubuntu).** R links the system `libblas.so.3`, which can be switched
without reinstalling R:

```bash
sudo apt install libopenblas-dev                       # or: intel-oneapi-mkl (Intel's apt repo) for MKL
sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu     # choose openblas (or mkl)
sudo update-alternatives --config liblapack.so.3-x86_64-linux-gnu
```

Restart R and check `sessionInfo()$BLAS`.

**Linux, R from dnf (Fedora / RHEL / Rocky).** R is linked to FlexiBLAS; install a backend and
select it (see section 3):

```bash
sudo dnf install openblas-openmp flexiblas
flexiblas default OPENBLAS-OPENMP
```

**Any platform, R from conda.** Pin the BLAS variant of the `libblas` package:

```bash
conda install -c conda-forge r-base "libblas=*=*openblas"    # or "libblas=*=*mkl"
```

**macOS, R from CRAN.** CRAN's R ships both the reference BLAS and Apple's Accelerate; switch by
pointing the BLAS symlink at the vecLib build (one command, reversible):

```bash
cd /Library/Frameworks/R.framework/Resources/lib
ln -sf libRblas.vecLib.dylib libRblas.dylib          # back to reference: ln -sf libRblas.0.dylib libRblas.dylib
```

**Windows.** R's `Rblas.dll` is the reference BLAS and there is no supported drop-in switch. The
practical routes are R from conda with MKL (above), or Linux under WSL.

**HPC clusters with modules.** Look for an R built against FlexiBLAS, OpenBLAS or MKL
(`module avail R`, `module avail openblas mkl flexiblas`); on the BU SCC, `R/4.5.2` is built against
FlexiBLAS, so nothing needs installing — only selecting (next section).

**Building R from source.** `./configure --with-blas=-lopenblas --with-lapack` (or `--with-blas="-lmkl_rt"`),
or `--with-blas=-lflexiblas` to be able to switch later.

## 3. Choosing which BLAS to use (FlexiBLAS)

Where R is linked to **FlexiBLAS** (the BU SCC, Fedora-family Linux, or a source build with
`-lflexiblas`), the backend is chosen at run time, per process, with no reinstall:

```bash
FLEXIBLAS=MklOpenMP Rscript my_script.R         # one run
export FLEXIBLAS=MklOpenMP                      # every R started from this shell / in this job script
flexiblas default MklOpenMP                     # a persistent user default, written to ~/.flexiblasrc
```

Backend names are case-insensitive (`MklOpenMP` is reported as `MKLOPENMP`). To see what is
available: `flexiblas list`, or list the configuration directory. On the BU SCC the `flexiblas`
command needs its own library on the path:

```bash
module load flexiblas/3.3.1
export LD_LIBRARY_PATH=$SCC_FLEXIBLAS_LIB:$LD_LIBRARY_PATH
flexiblas list
# or simply:  ls /share/pkg.8/flexiblas/3.3.1/install/etc/flexiblasrc.d/
```

> Do **not** run `flexiblas default` with no argument: it *removes* your user default.

On the SCC the backends are NETLIB (the default, reference), MklOpenMP, MklSerial,
OpenBLASOpenMP, BlisOpenMP, BlisSerial, AOCL, AOCL_Serial and ATLAS. **MklOpenMP is the one HyperCom was
benchmarked with**; OpenBLASOpenMP hung in R on this system and should be avoided there. Verify
inside R that the choice took effect:

```r
sessionInfo()$BLAS           # "FlexiBLAS MKLOPENMP"
```

A complete SCC job header for a HyperCom run:

```bash
module load gcc/12.2.0 R/4.5.2
export LD_PRELOAD="$(g++ -print-file-name=libstdc++.so.6)"   # the GLIBCXX fix used by this project's scripts
export FLEXIBLAS=MklOpenMP
export OMP_NUM_THREADS=1        # see "Threads" below
Rscript run_hypercom.R
```

Where R is **not** linked to FlexiBLAS, there is nothing to choose at run time: the BLAS is the one
R was linked to, and switching means section 2 (`update-alternatives`, the macOS symlink, or a
different conda build).

### Threads

An optimised BLAS runs the inverse on several cores. HyperCom sets the BLAS thread count itself
for the duration of `generate_infMat(..., threads = )` (through FlexiBLAS's, OpenBLAS's or MKL's
thread-setting function, restoring the previous value afterwards), so the usual environment
variables are not needed for it — and jobs can safely export `OMP_NUM_THREADS=1` for the sake of
forked workers, as this project's scripts do. If you call the inverse outside HyperCom, the
environment variables are `OMP_NUM_THREADS`, `MKL_NUM_THREADS` and `OPENBLAS_NUM_THREADS`.

### Troubleshooting

- `sessionInfo()$BLAS` still says `FlexiBLAS NETLIB` inside a batch job: the `FLEXIBLAS` variable
  was set in a shell that is not the job's (put the `export` in the job script itself).
- R hangs during `solve()` with an OpenMP backend after forking workers (`mclapply`): OpenMP and
  fork do not mix; set `OMP_NUM_THREADS=1` before forking, or use the serial backend in the workers.
- To see which BLAS R is actually linked against at the library level:
  `ldd $(R RHOME)/lib/libR.so | grep -i -E "blas|flexiblas|mkl"`.
- Results differ in the last digits between backends: expected. Optimised libraries sum in a
  different order; the HyperCom influence matrix agrees between NETLIB and MKL to ~1e-13 relative,
  and the LR priorities to ~1e-15, with identical rankings.
