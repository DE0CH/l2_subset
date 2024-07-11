## Instruction for running the code

1. Compile the source code
```bash
autorecof --install
mkdir build
cd build
../configure
make
cd ..
```

The run scripts are separated into tro parts, one with pre-computation some data needed for all trials. When compiling with flag `-DCOMPUTE_MODE=1`, it computes the matrix (documentation needed what that actually is) and store it to disk. When compiling with flag `-DCOMPUTE_MODE=2`, it doesn't do much, but there are some ideas such as pre-computing the weights (documentation needed) so speed up starting time. 

The second part is the running script for a trial, which has a random start, and it loads data from the the pre-computed data. More specifically, it actually mmap the data to conserve memory and improve startup time.

1. Set variables
```bash
export N=1000000
export M=50
export D=3
```

1. Generate the points
```bash
build/gen_points "$D" points.txt "$N"
```

1. Pre-compute. This only needs to be run once for each N, M and D.
```bash
build/l2_subset_compile_matrix points.txt $M points.p
```

1. Run the trial, the trial number also serves at the seed for the code
```bash
trial=1
python3 build/l2_subset_from_compiled_matrix points.p $trial
```

## FAQ

### Compilation Error

If you see error with missing `gsl` or `cblas`. Make sure you install GNU Scientific Library. It is `libgsl-dev` on Ubuntu. It includes `cblas` so you don't need to do anything to install `cblas`. But in other cases such as compiling and installing `gsl` from source, you might need to install `cblas` separately. 

Since the project uses `autoconf`, all the trick applies. Particularly I find [this trick with prefix](https://stackoverflow.com/questions/7561509/how-to-add-include-and-lib-paths-to-configure-make-cycle) useful when I have libraries installed in `~/.local`, which often happens when I am using a cluster that I don't have admin access to.
