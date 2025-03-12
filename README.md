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

1. Use restart.py
```bash
cd build
# make sure the cwd is build.
../src/restart.py ...
```

Run it without argument for explanations for the variables

## FAQ

### Compilation Error

If you see error with missing `gsl` or `cblas`. Make sure you install GNU Scientific Library. It is `libgsl-dev` on Ubuntu. It includes `cblas` so you don't need to do anything to install `cblas`. But in other cases such as compiling and installing `gsl` from source, you might need to install `cblas` separately. 

Since the project uses `autoconf`, all the trick applies. Particularly I find [this trick with prefix](https://stackoverflow.com/questions/7561509/how-to-add-include-and-lib-paths-to-configure-make-cycle) useful when I have libraries installed in `~/.local`, which often happens when I am using a cluster that I don't have admin access to.
