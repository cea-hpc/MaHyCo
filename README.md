# MaHyCo
Maquette hydrodynamique collaborative

# Pre-required

You should have installed *openmp*, *hwloc* and *kokkos*. 

## *openmp* and *hwloc*

To install *openmp* and *hwloc*, the recomanded way is to use your package manager.

For example on Debian/Ubuntu linux distributions :

```
sudo apt-get install libomp-dev hwloc libhwloc-dev
```

On Mac Os:

```
brew install hwloc
```

For Mac users, do not install *openmp*.

## *kokkos*

To install *kokkos*, first clone the repository:

```git clone https://github.com/kokkos/kokkos.git```

Then create a build directory and build *kokkos* in it :

```
mkdir build_kokkos
cd build_kokkos
cmake ../kokkos -DKokkos_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=/path/to/desired/kokkos/install/dir -DKokkos_ENABLE_OPENMP=On -DKokkos_ENABLE_HWLOC=On
make install
```

For Mac users, replace the `cmake` command with :

```
cmake ../kokkos -DKokkos_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=/path/to/desired/kokkos/install/dir -DKokkos_ENABLE_OPENMP=Off -DKokkos_ENABLE_HWLOC=On
```

# Download *MaHyCo*
Just clone the repository:

```
git clone https://github.com/hippo91/MaHyCo.git
```

# Compile *MaHyCo*
Create a directory in which *MaHyCo* will be built and change directory into it:

```
mkdir /tmp/build_mahyco
cd /tmp/build_mahyco
```

Then call *CMake* to generate the build files by specifying the path toward *kokkos* installation directory:

```
cmake -DCMAKE_BUILD_TYPE=Release -DKOKKOS_INSTALL_PATH=/path/to/kokkos/install/dir  /path/to/mahyco
```

And then just build it:

```
make -j 4
```

In the build directory the executable `eucclhydremap` is now available

## Add new integration test

To add a new integration test, in the `tests` directory, create a new directory which name will be referenced as `test_name`.

In this new directory adds :

- an `args.txt` file and fill it with the required arguments that have to be passed to `eucclhydremap` ;
- a `reference` directory which holds the expected results.

Then register the new test in the `CMakeLists.txt` file by adding it to the `TEST_NAMES` variable :

```
set( TEST_NAMES "SOD_X_bi_mat_Superbee_PenteBorne_mix"
                "SOD_X_bi_mat_MinMod_simple"
                "SOD_X_bi_mat_Superbee_PenteBorne_simple"
                "SOD_X_bi_mat_Superbee_simple"
                "SOD_Y_bi_mat_MinMod_simple"
                "SOD_Y_bi_mat_Superbee_PenteBorne_mixte"
                "SOD_Y_bi_mat_Superbee_PenteBorne_pure"
                "SOD_Y_bi_mat_VanLeer_PenteBorne_pure"
                "SOD_X_Mono_mat_MinModG_simple"
                "SOD_X_Mono_mat_MinMod_simple"
                "SOD_X_Mono_mat_SuperBee_simple"
                "SOD_X_Mono_mat_VanLeer_simple"
                "test_name" )
```

## Launch integration tests

To launch integration tests just run : 

```
ctest -j 4
```
 
in the directory where *MaHyCo* has been built.
