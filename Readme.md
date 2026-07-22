# MaHyCo

## Pré-requis

- **CMake >= 3.21**
- **C++17**
- **Arcane >= 3.16.0** (source + installation séparées)
- Compilateur C++ et, si GPU, `nvcc`

Sur le cluster Inti (CEA), l'environnement est chargé automatiquement par les scripts via `env_gcc12.3_cuda12.4_mpi4.1.7.sh` (modules cmake, gcc/12, cuda/12.4, openmpi/4.1.7, ...).

## Compilation

Le projet fournit deux scripts de build dans `bin/` : le premier compile et installe **Arcane**, le second compile **MaHyCo** contre cette installation.

### 1. Compiler Arcane

```bash
bin/do_build_Arcane.sh <PATH_TO_ARCANE_SRC> <PATH_TO_ARCANE_INSTALL_ROOT> [options]
```

Arguments et options principales :
- `<PATH_TO_ARCANE_SRC>` : chemin vers les sources Arcane
- `<PATH_TO_ARCANE_INSTALL_ROOT>` : chemin racine de l'installation. L'arborescence finale sera :
  `<ROOT>/arcane<VERSION><_SUFFIX>/${OS}/${BUILD_TYPE}/`
- `-build=release|debug|check` (défaut : `release`)
- `-acc=CUDA` (ou `ROCM`, `HIP`, `SYCL`) pour activer le support accélérateur
- `-suff=<hash_ou_label>` : suffixe ajouté au nom de version dans le chemin d'installation
- `-arc_tests` : lance les tests Arcane après compilation (long)
- `-v` : mode verbeux

**Exemple (cluster Inti, release avec CUDA)** :
```bash
bin/do_build_Arcane.sh /path/to/arcane/src /path/to/arcane/install -build=release -acc=CUDA
```

**Exemple (laptop, debug CPU)** :
```bash
bin/do_build_Arcane.sh /path/to/arcane/src /path/to/arcane/install -build=debug
```

### 2. Compiler MaHyCo

```bash
bin/do_build_Mahyco.sh <PATH_TO_MAHYCO_SRC> <PATH_TO_ARCANE_INSTALL> [options]
```

Arguments et options principales :
- `<PATH_TO_MAHYCO_SRC>` : chemin absolu vers les sources MaHyCo
- `<PATH_TO_ARCANE_INSTALL>` : chemin absolu vers l'installation **Arcane** (même racine que ci-dessus)
- `-build=release|debug|check` (défaut : `release`)
- `-acc=CUDA` : active le support GPU (nécessite qu'Arcane ait été compilé avec `-acc=CUDA`)
- `-nvtx=on|off` : active/désactive le profiling NVTX (défaut : `on`)
- `-v` : mode verbeux

Le script crée automatiquement le répertoire :
```
build_${OS}[_CUDA]/${BUILD_TYPE}/
```
et y lance `cmake` puis `make`.

**Exemple (cluster Inti, release avec CUDA)** :
```bash
bin/do_build_Mahyco.sh $(pwd) /path/to/arcane/install -build=release -acc=CUDA
```

**Exemple (laptop, debug CPU)** :
```bash
bin/do_build_Mahyco.sh $(pwd) /path/to/arcane/install -build=debug
```

## Exécution

L'exécutable se trouve dans le répertoire de build :

```bash
# Séquentiel
./build_.../src/Mahyco Donnees.arc

# MPI
mpiexec -n 4 ./build_.../src/Mahyco Donnees.arc

# Multi-thread avec l'API Accélérateur (CPU, sans MPI)
./build_.../src/Mahyco -A,T=4 Donnees.arc

# GPU (1 processus, 1 GPU)
./build_.../src/Mahyco -A,AcceleratorRuntime=cuda Donnees.arc

# Multi-GPU (1 sous-domaine MPI par GPU, 4 GPUs)
mpiexec -n 4 bin/wrapper_mgpu.bash ./build_.../src/Mahyco -A,AcceleratorRuntime=cuda Donnees.arc
```

Options à placer **avant** le fichier `.arc` :
- `-A,MaxIteration=$N` ou `-arcane_opt max_iteration <N>` : limite à N itérations
- `-arcane_opt continue` : reprise depuis une protection (checkpoint)

## Tests

Les tests sont des tests de **non-régression**. Ils comparent le contenu de `output/depouillement` avec des résultats de référence stockés dans chaque dossier de `NONREGRESSION/`.

Lancer tous les tests depuis le répertoire de build :
```bash
cd build_.../
ctest
```

Lancer une sous-catégorie :
```bash
# Tests séquentiels
ctest -R seq_

# Tests parallèles sur 4 cœurs
ctest -R para_4_

# Tests GPU (uniquement si compilé avec -acc=CUDA)
ctest -R cuda_
```

Si des différences sont attendues et legitimes, mettre à jour les références :
```bash
# Le fichier list_of_cases_to_change est généré par ctest en cas d'échec
./NONREGRESSION/bascule_ref.sh . build_<os>_<arch>/<BuildType>
```

Variables d'environnement pour `bascule_ref.sh` :
- `AFFICHE_DIFF=1` : affiche le diff
- `OUVRE_PARAVIEW=1` : ouvre ParaView pour comparer
- `BASCULE_FORCEE=1` : force la mise à jour sans confirmation

## Sorties

Les résultats d'exécution sont écrits dans le répertoire `output/` :
- `output/courbes/` : courbes par itération
- `output/depouillement/` : maillage et variables pour le post-traitement (format Ensight)

## Notes

- Sur le cluster Inti, le launcher MPI par défaut est `/usr/bin/ccc_mprun`.
- Si vous compilez MaHyCo avec `-acc=CUDA`, le script vérifie qu'Arcane fournit bien `libarcane_accelerator_cuda_runtime.so`.
- Sur le cluster, les tests non-régression s'exécutent dans des répertoires temporaires sous `${CCCSCRATCHDIR}/MAHYCO/NONREGRESSION/`.
- **Multi-GPU** : le script `bin/wrapper_mgpu.bash` bind chaque rang MPI à un GPU local via `CUDA_VISIBLE_DEVICES=$LOCAL_RANK`. Il détecte automatiquement `OMPI_COMM_WORLD_LOCAL_RANK` (OpenMPI) ou `SLURM_LOCALID` (SLURM).
- **Multi-thread CPU** : l'option `-A,T=${N}` exploite l'API Accélérateur d'Arcane sur CPU avec `N` threads, sans recourir à MPI.
- **MPI** : il n'est pas nécessaire de définir `ARCANE_PARALLEL_SERVICE=Mpi` ; Arcane détecte automatiquement le contexte MPI quand on lance via `mpiexec`/`ccc_mprun`.
