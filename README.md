## Installation

Pour utiliser Mahyco, recopier ce répertoire 'Mahyco'
chez vous.

Pour compiler Mahyco, se placer dans le répertoire correspondant
et lancer la commande
```bash
    cmake . -B "repertoire build" -DArcane_ROOT="/your/path/to/arcane/"
```

Voir le [guide d'installation complet](https://github.com/cea-hpc/MaHyCo/wiki/Installation).

## Lancer MahyCo

Pour tester Mahyco , se placer dans le répertoire `repertoire_build/src` et lancer l'exécutable en
spécifiant en argument le jeu de données (extension `.arc`).

Par exemple :
```bash
    cd "repertoire build"/src/
    make
    ./Mahyco Mahyco.arc
```

Voir le [guide d'utilisation complet](https://github.com/cea-hpc/MaHyCo/wiki/Utilisation).

## Tests
Pour tester les évolutions avant rangement :
```bash
    cd "repertoire build"/src
    make
    ctest
```
Si tout est OK (aucune différence entre les passages, Tests tous marqués "Passed"), vous pouvez ranger...

## Visualisations et sorties

Pour les sorties, elles sont dans le repertoire `output`. Dans
ce répertoire, un répertoire `courbes` contient les courbes
par itérations et le répertoire `depouillement` le maillage et
les variables pour le post-traitement.

## Options en ligne de commande

Vous pouvez ajouter les options suivantes pour chaque exemple. Les
options doivent être ajoutées avant le jeu de données (qui doit
toujours être le dernier argument).
```bash
    -arcane_opt max_iteration N  
```
avec N le nombre d'iterations a effectuer
```bash
    -arcane_opt continue
```
pour faire une reprise: continuer une exécution précédente.

## Lancer MaHyCo en parallèle

Pour lancer un cas en parallèle, il faut specifier le service
de parallélisme via la variable d'environnement `ARCANE_PARALLEL_SERVICE`.
Les valeurs possibles sont: `Mpi`, `Thread`, `MpiThread`.
Dans le cas ou on utilise des threads, il faut spécifier leur nombre
via la variables d'environnement `ARCANE_NB_THREAD`.

Par exemple, pour 3 process MPI:
```bash
    export ARCANE_PARALLEL_SERVICE=Mpi
    mpiexec -n 3 ./Mahyco Mahyco.arc
```
pour 4 threads:
```bash
    export ARCANE_PARALLEL_SERVICE=Thread
    export ARCANE_NB_THREAD=4
    ./Mahyco Mahyco.arc
```
pour 3 process MPI et 4 threads (soit 12 sous-domaines)
```bash
    export ARCANE_PARALLEL_SERVICE=MpiThread
    export ARCANE_NB_THREAD=4
    mpiexec -n 3 ./Mahyco Mahyco.arc
```
