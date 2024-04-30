# Maquette Hydrodynamique Collaborative

## Compilation de Mahyco

Pour utiliser Mahyco, recopier ce repertoire 'Mahyco'
chez vous.

Pour compiler Mahyco, se placer dans le repertoire correspondant
et lancer la commande

```sh
  cmake . -B "repertoire build" -DArcane_ROOT="/your/path/to/arcane/"
```

Pour tester Mahyco , se placer dans le repertoire "repertoire build"/src et lancer l'executable en
specifiant en argument le jeu de donnees (extension .arc).

Par exemple :

```sh
  cd "repertoire build"/src/
  make
  ./Mahyco Mahyco.arc
```

## Non régression

Pour tester les évolutions avant rangement :

```sh
  cd "repertoire build"/src
  make
  ctest
```

Si tout est OK (aucune différence entre les passages, Tests tous marqués "Passed"), vous pouvez ranger...

## Sorties

Pour les sorties, elles sont dans le repertoire 'output'. Dans
ce repertoire, un repertoire 'courbes' contient les courbes
par iterations et le repertoire 'depouillement' le maillage et
les variables pour le post-traitement.


## Options de lancement d'un calcul

Vous pouvez ajouter les options suivantes pour chaque exemple. Les
options doivent etre ajoutees avant le jeu de donnees (qui doit
toujours etre le dernier argument).

```sh
 -arcane_opt max_iteration N  
```
     avec N le nombre d'iterations a effectuer

```sh
 -arcane_opt continue
```
     pour faire une reprise: continuer une execution precedente.


Pour lancer un cas en parallele, il faut specifier le service
de parallelisme via la variable d'environnement ARCANE_PARALLEL_SERVICE.
Les valeurs possibles sont: 'Mpi', 'Thread', 'MpiThread'.
Dans le cas ou on utilise des threads, il faut specifier leur nombre
via la variables d'environnement ARCANE_NB_THREAD.

Par exemple, pour 3 process MPI:

```sh
export ARCANE_PARALLEL_SERVICE=Mpi
mpiexec -n 3 ./Mahyco Mahyco.arc
```

pour 4 threads:

```sh
export ARCANE_PARALLEL_SERVICE=Thread
export ARCANE_NB_THREAD=4
./Mahyco Mahyco.arc
```

pour 3 process MPI et 4 threads (soit 12 sous-domaines)

```sh
export ARCANE_PARALLEL_SERVICE=MpiThread
export ARCANE_NB_THREAD=4
mpiexec -n 3 ./Mahyco Mahyco.arc
```