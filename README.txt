Pour utiliser Mahyco, recopier ce repertoire 'Mahyco'
chez vous.

Pour compiler Mahyco, se placer dans le repertoire correspondant
et lancer la commande

  cmake .

Pour tester Mahyco , se placer dans le repertoire src et lancer l'executable (extension .exe) en
specifiant en argument le jeu de donnees (extension .arc).

Par exemple :

  cd src/
  make
  ./Mahyco Mahyco.arc

Pour tester les évolutions avant rangement :

  cd src
  make
  ctest

Si tout est OK (aucune différence entre les passages, Tests tous marqués "Passed"), vous pouvez ranger...

Pour les sorties, elles sont dans le repertoire 'output'. Dans
ce repertoire, un repertoire 'courbes' contient les courbes
par iterations et le repertoire 'depouillement' le maillage et
les variables pour le post-traitement.

Vous pouvez ajouter les options suivantes pour chaque exemple. Les
options doivent etre ajoutees avant le jeu de donnees (qui doit
toujours etre le dernier argument).

 -arcane_opt max_iteration N  

     avec N le nombre d'iterations a effectuer

 -arcane_opt continue

     pour faire une reprise: continuer une execution precedente.


Pour lancer un cas en parallele, il faut specifier le service
de parallelisme via la variable d'environnement ARCANE_PARALLEL_SERVICE.
Les valeurs possibles sont: 'Mpi', 'Thread', 'MpiThread'.
Dans le cas ou on utilise des threads, il faut specifier leur nombre
via la variables d'environnement ARCANE_NB_THREAD.

Par exemple, pour 3 process MPI:

 export ARCANE_PARALLEL_SERVICE=Mpi
 mpiexec -n 3 ./Mahyco Mahyco.arc

pour 4 threads:

 export ARCANE_PARALLEL_SERVICE=Thread
 export ARCANE_NB_THREAD=4
 ./Mahyco Mahyco.arc

pour 3 process MPI et 4 threads (soit 12 sous-domaines)

 export ARCANE_PARALLEL_SERVICE=MpiThread
 export ARCANE_NB_THREAD=4
 mpiexec -n 3 ./Mahyco Mahyco.arc
