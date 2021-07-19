Pour utiliser Mahyco, recopier ce repertoire 'Mahyco'
chez vous.

Pour compiler Mahyco, se placer dans le repertoire correspondant
et lancer la commande

  cmake . -B "repertoire build"

Pour tester Mahyco , se placer dans le repertoire "repertoire build"/src et lancer l'executable en
specifiant en argument le jeu de donnees (extension .arc).

Par exemple :

  cd "repertoire build"/src/
  make
  ./Mahyco Mahyco.arc

Pour tester les évolutions avant rangement :

  cd "repertoire build"/src
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


---------------------------
MEMO POUR LES ACCELERATEURS
---------------------------

Pour le support des accélérateurs dans Arcane, on suppose que l'on compile et installe Arcane dans un répertoire dédié "Arcane_ROOT".
Puis on compile Mahyco via cmake avec les commandes suivantes : 

rm -rf build
mkdir build
cd build
# Build with makefiles in parallel
cmake ..  -DWANT_CUDA=TRUE -DWANT_PROF_ACC=FALSE -DCMAKE_BUILD_TYPE=Debug -DArcane_ROOT=/your/path/to/installed/arcane -DCMAKE_CUDA_COMPILER=/chemin/vers/bin/nvcc
[Pour profiling avec nvtx : -DWANT_PROF_ACC=TRUE]
#cmake ..  
cmake --build . -- -j16
# Build with Ninja (natively parallel)
#cmake .. -G Ninja -DARCANE_WANT_CUDA=ON
#cmake --build . 
cmake --build . --target test

Execution sur une machine avec accélérateur :
 /chemin/vers/build/src/Mahyco -A,AcceleratorRuntime=cuda Donnees.arc

Execution avec instrumentation nsys [si projet configuré avec
-DWANT_PROF_ACC=TRUE, prise en compte des points d'entrée] :
 nsys profile --stats=true --force-overwrite true -o mahyco /chemin/vers/build/src/Mahyco -A,AcceleratorRuntime=cuda Donnees.arc

Exécution avec instrumentation nvprof (sortie ASCII dans mahyco.lognvprof) 
[si projet configuré avec -DWANT_PROF_ACC=TRUE, prise en compte des points d'entrée] :
 nvprof --print-api-trace --print-gpu-trace --log-file mahyco.lognvprof /chemin/vers/build/src/Mahyco -A,AcceleratorRuntime=cuda Donnees.arc

