# Lancement de la non-régression

Dans le répertoire de build :

```sh
cd build
cmake -- build .
ctest .
```

La configuration des tests à lancer est définie dans `src/CMakeLists.txt`.

Dans ce scénario, on définit une liste de cas à lancer en séquentiel `TEST_NAMES`, parallèle (4 et 8 procs `TEST_NAMES_PARA_4` et `TEST_NAMES_PARA_8`) avec protection-reprise `TEST_NAMES_PR`.

Pour chaque test des listes, on appelle le script ../NONREGRESSION/launch_test.sh qui se charge de lancer le cas en fonction des arguments qui sont donnés, et compare les résultats aux références versionnées.

Le lancement se fait via les fonctions 
+ `launch_computation`
+ `launch_computation_pr`
+ `launch_computation_para_4`
+ `launch_computation_para_8`

Puis comparaison des résultats avec la fonction compare_results.

`ctest` produit deux fichiers via l'appel du script launch_test.sh :
- `list_of_pb_exec` : contient la liste des cas qui ont eu des erreurs lors de l'exécution (erreur Arcane, ...)
- `list_of_cases_to_change` : contient la liste des cas qui changent les résultats (exécution ok mais comparaison ratée avec la référence)

Cette distinction permet d'éviter de supprimer accidentellement des références lorsque les cas n'ont pas tourné alors qu'ils figurent sur la liste des cas à mettre à jour.

Attention, cette liste ne contient que le nom des cas. Se référer au log de `ctest` pour connaître le contexte d'exécution des cas tests (seq, parallèle, ...).

Les cas de `list_of_cases_to_change` doivent ensuite être vérifiés (et mis à jour) avec le script `bascule_ref.sh`.

```sh
cd ..
./src/bascule_ref.sh .
```

## TODO

- [ ] Corriger les cas qui ne sont pas répétables d'une exécution à l'autre
- [ ] Vérifier les changements de résultats
- [ ] Gérer dans ces deux scripts les cas qui tournent dans plusieurs modes d'exécution (seq, para, pr)
- [ ] Paramétrer dans la ligne de commande lamise à jour auto des résultats, l'affichage ou non du diff et l'ouverture de paraview
