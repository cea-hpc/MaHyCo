# Etat des lieux sur l'IPG dans Mahyco -- 25/04/2024

(fichier versionné dans la branche ipg de Mahyco)

Etat des lieux des travaux IPG dans Mahyco. 
Pour certaines thématiques, des idées d'implantation sont suggérées.

## Rappel de l'existant

1 module Ipg qui :
- crée 3 particules avec une vitesse constante orientée selon 2 * k * Pi / 3, avec k = localId de la particule
- calcule les positions à chaque pas de temps des particules, par incrément v * delta t
- déclenche l'écriture de sorties via l'appel à une interface

1 package ipg_output
- 1 interface implémentée par un service qui écrit les sorties au format vtk (fichier ascii)

Pour jouer, j'ai créé 5 variables particules :
- "ParticleCoordinates"
- "ParticleVelocity"
- "ParticleRadius"
- "ParticleWeight" 
- "ParticleTemperature"
qui sont initialisées n'importe comment. Les valeurs à l'intérieur ne veulent rien dire, c'est pour vérifier que j'arrive à les sortir.

1 cas test avec des particules

## Sorties des particules

Actuellement, écriture d'un fichier ascii (vtk) par temps de sortie. 
En séquentiel uniquement.
Ca dépanne mais pas sûre que ça soit efficace pour beaucoup de particules...

### Piste d'amélioration

- Créer "un film" qui chargerait toutes les particules.
- Coder le rythme de sorties en temps = synchroniser les sorties Particules avec les sorties Arcane.
- Rendre les sorties parallèles
- Réfléchir à la pertinence d'une sortie au format binaire.

## Création de particule

Actuellement, particules créées en dur dans le module IPG, pour test de faisabilité.
Attention, le point d'entrée est actuellement sur start_init. Il faut peut-être envisager de le déplacer dans compute si la modélisation le requiert.

### Piste d'amélioration

- Créer une interface avec des services pour la création des particules
    - Générer des particules au hasard dans le gaz
    - Avec une loi d'éjection à partir de l'interface (voir papier LLNL)

- Prévoir conservation de la masse et autres quantités thermo


## Trajectographie (à créer)

- Pouvoir localiser les particules dans une maille
    - Ecrire un algo géométrique qui localise une particule dans une maille en fonction de ses coordonnées
    - Affecter la maille à la particule (voir doc Arcane)

- Gérer la trajectographie en parallèle
    - Changement de maille
    - Changement de sous-domaine

- Equations de mouvement, avec et sans rétro-action du gaz

[Doc Arcane de la classe Particle](https://arcaneframework.github.io/arcane/userdoc/html/dd/d6c/classArcane_1_1Particle.html)

## Reconstruction d'interface (à créer)

Actuellement, pente-borne avec étalement de l'interface sur plusieurs mailles. A + AB + AB + AB + B

Interface plus ou moins étalée suivant le limiteur choisi.

- Créer un nouveau service de projection (à côté de ADI + pente-borne) : ADI + Youngs. Objectif : A + A + A/B + B + B

Intérêt pour la reconstruction de la frontière gaz / lourd qui sert à l'éjection et au rattrapage.

En attendant, c'est toujours possible de faire des cas tests Lagrange, avec ou sans mailles mixtes.

## Initialisation de cas tests via un fichier de maillage

Actuellement, l'initialisation "Extenals" existe, mais est très limitée pour les applications. Monobloc + pas de mailles mixtes.

- Nécessite de pouvoir prendre en compte les mailles mixtes
- Nécessite une lib de projection pour projeter le maillage matériau sur le maillage de calcul (selon Gilles Grospellier : prévu prochainement).

### Branchement

- Nouveau service de l'initialisation dans casTest pour coder la fonctionnalité avec appel à la lib de projection

## Rattrapage (à créer)

- Suppression des particules qui sont rattrapées
- Prévoir conservation de la masse (et autres quantités thermo)

### Branchement
- Nouveau package avec interface et service
- Ajouter un point d'entrée dans le module Ipg.

## Calcul et sortie des nombres caractéristiques 
Mach, Weber, Knudsen, ...

### Branchement
- Créer les variables arcane correspondantes
- Calculer ces variables à chaque pas de temps là où c'est nécessaire
- Ajouter éventuellement ces variables à la liste des sorties

## Prise en compte du taux de vite dans les équations de conservation (à créer)
A spécifier