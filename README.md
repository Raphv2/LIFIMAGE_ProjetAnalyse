# Projet

## Prérequis

Avant de commencer, assurez-vous d'avoir les éléments suivants installés sur votre système :

- [CMake](https://cmake.org/download/)
- [Make](https://www.gnu.org/software/make/)
- Un compilateur compatible (GCC, Clang, MSVC, etc.)

## Compilation et exécution

Un script `build.sh` est fourni pour faciliter la compilation. Suivez ces étapes pour compiler et exécuter le projet :

### 1. Cloner le dépôt
```sh
git clone <url_du_repository>
cd <nom_du_projet>
```

### 2. Construire le projet

Exécutez le script `build.sh` pour configurer et compiler le projet :
```sh
./build.sh
```

Si la commande `build.sh` n'est pas reconnue ou génère une erreur d'accès refusé, essayez :
```sh
chmod +x build.sh
./build.sh
```

Si la compilation réussit, les exécutables seront générés dans le répertoire `bin`.

### 3. Exécuter les programmes

Une fois la compilation terminée, vous pouvez exécuter les programmes ainsi :
```sh
./bin/synthese
./bin/analyse
```

## Nettoyage du projet

Si vous souhaitez nettoyer les fichiers générés lors de la compilation, utilisez la commande :
```sh
rm -rf build
```
Cela supprimera le répertoire `build` et tout son contenu.

## Dépannage

Si la compilation échoue :
- Vérifiez que CMake et Make sont bien installés.
- Assurez-vous que votre compilateur est bien configuré et accessible depuis le terminal.
- Consultez les messages d'erreur pour identifier d'éventuels problèmes dans le code ou la configuration.
