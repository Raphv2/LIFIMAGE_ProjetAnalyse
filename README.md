# LIFIMAGE – Synthèse & Analyse d’image

Projet réalisé dans le cadre du module **LIFIMAGE**. Il se compose de deux volets :

- `synthese.cpp` → un moteur de lancer de rayons (ray tracing)
- `analyse.cpp` → un outil d’analyse et traitement d’image

## 📁 Structure du projet

```
/bin/                → exécutables générés
/data/               → modèles .obj et fichiers image
/src/                → code source principal
include/             → en-têtes (gKit)
lib/                 → gKit3
```

## Prérequis

- CMake ≥ 3.22
- Un compilateur C++ compatible C++20
- OpenCV (pour les effets post-traitement)

## Compilation

Un script est fourni pour simplifier la compilation :

```bash
chmod +x build.sh
./build.sh
```

Les exécutables sont générés dans `bin/`.

### En cas de conflit CMake

Si vous avez cette erreur :
```
CMake Error: The source "..." does not match the source "..."
```

Cela signifie que vous avez lancé cmake dans un dossier build/ déjà utilisé par un autre projet.

Dans ce cas, supprimez le dossier build/ et recommencez :

```
rm -rf build/
mkdir build
cd build
cmake ..
make
```

## Exécution

```bash
./bin/synthese   # lance le moteur de ray tracing
./bin/analyse    # lance les fonctions d’analyse d’image
```

## Nettoyage

```bash
rm -rf build
```

## Auteurs

- **Yanis LAASSIBI**
- **Raphaël GOSSET**