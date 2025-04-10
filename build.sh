#!/bin/bash

# Répertoire de build
BUILD_DIR="build"

# Créer le répertoire de build s'il n'existe pas
if [ ! -d "$BUILD_DIR" ]; then
    mkdir "$BUILD_DIR"
fi

# Accéder au répertoire de build
cd "$BUILD_DIR"

make clean

cmake ..

# Construire le projet
make

# Vérifier si la construction a réussi
if [ $? -eq 0 ]; then
    echo "La construction a réussi !"
else
    echo "Erreur lors de la construction."
fi

# Retour au répertoire d'origine
cd ..