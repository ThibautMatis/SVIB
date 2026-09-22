# SVIB v2 — prototype de recherche (non diagnostique)

## Installation sur GitHub Pages
Déposer les fichiers `index.html`, `style.css`, `app.js` à la racine du dépôt `SVIB`. Ajouter également le fichier **exact** `BRCA2_NM_000059.4.fasta` téléchargé depuis NCBI (NM_000059.4, 11954 nucléotides) à la racine. Le chargement de la référence est automatique après rechargement du site. Si le fichier est absent, l’import manuel reste possible.

Le champ CDS doit être renseigné d’après les annotations CDS de l’enregistrement GenBank de la même version ; ne pas inférer l’ATG en recherchant le premier ATG de la séquence. Les substitutions affichées restent indicatives.

## Navigation
Glisser horizontalement le chromatogramme ; double-cliquer pour zoomer ; bouton Dézoomer ou Vue complète ; molette pour ajuster. Cliquer sur une base pour ses coordonnées approximatives.

## Recherche d’indels
Le bouton « Rechercher des régions candidates » effectue un **classement heuristique exploratoire**, sur les intensités des pics et des décalages de quatre bases, indépendamment pour chaque orientation. Les positions ne sont **pas** des appels d’indels ni des nomenclatures HGVS ; la déconvolution complète des deux allèles et la validation analytique ne sont pas implémentées. Les scores ne sont pas des probabilités.

## Confidentialité
Les fichiers .ab1 sont lus dans le navigateur, sans upload programmé. Le FASTA public est téléchargé depuis le même site GitHub Pages. Ne jamais committer de fichiers patients dans le dépôt public. Pour des données cliniques, auditer le code, les dépendances et l’environnement avant utilisation.
