# SVIB V4 — analyse exploratoire de délétion de 4 nt

Application statique pour GitHub Pages : remplacer `index.html`, `app.js`, `style.css` à la racine du dépôt SVIB. Ne pas publier de fichiers patients.

## Nouveautés
- Recherche de candidats 4 nt sur chaque lecture avec comparaison d'un modèle WT seul à un modèle de mélange WT/délétion sur les intensités au sommet des pics.
- Comparaison des deux orientations à l'aide de leurs coordonnées sur le transcrit ; régions concordantes mises en évidence.
- Tableau de candidats cliquable pour inspecter le chromatogramme et exporter les résultats exploratoires.
- La référence et ses annotations CDS sont chargées à la demande depuis NCBI (requêtes publiques uniquement).

## Limites importantes
La recherche est heuristique : les pics peuvent être décalés, les bases appelées après un indel peuvent être incorrectes et l'alignement local peut dériver. Un score ou une concordance F/R ne constitue pas une preuve de variant. Cette version ne détermine **pas** les quatre bases délétées et ne produit **pas** de nomenclature HGVS d'indel. L'outil n'est pas validé pour le diagnostic. Utiliser des données anonymisées pour les essais. Les chromatogrammes ne sont pas envoyés à NCBI ; vérifier néanmoins le code et l'environnement avant utilisation clinique.
