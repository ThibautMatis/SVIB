# SVIB V7 — prototype de recherche (non validé)

Cette version remplace l'adaptateur Tracy incomplet et l'explorateur V6 par un ajustement indépendant de chaque chromatogramme ABI sur une référence NM_ (1 à 20 nt de décalage). Les quatre canaux sont analysés localement ; la référence publique est récupérée auprès de NCBI.

**Limites :** ce n'est PAS Tracy ni une déconvolution validée. Les séquences secondaires sont des estimations guidées par la référence, sans alignement allèle par allèle ni appel HGVS fiable. La recherche peut manquer des variants ou produire des artefacts. Ne pas utiliser pour le diagnostic. La duplication BRCA2 c.1813dup n'a pas été confirmée par ce moteur sur les ABI fournis.

## Installation
Téléverser `index.html`, `app.js`, `style.css` et `v7-engine.js` à la racine de GitHub Pages, puis actualiser la page. Charger la référence, importer les ABI, cliquer « Lire et analyser localement » puis « Explorer les signaux F et R ». Les fichiers ABI ne sont pas envoyés à un serveur.

## Vérification technique
Syntaxe JS contrôlée avec `node --check`. Pas de test navigateur de bout en bout ni de validation de sensibilité/spécificité. Les résultats ne constituent pas des variants identifiés.
