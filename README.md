# SVIB V5 — prototype de modèles bi-alléliques

Déployez `index.html`, `app.js`, `style.css` à la racine du dépôt GitHub Pages. Chargez un NM_ puis les ABI, cliquez sur Lire et analyser localement, puis sur Tester les modèles bi-alléliques.

## Limites scientifiques importantes

Ce moteur est une implémentation originale de **comparaison de modèles de mélange guidés par la référence**, pas Tracy ni une déconvolution générale validée. Il teste des décalages de 1–12 nt et produit des fenêtres de séquences *modélisées*. Pour une insertion, les bases insérées sont représentées par N : elles ne sont pas reconstruites. Les substitutions sont des discordances indicatives issues de l’alignement initial, non des appels hétérozygotes fiables. Le moteur ne calcule pas d’HGVS d’indel. Les coordonnées peuvent être fausses si l’alignement primaire échoue après un indel ; la concordance Forward/Reverse n’est pas validée. Un tableau vide n’exclut pas de variant. Aucun usage diagnostique.

## Confidentialité

Les fichiers ABI sont traités en JavaScript dans le navigateur ; les recherches de références publiques sollicitent NCBI. N’envoyez pas de fichiers de patients sur GitHub.
