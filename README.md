# SVIB V6 — moteur JavaScript exploratoire, sans Tracy

Installation : déposer `index.html`, `style.css`, `app.js` et `v6-engine.js` à la racine du dépôt GitHub Pages. Ne pas publier les fichiers AB1. Recharger la page (Ctrl+F5), charger un transcrit NM_ puis importer F et R et cliquer « Lire et analyser localement », puis « Analyser F et R séparément ».

Cette version remplace le pont Tracy inopérant par un moteur JavaScript exécutable localement. Chaque lecture est alignée séparément sur le NM_ puis le module explore des transitions et des décalages de 1–20 bases sur les quatre canaux. Les candidats restent des **régions exploratoires**, sans appel de variant ni reconstruction complète des deux haplotypes ; les nucléotides insérés ne sont pas identifiés. Une table vide ne permet pas d'exclure un variant. Ne pas utiliser pour le diagnostic.

Contrôles réalisés : syntaxe JavaScript via `node --check`. Le moteur n'a pas été validé en navigateur ni sur les fichiers BRCA2 transmis. Aucun résultat HGVS ne doit être déduit des positions approximatives.

Confidentialité : traitement des ABI en mémoire dans le navigateur. Les requêtes NCBI servent exclusivement au téléchargement des références publiques. Le logiciel ne téléverse pas les chromatogrammes ; vérifier néanmoins les politiques de sécurité et d'hébergement avant un usage clinique.
