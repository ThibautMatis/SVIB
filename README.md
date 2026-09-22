# SVIB V5 — adaptateur Tracy WebAssembly (NON FONCTIONNEL sans compilation)

Cette archive **ne contient pas de moteur Tracy**. Elle conserve SVIB V5 et ajoute un adaptateur JS qui charge exclusivement `tracy/tracy.js` et `tracy/tracy.wasm` depuis le même site, écrit la référence et les ABI dans le système de fichiers virtuel de WebAssembly, puis lance **séparément** `tracy decompose -r /reference.fa -o /forward /Forward.ab1` et son équivalent Reverse. Les sorties `.align1` et `.align2` sont affichées si disponibles. Aucun variant ni HGVS n'est annoncé par ce prototype.

## État du portage

Le code source Tracy n'a pas pu être récupéré depuis l'environnement de construction (accès réseau GitHub indisponible) et Emscripten n'y est pas installé. **Aucun fichier `.wasm` n'a été produit, aucune compilation ni exécution Tracy n'a été validée.** Le bouton de diagnostic indique donc explicitement « Tracy non disponible » jusqu'à présence d'un port fonctionnel. Ne pas publier ceci en remplacement de SVIB V5 pour obtenir du variant calling.

## Ce qu'il reste à faire avant intégration

1. Récupérer le dépôt officiel Tracy avec ses sous-modules, conserver les mentions BSD-3-Clause et établir une compilation native de référence.
2. Installer Emscripten et porter ses dépendances Boost, HTSlib et SDSL-lite ; adapter les appels POSIX éventuels. Compiler une factory ES module avec système de fichiers virtuel, `callMain` exporté, `INVOKE_RUN=0`, sans téléchargement ni réseau. Le nom de sortie attendu par cet adaptateur est `tracy/tracy.js` et `tracy/tracy.wasm`.
3. Tester le module en navigateur sur des fichiers **non cliniques**, puis sur les ABI BRCA2 dans un environnement autorisé. Vérifier que les sorties `.align1/.align2` existent, que la référence est correctement ancrée et que les deux orientations donnent des résultats cohérents. Adapter l'interface à la forme exacte des sorties de la version compilée.
4. Implémenter ensuite le rapprochement des variants F/R et la normalisation HGVS locale. Ne pas confondre `.align1/.align2` avec un variant validé.

La documentation de la commande native est : `tracy decompose -r sequence.fa -o outprefix input.ab1` (https://www.gear-genomics.com/docs/tracy/cli/). Dépôt : https://github.com/gear-genomics/tracy .

**Confidentialité :** ce module n'envoie pas d'ABI vers un serveur. La recherche de transcrits existante continue d'interroger NCBI pour des références publiques. Ne pas publier d'ABI sur GitHub.
