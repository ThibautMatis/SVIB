# SangerVariant — prototype web sans serveur

Application statique : `index.html`, `style.css`, `app.js`. Les chromatogrammes sont lus en mémoire par JavaScript, sans requête réseau applicative ni dépendance distante. Hébergement possible sur GitHub Pages ; aucune donnée de patient ne doit être déposée dans le dépôt.

## Publication depuis Windows (sans installer de logiciel)

1. Créez un dépôt GitHub **public** `SangerVariant` (le code sera public ; jamais les données de patients).
2. Dans le dépôt, `Add file` → `Upload files` ; glissez **les trois fichiers** `index.html`, `style.css`, `app.js` (et ce README), puis `Commit changes`.
3. `Settings` → `Pages` → `Build and deployment` → `Deploy from a branch` → `main` et `/(root)` → `Save`.
4. Ouvrez `https://VOTRE_IDENTIFIANT.github.io/SangerVariant/` une fois le déploiement terminé.
5. Chargez un FASTA **vérifié** de `NM_000059.4` depuis NCBI RefSeq ; le fichier n'est pas fourni. Pour une annotation c. indicative, renseignez la position du A de l'ATG initiateur dans la séquence chargée après vérification de la CDS officielle.
6. Importez les fichiers `.ab1` directement dans la page. **Ne les importez jamais sur GitHub.**

## Limites importantes

- Le programme lit ABIF (PBAS, PLOC, FWO_, DATA9–12), affiche les signaux et produit un alignement local simple des bases appelées.
- La table liste des **discordances indicatives** de substitutions simples, pas des variants confirmés. Les positions c. ne sont calculées que si l'utilisateur fournit la position CDS ; aucune validation HGVS externe, aucune nomenclature protéique.
- **Aucune détection fiable des indels hétérozygotes**, aucune déconvolution des pics superposés, aucune normalisation HGVS 3′, aucune analyse des exons/introns sur la référence transcriptomique. Le résultat attendu sur les deux fichiers BRCA2 avec délétion de quatre bases **n'est pas encore implémenté**.
- Un alignement local peut produire des discordances artificielles, notamment après un indel ou aux extrémités. Ne pas utiliser en diagnostic ; valider les résultats par une méthode indépendante.
- Une application chargée depuis GitHub Pages reçoit son code du serveur GitHub : le traitement est local, mais cela n'établit pas à lui seul une garantie de confidentialité. Examiner le code, les dépendances et la politique de sécurité du poste avant d'y charger des données cliniques.
