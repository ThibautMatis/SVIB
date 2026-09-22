# SVIB V3 — prototype de recherche

Application statique pour GitHub Pages. Remplacer `index.html`, `style.css`, `app.js` à la racine du dépôt SVIB. Aucun FASTA à héberger.

Rechercher un symbole de gène humain, choisir un NM_ versionné, puis charger sa référence NCBI (séquence + annotation CDS). Possibilité de saisir directement un accession.version. La recherche NCBI est limitée à 200 résultats ; elle ne garantit pas une liste exhaustive de toutes les versions historiques. La récupération dépend du réseau, de CORS et de la disponibilité de NCBI ; un FASTA manuel reste possible sans annotation CDS automatique.

Les fichiers AB1 sont lus localement ; seules les références publiques sont demandées à NCBI. Les requêtes de référence peuvent révéler le gène/transcrit recherché au fournisseur réseau/NCBI. Aucun fichier patient ne doit être ajouté au dépôt GitHub.

**Limitations :** alignement heuristique et discordances SNV indicatives ; les délétions hétérozygotes ne sont pas appelées ni normalisées en HGVS. Le module indel explore des régions candidates uniquement. Le calcul c. des SNV est limité aux positions de CDS et reste indicatif. Pas de validation clinique ni de garantie de sécurité pour données patient.
