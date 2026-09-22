# SVIB V4.1 — correctif de diagnostic (pas de moteur Tracy)

Cette version ne constitue PAS un variant caller. Elle conserve la comparaison exploratoire V4 et ajoute :
- messages explicites si la référence ou les alignements manquent ;
- remise à zéro des anciens alignements quand la référence change ;
- remise à zéro des anciens candidats quand les ABI sont rechargés ;
- affichage des exceptions au lieu d’un bouton silencieux ;
- version du script pour limiter les problèmes de cache GitHub Pages.

Ordre : charger NM_ → importer F/R → Lire et analyser localement → Comparer.

Aucune nomenclature d’indel n’est produite. Ne pas utiliser à des fins diagnostiques. Ne jamais publier de chromatogrammes de patients sur GitHub.
