# Campagne Jurin AMR — VM Linux

Ce dossier est une copie prête à transférer. Les essais locaux et leurs temps de
calcul ne sont pas inclus. Les simulations créent leurs résultats dans campaign/.

## Prérequis de la VM

Bash, GCC, Basilisk avec qcc, Python 3 avec NumPy et Matplotlib, ainsi que les
outils Linux setsid et sha256sum. Le script les vérifie avant le premier calcul.
Aucun logiciel système n'est installé automatiquement. Jupyter sert uniquement
à l'analyse ; il n'est pas obligatoire pour lancer les simulations.
/usr/bin/time est facultatif : s'il est présent, le pic de RAM est aussi mesuré.

Si qcc n'est pas dans le PATH et n'est pas installé dans ~/basilisk/src :

```bash
export QCC=/chemin/vers/basilisk/src/qcc
```

## Extraire et démarrer

Dans le répertoire contenant Jurin_AMR_VM.zip :

```bash
python3 -m zipfile -e Jurin_AMR_VM.zip .
cd jurin_amr_vm
sha256sum -c SHA256SUMS
bash run_campaign.sh pilot
```

Le pilote effectue quatre calculs : MAXLEVEL 8 et 9, U_TOL 1e-2 et 1e-3 m/s.
Chaque calcul est compilé puis exécuté en série, sur un seul cœur, sans vidéo.
Les journaux sont dans campaign/Lxx_Uxxx/ : compile.log, simulation.log,
validation.log. status doit indiquer complete après vérification des sorties.
Consulter temps.dat pour le temps réel (secondes) et le pic de RAM (kio, ou nan
si /usr/bin/time est absent). Les niveaux 11 et 12 peuvent être coûteux.
Si la VM utilise un ordonnanceur, respecter les règles de soumission de Sorbonne.

## Référence et campagne complète

Après avoir examiné le pilote, lancer ces deux commandes successivement :

```bash
bash run_campaign.sh reference
bash run_campaign.sh full
```

reference : MAXLEVEL=12, U_TOL=1e-5 ; résultats dans campaign/reference/.
full : 35 configurations, niveaux 6 à 12 et U_TOL=1e-2, 3e-3, 1e-3, 3e-4, 1e-4.
Les calculs du pilote déjà terminés et vérifiés sont réutilisés.
Paramètres fixes : MINLEVEL=5, F_TOL=1e-3, T_END=0.8 s, T_COMPARE=0.15 s.
Une seule commande de campagne doit fonctionner à la fois dans ce dossier.

Pour conserver une session de calcul après une déconnexion SSH, on peut lancer
la référence puis la campagne avec :

```bash
nohup bash -c 'bash run_campaign.sh reference && bash run_campaign.sh full' > campagne.log 2>&1 &
```

## Interrompre et reprendre

En premier plan, Ctrl+C arrête le calcul en cours et le marque interrupted.
Relancer la même commande pour reprendre : les calculs terminés sont vérifiés
et conservés ; le calcul interrompu ou échoué est recommencé depuis t=0.
Un échec arrête la campagne et reste identifiable dans status et les journaux.
Conserver les sources pendant toute la campagne : un changement du code AMR
invalide la réutilisation des anciens résultats.

## Renvoyer les résultats et les analyser

Depuis jurin_amr_vm/, après les calculs :

```bash
python3 -m zipfile -c resultats_jurin.zip campaign
```

Renvoyer resultats_jurin.zip. Le dossier campaign/ contient tous les fichiers
nécessaires à l'analyse : hauteurs, vitesse axiale, cellules, paramètres,
temps, empreinte du code et statut. Conserver l'arborescence des dossiers.
Extraire les résultats à côté de analyse_campagne.ipynb, puis lancer :

```bash
jupyter notebook analyse_campagne.ipynb
```

Le notebook ne lance pas Basilisk. Les calculs absents ou invalides restent
manquants ; les erreurs ne sont calculées qu'avec une référence valide.
La référence est provisoire : F_TOL=1e-3, MAXLEVEL=12 et l'échantillonnage fixe
1024 × 64 peuvent limiter la précision. Sa convergence reste à vérifier.

## Fichiers fournis

- jurin_axi_amr.c : simulation et exports numériques, en unités SI.
- jurin_amr_video.h : conservé car qcc parcourt cette inclusion optionnelle ;
  la campagne désactive son exécution et n'a besoin ni de FFmpeg ni des scripts vidéo.
- run_campaign.sh : compilation, calculs séquentiels et reprise.
- analyse_campagne.ipynb : trois figures et tableau récapitulatif.
- SHA256SUMS : vérification de l'intégrité de cette copie.
