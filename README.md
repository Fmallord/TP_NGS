# TP_NGS
-------------------------------------------
## Objectif

Ce TP avait pour but final de déterminer comment les mutations délétères sont associés dans le génome humain, et ce en essayant de recréer l'article "Negative selection in humans and fruit flies involves synergistic epistasis", de Sohail et al (2017).      
Cet article concluait que les mutations délétères se repoussaient au niveau individuel, car l'association des mutations délétères provoquait une chute encore plus significative de fitness (= synergie épistatique, opposée à l'épistasie antagoniste où les mutations délétères associées ont un impact négatif global atténué, et opposée aussi à la neutralité où les effets des mutations délétères s'additionnent).  
Pour tester cette hypothèse, nous avons donc décidé d'analyser les génomes de nombreux individus, afin de déterminer le type d'association présent entre leurs mutations délétères.   

-------------------------------------------
## Première semaine: test des fichiers vcf (= fichiers comportant les SNPs génomiques)

Durant la première semaine de TP, nous avons tout d'abord extrait le génome de 3 individus vietnamiens formant un pedigree père-mère-fille sous forme de raw data (=reads du 1000G Project). Ensuite, nous avons aligné ces reads sur un génome de référence, puis nous avons effectué le mapping et la correction de ces alignements via différents programmes présentés ci-dessous, ainsi que le variant calling final des 3 génomes.
L'objectif était de savoir si le fichier vcf file obtenu pour la fille comportait un nombre d'erreurs cohérent avec ce qui était biologiquement attendu par rapport à ses parents (= 100 mutations de différence).
Afin de réduire la charge de mémoire, nous n'avons comparé et aligné que les séquences de l'exome du chromosome 20 des 3 individus.
Le résultat était biologiquement cohérent en termes d'erreur, ce qui a permis de valider les vcf files du 1000G Project comme produites de manière cohérente et proposant une correction (faite comme nous nous l'avions faite, sur le principe) réaliste. 
Nous pouvions donc nous permettre d'utiliser les fichiers vcf de tous les individus du 1000 Genome Project pour essayer de répondre à notre question en deuxième semaine.

-------------------------------
Programmes utilisés (semaine 1)  (dans l'ordre)
-------------------------------
installation.sh : permet d'installer des programmes qu'on utilisera par la suite (bedtools pour corriger les alignements par ex).
mapping.sh : permet d'aligner les reads de chaque individu contre le génome de référence, afin d'en sortir un fichier .sam au final (format de séquence ADN (et aussi d'acides aminés) visualisable dans le logiciel IGV).
trio-analysis.sh : permet de phaser les séquences des individus et d'obtenir les vcf files.
variant-calling.sh : permet de corriger les alignements et de réajuster les variants en les comparant aux variants connus et aux variants des parents (pour la fille); fait les vcf files finaux.

-------------------------------------------
## Deuxième semaine: analyse des données exomiques du 1000GP

Durant la deuxième semaine du TP, nous avons analysé de multiples fichiers vcf (tout ceux du 1000GP, grâce aux analyses de Thibault), sur l'exome de tout le génome.
Afin d'avoir une idée de l'effet délétère de chaque mutation, nous les avons classé par type (synonyme, non-synonyme, non sens) considérés respectivement de plus en plus délétère.
Nous avons ensuite calculé pour chaque catégorie de SNPs le rapport sigma²/Va, avec sigma²=la variance de la somme des SNPs délétères par individu, et Va=la somme des variances par SNP délétère, dans la population considérée.
La différence non significative entre Va et sigma² permet de conclure à l'indépendance de l'effet des mutations délétères. Au contraire, la dépendance entre les effets négatifs des mutations est signalée par un écart significatif entre Va et sigma² (voir début de la partie "Résultats de la deuxième semaine").
Pour ce faire, nous avons trié les individus par population, leurs mutations par les catégories indiquées plus haut, et nous avons calculé les statistiques sigma² et Va pour toutes les populations et méta-populations du 1000GP.
Les statistiques ont été calculées pour des seuils fixes de nombre d'individus porteurs de mutations (= considération de toutes les mutations présentes au plus 1 fois dans la population, ou 5 fois, ou 50 fois par ex).

-------------------------------
Programmes utilisés (semaine 2)  (dans l'ordre)
-------------------------------
extract_pop.py : extraire les génomes par population.
analysis.sh : le principal fichier de la deuxième semaine, il permet de passer des fichiers vcf individuels à des fichiers de SNPs de l'exome classés dans des fichiers par population, et par type de mutation.
gtf_to_bed.py : permet de ne conserver que les exons non vides, référencés et "vrais" (= nombre de nucléotides multiples de 3).
vcf_coding_polymorphism.py : permet de classer les mutations par type de polymorphisme (syn, non-syn, ou LoF=stop pour nous).
vcf_analysis.py : le programme qui calcule Va, sigma², et leur ratio et détermine s'il est significativement différent de 1.
vcf_meta_analysis.py : script qui automatise les analyses populationnelles du ratio sigma²/Va. Cela fait des bootstrap sans remise des SNPs utilisés par populations, par type de mutation, avec tirage sans remise de 10 SNPs à chaque fois (pour les mutations syns et non-syns seulement) et calcul du ratio sigma²/Va à chaque fois. Ces ratios sont représentés graphiquement par des grpahes faits par le même programme.


## Résultats de la deuxième semaine

Nous avons pu constater que la synergie épistatique (sigma²/Va<1), l'épistasie antagoniste (sigma²/Va>1) et l'hypothèse de neutralité de répartition des mutations délétères (sigma²/Va=1) pouvaient toutes se retrouver, avec une nette dominance de la neutralité.
Cependant, les résultats dans une même population dépendaient du seuil de mutation choisi; les effets non neutres prédominaient aux petits seuils (là où se trouvent les mutations susceptibles d'être fortement délétères) et aux gros seuils en considérant tous les SNPs (probablement des effets d'haplotypes, qui entraîne un biais dans l'indépendance des mutations considérées).
Notamment, les populations sud-américaines présentaient de manière assez significatives et persistantes un effet de synergie.
Cette différence étrange entre les processus de sélection en Amérique et ailleurs serait plus probablement un biais, à chercher du côté de la population américaine actuelle, qui est profondément métissée, avec des apports d'Europe, d'Amérique et d'Afrique.
Ces effets de synergie pourraient être des effets comparables au biais qui consiste à choisir des individus de différentes populations pour l'analyse; on risque de trouver de la synergie car des mutations populations spécifiques ne se retrouvent pas dans les 2 populations.
La profonde divergence entre les 3 sources de la population amérindienne en fait une population plus profondément diverse que toutes les autres, et se réduire à une population amérindienne ne supprime pas le fait que des individus dans cette population ont des haplotypes d'afrique, d'europe ou d'amérique population-spécifiques de ces sources, ce qui crée donc un effet d'exclusion artéfactuelle entre mutations entre individus d'une même population, mais partageant des haplotypes d'origines profondément différentes au mêm locus.
Peut être refaire ces analyses sur de l'ADN ancien pré-colombien permettrait de trancher sur cet effet américain.

-------------------------------------------
## Conclusion

En tout cas, la neutralité est très majoritaire (>80% des cas), contredisant le papier de Sohail en se basant sur des données populationnelles plus fournies, avec des seuils de mutations considérés plus variés.
