********
AlphaMate
********

VERSION 1.00
============

.. contents:: Table of Contents
   :depth: 5

Introduction
============
Alphamate is a software package for optimizing genetic contributions TODO ...


List of contributors to the development and testing of AlphaMate
===============================================================
Gregor Gorjanc, John M Hickey, Chris R Gaynor.


Funding
=======
TODO


Availability
============
AlphaMate is available from: http://www.alphagenes.roslin.ed.ac.uk/software-packages/alphamate/

Material available includes the compiled programs for 64 bit Linux and Mac OSX machines, together with a User Manual.

Please report bugs or suggestions on how the program / user interface / manual could be improved or made more user friendly to John.Hickey@roslin.ed.ac.uk. TODO


Conditions of use
=================
AlphaMate is available to the scientific community free of charge. Users are required, however, to credit its use in any publications. Commercial users should contact John Hickey (John.Hickey@roslin.ed.ac.uk). TODO

Please cite the AlphaMate paper TODO


Disclaimer
==========
While every effort has been made to ensure that AlphaSim does what it claims to do, there is absolutely no guarantee that the results provided are correct. Use of AlphaMate is entirely at your own risk!


Advertisement
=============
AlphaMate is part of the AlphaSuite collection of software programs that we have developed. The AlphaSuite collection can perform many of the common tasks in animal breeding, plant breeding, and human genetics including genomic prediction, breeding value estimation, variance component estimation, GWAS, imputation, phasing, optimal contributions, simulation, field trial designs, and various data recoding and handling tools.

The AlphaSuite [#f1]_ is available at this link: http://www.alphagenes.roslin.ed.ac.uk/software-packages/


Description of the method
=========================
The method implemented in AlphaMate to TODO

More details about the methods in AlphaMate are provided in section `AlphaMateSpec.txt`_ below, which describes the parameters file.


Using AlphaMate
==============


Input files
-----------
Running AlphaMate requires at least file providing the program parameters (`AlphaMateSpec.txt`_), a file with breeding values, and a file with relationship matrix. In addition a file with gender information can be supplied. TODO PAGE, individual and mate-pair files


.. _`AlphaMateSpec.txt`:

AlphamateSpec.txt (mandatory)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
An example of *AlphaMateSpec.txt* is shown in `Figure 1`_. Text to the left of the comma should not be changed. The program is controlled by changing the input right of the first comma::

    Mode                                    ,MinThenOpt
    RelationshipMatrixFile                  ,RelMat.txt
    BreedingValueFile                       ,Ebv.txt
    GenderFile                              ,Gender.txt
    NumberOfIndividuals                     ,200
    NumberOfMatings                         ,50
    NumberOfParents                         ,20
    NumberOfMaleParents                     ,10
    NumberOfFemaleParents                   ,50
    EqualizeParentContributions             ,No
    EqualizeMaleParentContributions         ,No
    EqualizeFemaleParentContributions       ,Yes
    LimitParentContributions                ,No,5,20,0.01
    LimitMaleParentContributions            ,Yes,5,20,0.01
    LimitFemaleParentContributions          ,Yes,1,1,0.01
    AllowSelfing                            ,No,0.1
    OldCoancestry                           ,Unknown
    TargetedRateOfPopulationInbreeding      ,0.01,1,Above
    IndividualInbreedingPenalty             ,0.1
    EvaluateFrontier                        ,No,3,0.000001,0.00001,0.0001
    EvolutionaryAlgorithmIterations         ,100,20000,500,1000,0.001,50
    EvolutionaryAlgorithmParameters         ,0.4,0.2,0.1,1.0,4.0
    Seed                                    ,19791123

.. _`Figure 1`:

**Figure 1**. Example of AlphaMateSpec.txt

**Mode** specifies the mode of running AlphaMate.
* When mode is ``Min``, the program optimises contributions that would give minimum possible inbreeding. Once the solution is found, the specified level of old coancestry (see **OldCoancestry**) is corroborated.
* When mode is ``Opt``, the program optimises contributions that would give maximum genetic gain under **TargetedRateOfPopulationInbreeding**.
* When mode is ``MinThenOpt``, the program runs first the ``Min`` mode and then the ``Opt`` mode. The ``MinThenOpt`` mode requires the least amount of information from the user and should be the default mode if you are new to AlphaMate.

**RelationshipMatrixFile** specifies the file holding a relationship matrix. The first column must be individual identification and the subsequent columns are relationship coefficients. Relationship coefficients need to be "proper" - they should reflect probability of sharing genetic material that is identical, either by descent or by state. The relationship matrix can be built from either pedigree or genotype information. We suggest our AlphaAGH program for this task. When using genotype data, make sure to use the Nejati-Javaremi version of relationship matrix and no addition to the diagonal elements of the matrix. This is required to have proper connections between relationship and inbreeding coefficients and the rate of inbreeding within the AlphaMate.

**BreedingValueFile** specifies a file holding (estimated) breeding values. The first column must be individual identification and the second column breeding values.

**GenderFile**  specifies a file holding gender information. The first must be individual identification and the second column gender coded as 1 for males and 2 for females. When gender is not relevant, i.e., individuals can act both as male or female parents, use keyword *None* and no file need to be provided. Wether gender information is relevant, the program performs optimisation in a different way.

**NumberOfIndividuals** specifies how many individuals from the **RelationshipMatrixFile** and **BreedingValueFile** will be read by AlphaMate and used in the optimisation.

**NumberOfMatings** specifies the number of matings that the user wants to make.

**NumberOfParents** specifies the number of parents that the user can have. This can either be all **NumberOfIndividuals** or a smaller number. The latter would indicate that say a breeder can only afford to maintain a smaller set of individuals as parents.

**NumberOfMaleParents** is the same as **NumberOfParents**, but specific for males when gender must be considered.

**NumberOfFemaleParents** is the same as **NumberOfParents**, but specific for females when gender must be considered.

**EqualizeParentContributions** enforces that all parents (selected individuals with contributions larger than zero) have equal contributions.

**EqualizeMaleParentContributions** is the same as **EqualizeParentContributions**, but specific for males when gender must be considered.

**EqualizeFemaleParentContributions** is the same as **EqualizeParentContributions**, but specific for females when gender must be considered.

**LimitParentContributions** specifies minimum and maximum number of contributions allowed. This option is ignored when the first keyword is ``No``. When the first keyword is ``Yes``, the user needs to provide three more values/keywords: i) minimum number of contributions, ii) maximum number of contributions, and iii) a penalty when such limits can not be imposed. The penalty is added to the optimised objective for each individual that violates the constraints. See the section TODO to understand what level of penalty values should you use.

**LimitMaleParentContributions** is the same as **LimitParentContributions**, but specific for males when gender must be considered.

**LimitFemaleParentContributions** is the same as **LimitParentContributions**, but specific for females when gender must be considered.

**AllowSelfing** specifies whether to allow selfing (the first keyword must be ``Yes``) or not (the first keyword must be ``No``). Selfing is only possible when gender need not be considered. When selfing is not allowed, the second keyword must also be given to specify a penalty - AlphaMate tries to avoid selfing, but sometimes this is not possible for every explored solution during the optimisation and in such a case a penalty is added to the optimised objective for each selfed mating encountered. See the section TODO to understand what level of values for penalty should you use.

**OldCoancestry** specifies average coancestry among the parents of individuals provided to AlphaMate. This can be either a number, say 0.125, or a keyword ``Unknown``. When keyword ``Unknown`` is specified, AlphaMate attempts to estimate the coancestry from average inbreeding coefficients of individuals provided to AlphaMate. This information is required to compute the rate of inbreeding (the population inbreeding part of objective in AlphaMate). When information about the **OldCoancestry** is not available or imorecise, we propose to use keyword ``Unknown`` and set the **Mode** to ``Min`` or ``MinThenOpt``.

**TargetedRateOfPopulationInbreeding** specifies the rate of inbreeding that AlphaMate should achieve when optimising genetic contributions of selected individuals. The second value specifies the penalty when the rate of inbreeding is either too low or too high.
See the section TODO to understand what level of values for penalty should you use. this penalty is very important and values around 1.0 guide optimisation to the desired rate of inbreeding (lower values allow for more exlporation, but also less weight on inbreeding versus genetic gain). The third value specifies if penalty should be applied when the achieved rate of inbreeding is larger then targeted ``Above`` or also when the the achieved rate of inbreeding is lower than targeted ``AboveAndBellow``. The latter is always switched on when evaluating the frontier.

**IndividualInbreedingPenalty specifies a penalty for inbreeding in future progeny given the proposed mating list. While the option **TargetedRateOfPopulationInbreeding** enables control of the genetic diversity of the whole population, this option enables penalising mating lists that give high inbreeding of individual progeny, perhaps with the aim to mitigate inbreeding depression. See the section TODO to understand what level of values for penalty should you use. This penalty is likely to be specific for each example - depending on the level of inbreeding values. If this penalty is set to zero, the produced matings are essentialy random matings.

**EvaluateFrontier** specifies if AlphaMate should evaluate frontier of Pareto equivalent solutions (the first keyword must be ``Yes``) or no. If ``Yes``, then the second keyword specifies the number of points on the Pareto curve, followed by the rates of inbreeding for each point on the curve.

**EvolutionaryAlgorithmIterations** controls the evolutinary algorithm. The first keyword specifies how many solutions should the algorithm test every iteration - values around 100 seems to be ok. The second keyword specifies the maximum number of iteration the algorithm should run. The third keyword specifies the number of so called "burn-in" iterations, where the aglgorithm is exploring solutions with "large/bold steps". The fourth and fifth keywords specifies the number of iterations to stop optimisation if objective is not improving for a given precision. The sixth keyword specifies an interval to print optimisation results.

**EvolutionaryAlgorithmParameters** control the evolutinary algorithm "engine" (these are advanced options and should only be changed with great caution). The keywords are: cross-over parameter in the "burn-in" phase and after it, the base mutation parameter and two more values for the mutation parameter to enable large/bold jumps in the optimisation - can help with local optima.

**Seed** specifies the seed value for random number generation. If ``None`` is given, a seed value is generated by the program.


Output files TODO STOPPED HERE
------------------------------
The output of AlphaMate is organised in three directories (*Chromosomes*, *Selection* and *SimulatedData*). A description of what is contained within each of these directories is given below.

In addition to output files, AlphaSim prints statements in the terminal. This output contains information about the process, some synthetic results, error statements and notes. Error statements stop the simulating process. In contrast, notes do not stop the simulating process; however, they can be useful to track an error that would not have been indicated by the program.


.. _`Chromosomes`:

Files created by AlphaSim
"""""""""""""""""""""""""
* *BaseIndividualsGameteTracking.txt* includes five columns: identifier of the individual, identifier of each of the two paternal haplotypes, and identifier of each of the two maternal haplotypes. The file is filled in for the individuals of the first generation only (i.e., the founders).
* *Chip"e"GenotypeGeneration"g".txt* includes the genotypes of each individual of generation *g* on chip *e*. The files have as many rows as there are individuals in the  *g*\ :sup:`th`\  generation, and as many columns as there are SNP on the eth chip plus one, for the individual ids (first column).
* *Chip"e"PhaseGeneration"g".txt* includes the sequence phases of each individual of generation *g* on chip *e*. The files have twice as many rows as there are individuals in the  *g*\ :sup:`th`\  generation (two lines per individual, one for each haplotype), and as many columns as there are SNP on the eth chip plus one, for the individual identifiers (first column).
* *Chip"e"SnpInformation.txt* includes information about the SNP present on chip *e*. There are four columns: line identifier, SNP identifier, its minor allele frequency, and its physical position on the chromosome.
* *Chip"e"Summary.txt* includes only one number, which is the mean heterozygosity of the SNP found on chip *e*. The hererozygosity is obtained as follows:

  .. math:: Heterozygosity = \dfrac{\sum_{j=1}^{n_{SNP}}\left[1-\left(p_j^2+(1-p_j)^2\right)\right]}{n_{SNP}}

  where *p* is the frequency of the non-zero allele at SNP *j*.

* *GenomeGeneration"g".txt* stores the genomic data of generation *g*. The data are under the form of a three-dimensional array, the three dimensions corresponding to the individual, the segregating site and the haplotype haplotype, respectively. The file is unformatted, it is readable by AlphaSim only.
* *MacsHaplotypesInBlocks.txt* stores the data of each haplotype expressed in binary numbers (see details in Description of the method for block definition). This file is unformatted, it is readable by AlphaSim only.
* *MacsHaplotypesIncluded.txt* includes the sequences of each of the haplotypes generated by MaCS. The file has as many rows as there are haplotypes and as many columns as there are segregating sites that are included (see details in Description of the method for "included" segregating sites).
* *MinorAlleleFrequency.txt* includes two columns: the identifier of the segregating site and the allele frequency.

  The minor allele frequency (:math:`MAF`) is computed for each segregating site based on the sequence phases in the first generation (i.e., the founder’s haplotypes) as follows:

  .. math:: MAF = \dfrac{\sum_{i=1}^{n_{indiv}} Phase_1 + Phase_2}{2n_{indiv}}

  where :math:`Phase_1` and :math:`Phase_2` are the alleles on the first and second haplotypes, and :math:`n_Indiv` is the number of founders. If :math:`MAF` is smaller than :math:`0.5`, then :math:`MAF` is substituted by :math:`1-MAF`.


An example of AlphaMate use
---------------------------
The provided parameters file can be used as an example of AlphaMate use. In the example, we TODO


Background reading
==================

#. TODO

.. rubric:: Footnotes

.. [#f1] People often ask as to the origins on the AlphaSuite naming convention. The original program in the AlphaSuite was AlphaBayes, which was named after the “Bayesian Alphabet” that it implemented. The name was suggested by Brian Kinghorn. Subsequent programs evolved via mutation with a conserved primer!

