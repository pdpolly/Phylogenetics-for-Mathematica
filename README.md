# Phylogenetics-for-Mathematica
This add-in package for Mathematica performs basic phylogenetic functions, including reading and drawing Newick format trees, calculating phylogenetically independent contrasts, reconstructing ancestral values for continuous traits, performing random walks, and simulating continuous traits on phylogenetic trees.

<b>Installation.</b> The file is in Mathematica's ".m" format, which can be imported into Mathematica 6.0 and later (functions do not work in earlier versions of Mathematica). Install using the "Install" item on the "File" menu. Once installed, you must load the package like any other with the line <<PollyPhylogenetics, using either this suggested name or another.

<b>Changes in Version 6.x.</b>  Many new functions.  PhylogeneticRegression performs a PGLS linear regression, MapTraitsOntoTree produces a colored graphic showing estimated trait values on a phylogenetic tree, PruneTree extracts a tree for a selected subset of taxa in a larger tree, Kappa estimates Blomberg's K for a univariate trait, KappaMultivariate estimates Adams' K for a highly multivariate geometric morphometric data set, and SharedPhylogeneticHistory and MeanTreeDistance report statistics useful in the study of phylogenetic community assembly.  New and more efficient functions are also available for estimating rates and ancestral node values for traits.  Changed all "Modules" to "Blocks" for improved efficiency.

<b>Changes in Version 6.6</b> Updated <i>Pagel's Lambda[]</i> and added new function <i>LambdaMultiplication[]</i>.

Funding for development for this package has been provided by grants NSF EAR 1338298, the Robert R. Shrock fund at Indiana University, the Yale Institute for Biospheric Studies, and the Lilly Endowment through its support for the Indiana University Pervasive Technology Institute and the Indiana METACyt Initiative.

<b>Cite as:</b> Polly, P.D. 2022. Phylogenetics for Mathematica. Version 6.6. https://github.com/pdpolly/Phylogenetics-for-Mathematica
