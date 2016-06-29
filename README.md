# Quantitative comparison of DNA methylation assays for biomarker development and clinical applications 


**Coordination & analysis team**: Christoph Bock, Florian Halbritter, Francisco J. Carmona, Sascha Tierling, Paul Datlinger

**Assay contributors (alphabetical order)**: Yassen Assenov, María Berdasco, Anke K. Bergmann, Keith Booher, Florence Busato, Mihaela Campan, Christina Dahl, Christina M. Dahmcke, Dinh Diep, Agustín F. Fernández, Clarissa Gerhauser, Andrea Haake, Katharina Heilmann, Thomas Holcomb, Dianna Hussmann, Mitsuteru Ito, Ruth Kläver, Martin Kreutz, Marta Kulis, Virginia Lopez, Shalima S. Nair, Dirk S. Paul, Nongluk Plongthongkum, Wenjia Qu, Ana C. Queirós, Frank Reinicke, Guido Sauter, Thorsten Schlomm, Clare Stirzaker, Aaron Statham, Ruslan Strogantsev, Rocío G. Urdinguio, Kimberly Walter, Dieter Weichenhan, Daniel J. Weisenberger

**Laboratory heads (alphabetical order)**: Stephan Beck, Susan J. Clark, Manel Esteller, Anne C. Ferguson-Smith, Mario F. Fraga, Per Guldberg, Lise Lotte Hansen, Peter W. Laird, José I. Martin-Subero, Anders O. H. Nygren, Ralf Peist, Christoph Plass, David S. Shames, Reiner Siebert, Xueguang Sun, Jörg Tost, Jörn Walter, Kun Zhang *for the BLUEPRINT consortium*

## Scientific Abstract:

DNA methylation patterns are altered in numerous diseases and often correlate with clinically relevant information such as disease subtypes, prognosis, and drug response. With suitable assays and after validation in large cohorts, such associations can be exploited for clinical diagnostics and personalized treatment decisions. Here we describe the results of a community-wide benchmarking study comparing the performance of all widely used methods for DNA methylation analysis that are compatible with routine clinical use. We shipped 32 reference samples to 18 laboratories in 7 different countries. These laboratories collectively contributed 21 locus-specific assays for an average of 27 predefined genomic regions, as well as 6 global assays. We evaluated assay sensitivity on low-input samples and assessed the assays’ ability to discriminate between cell types. Good agreement was observed across all tested methods, with amplicon bisulfite sequencing and bisulfite pyrosequencing showing the best all-round performance. Our benchmarking analysis can help with the selection, optimization, and use of DNA methylation assays in large-scale validation studies, biomarker development, and clinical diagnostics.

## Repository description:

Supplementary code repository accompanying the BLUEPRINT benchmark study of locus-specific DNA methylation assays.

## Repository structure:

There are only two directories:

+ ``/data``: tab-delimited text files containing all relevant assay data and metadata
+ ``/analysis``: a collection of R scripts, organized per figure of the article. ``methBench_analysis.R`` is an entry script that loads all required libraries and consecutively executes each step of the analysis

## Links:

+ Full article: [Nature Biotechnology](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3605.html)
+ PubMed abstract: [27347756](http://www.ncbi.nlm.nih.gov/pubmed/27347756)
+ Illumina 450k data deposited at the Gene Expression Omnibus (GEO): [GSE77965](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77965)
+ BLUEPRINT consortium homepage: [Blueprint Epigenome](http://www.blueprint-epigenome.eu/)

