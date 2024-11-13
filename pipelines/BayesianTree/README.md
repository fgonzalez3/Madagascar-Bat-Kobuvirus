# BayesianTree Pipeline Description - Building a Bayesian Time Tree using BEAST2 

The goal of this analysis was to estimate the timing of the most recent common ancestor (MRCA) of all Kobuviruses, as well as to estimate the time of divergence of the Madagascar bat Kobuvirus from other Aichivirus A lineages defined in our paper

# Selecting Representative Sequences 

To do this, I subselected from the kobuvirus ML phylogeny using Parnas to ensure adequate coverage of tree diversity. This yielded a subselection of kobuvirus sequences, including OP287812 and a Rabovirus outgroup which was later added 

# Alignment and Model Selection

After compiling sequences for this tree, we renamed the sequences in BEAST format, which can be found in the scripts directory. This included the accession number and the collection date. For cases where only a collection year was reported, we set the middle of the year (July 15th) as the corresponding date

# Building a Phylogenetic Tree Using BEAST2

The BEAST community maintains a number of helpful tutorials that you should practice before getting BEAST2 running, found [here](https://taming-the-beast.org/tutorials/Introduction-to-BEAST2/)

The first step in this process requires the generation of a .xml file for input into BEAST2, using the program BEAUTi. Some substitution models are easily specified in BEAUTi. For specifying less common models in BEAUTi, see this [blog post](https://justinbagley.rbind.io/2016/10/11/setting-dna-substitution-models-beast/), the recommendations [here](https://groups.google.com/g/ggplot2/c/H50aGubqt2U), [here](http://www.iqtree.org/doc/Substitution-Models), or [here](https://groups.google.com/g/beast-users/c/FH8OG_taajw). It is also possible to generate .xml files outside of BEAUTi, though this approach was not needed for this scenario

To prepare the .xml file, we used the following parameters in the tab inputs at the top of the screen in BEAUTi:

- **Tip Dates**: We used the date of the sample collection as the "Tip Date". For any sample from GenBank that only listed the year of collection, we set the tip date to July 15 of the year of collection. Each alignment was uploaded to BEAUTi with sequence names arranged so as to easily decipher the date. To specify the date, click "as dates with format" and leave as YYYY-MM-DD. 
- **Site Model**: Following outputs from ModelTest-NG, we selected a 'Gamma Site Model' with Gamma Category 4 and estimated the proportion of invariant sites by clicking the 'estimate' parameter. Proportion invariance was set to 0.001 following Cara's suggestion
- **Substitution Model**: The model suggested to us by ModelTest-NG was GTR+I+G4, which was luckily easy to select for within BEAUTi and done when we specified our site model. For more complex substitution models, such as a TPM2 model for example, you would need to link AC-AT, CG-GT, and AG-CT, while keeping the frequencies as 'estimated' and selecting the 'estimate' box next to the shape parameter to completely specify a TPM2uf+G4 model
- **Clock Model**: We built phylogenies using a strict molecular clock, specified within a non-parametric Bayesian Skyline Coalescent model, following previous approaches for bat coronaviruses analyses, [Lau et al.2020](https://journals.asm.org/doi/full/10.1128/JVI.02219-09). Make sure to match the clock rate value here to the lognormal distribution value in the next step (0.001). 
- **Priors**: We used a Bayesian Skyline Coalescent model. The clock rate prior was set to a lognormal distribution with a mean of 0.001, following published values for RNA viruses [Jenkins et al.2014](https://link.springer.com/article/10.1007/s00239-001-0064-3). All other priors were left at default values specified in BEAUTi
- **MCMC**: We used an MCMC chain of 700,000,000 iterations and set tracelog and treelog every 10,000 iterations. Additionally, the header for the output file for the treelog and tracelog should also be the same (e.g. kobu.log for trace log and kobu.trees for treelog). All other MCMC metrics were left at default
- **Population Size**: The Bayesian Skyline Coalescent model by default assumes that the population size of the dataset will change 5 times space evenly across the course of this sampling period. Because our available samples were limited and spanned a large geographic area, we edited this parameter to permit only one population size. You can make this edit in BEAUTi by clicking 'View' in the top panel, selecting 'Show Initialization Panel', and specifying the dimension of 'bPopSizes' and 'GroupSizes' both to 1 instead of 5. Also, enter a value >0 for the proportion invariant estimate. This should match our values in the Clock Model and Priors (0.001) rather than being 0

# Visualizing Bayesian TimeTree

The BEAST2 output was then manipulated preceding visualization. The initial 10% of MCMC iterations were removed as burn-in. Parameter convergence was assessed visually using [Tracer](https://www.beast2.org/tracer-2/). We used TreeAnnotator to average across the BEAST2 tree output by (a) setting the burn-in to 10% to discard the first 10% of trees in the log file, (b) posterior probability limit at default, (c) leaving target tree type at default maximum clade credibility, and (d) selecting mean heights for node height. This was then visualized in FigTree. After checking for basic alignment with parametric phylogenies generated from RAxML, we converted the Bayesian tree which is output in Nexus format to Newick format by exporting from FigTree. We then imported the resulting Newick file of the average tree in R and visualized it


