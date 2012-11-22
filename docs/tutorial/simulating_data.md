# Simulating non-homogeneous TdG12 models

In the simplest case, say you have two clades for which you want to generate synthetic data. The program uses the first two alpha-numeric characters of the taxon name to designate the grouping of the taxa. Imagine we want to generate simulated data for the influenza PB2 protein, with avian and human virus each having their own substitution models. Each taxa in the PB2 tree (etc/PB2_FMutSel0.tree) starts with either "Av" or "Hu", indicating that the sequence is either an avian or human virus protein, respectively.

1. *Preparing the tree*. The tree needs to be fully annotated. We do this by running:

	`java cp dist/tdg12.jar tdg.tree.TreeNodeLabeler PB2_FMutSel0.tree`
	
	The TreeNodeLabeler program reads the original tree and writes a new tree (with suffix ".out") in which every internal node has been labeled either "Av" or "Hu". Additionally, there will be a node labeled "Av_HS", indicating that this internal node is connected to a child node with a different label (i.e. we're going from "Av" to "Hu"). The non-homogenous model will use the midpoint of this branch as the point at which substitution models will change. It is recommended that you open the new labeled tree in a tree viewer program such as Figtree and check all the nodes are labeled as you expect. If they aren't, you can manually edit the tree or original Newick file to split the tree as you like.
	
2. *Simulating a single set of fitness*. Say that you want to simulate a single site 100 times having residues P, S, T and H. In the avian branches, the fitness of these amino acids is 0,-1,-0.2,-1.9, whilst in the human branches they are 0,-0.2,-3.1,-3.5. In this case, we would run:

	`java -cp dist/tdg12.jar tdg.sim.AlignmentSimulator -tree PB2_FMutSel.tree.out -sites 100 -heteroclades Av,Hu -fitness 0,-1,-0.2,-1.9 -fitness 0,-0.2,-3.1,-3.5 -characters P,S,T,H -tau 1e-2 -kappa 6.5 -pi 0.25,0.25,0.25,0.25 -mu 2.0 -gc standard -output sitesim.phyl`
	
	where:
	
	**-heteroclades** gives the 2-letter labels that we are using the in tree (the root of the tree is always assumed to be the first group in this list).
	
	**-fitness** (as many -fitness parameters as there are groups). The order in which the -fitness parameters are given should match the order that you specified the group names -heteroclades.
	
	**-characters** specifies the corresponding amino acids for the list of fitness values.
	
	**-tau** **-kappa** **-mu** **-pi** are the nucleotide mutational model parameters.
	
	**-gc** is the genetic code (standard, vertebrate_mit or plastid).
	
	**-output** specifies the filename of the simulated data.
	
	After completion, you should see:
	
    ```
tdg.cli.GeneticCodeConverter - Using standard genetic code.
Av has fitness [0.0, -1.0, -0.2, -1.9] for residues [14, 15, 16, 8].
Amino acid frequencies are:
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.031, 0.0, 0.0, 0.0, 0.0, 0.0, 0.41, 0.23, 0.33, 0.0, 0.0, 0.0
Hu has fitness [0.0, -0.2, -3.1, -3.5] for residues [14, 15, 16, 8].
Amino acid frequencies are:
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007, 0.0, 0.0, 0.0, 0.0, 0.0, 0.44, 0.54, 0.02, 0.0, 0.0, 0.0
```

3. *Simulating an alignment of different sites*. If you have many sites you wish to simulate you can put the fitness parameters into files. You need as many files as groups (2 in this example). Each line in the file should have a space-separated list of 20 fitness values. They are expected to be in the canonical amino acid order (i.e. ARNDCQEGHILKMFPSTWYV). You can then synthesize an alignment by running:

	`java -cp dist/tdg12.jar tdg.sim.AlignmentSimulator -tree PB2_FMutSel0.tree.out -heteroclades Av,Hu -fitnessfile fitness.av.txt -fitnessfile fitness.hu.txt -tau 1e-2 -kappa 6.5 -pi 0.25,0.25,0.25,0.25 -mu 2.0 -gc standard -output alignmentsim.phyl`
	
	where the options are as above, but:

	**-fitnessfile** (as many -fitnessfile parameters as there are groups). The order in which the -fitnessfile parameters are given should match the order you specified the group names in -heteroclades.


