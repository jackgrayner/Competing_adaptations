Fitness values should be changed per scenario as follows. Both mutations only have fitness effects in males, so the dominance coefficient does not apply to the X-linked mutation which is always expressed.
- additive fitness benefits (i.e., individuals expressing both have greater fitness):
	 - fullcomb=1.6; where both mutations are fully expressed
	 - heterozygouscomb=1.3+(0.3*dominancecoef); where autosomal mutation is heterozygous, co-expressed with the X-linked mutation
	 - fullsingle=1.3; where a single mutation is fully expressed
	 - heterozygousauto=1.0+(0.3*dominancecoef); where autosomal mutation is heterozygous, and X-linked mutation absent
- non-additive fitness benefits (i.e., individuals expressing both have no additional fitness benefit vs. those expressing one):
	 - fullcomb=1.3
	 - heterozygouscomb=1.3
	 - fullsingle=1.3
	 - heterozygousauto=1.0+(0.3*dominancecoef)
- negative fitness benefits (i.e., individuals carrying both have reduced fitness vs. those expressing one):
	 - fullcomb=1.2
	 - heterozygouscomb=1.3-(0.1*dominancecoef)
	 - fullsingle=1.3
	 - heterozygousauto=1.0+(0.3*dominancecoef)
- for single-mutation sims just comment out lines 39 & 40 or 41 & 42, respectively, and use:
	- fullsingle=1.3
 	- heterozygousauto=1.0+(0.3*dominancecoef)

To treat both mutations as autosomal, remove the modifyChild callback.
```
initialize() {
	initializeMutationRate(0);
	initializeMutationType("m1", 1.0, "f", 0.0);// auto
	initializeMutationType("m2", 1.0, "f", 0.0);// X marker
	initializeMutationType("m3", 1.0, "f", 0.0);// Y marker
	initializeMutationType("m4", 1.0, "f", 0.0);// adaptive auto (Cw)
	initializeMutationType("m5", 1.0, "f", 0.0);// adaptive X (Fw)
	m4.convertToSubstitution = F;//prevent conversion as this will remove epistatic effects once either is fixed
	m5.convertToSubstitution = F;
	initializeGenomicElementType("g1", m1, 1.0); 
	initializeGenomicElementType("g2", m2, 1.0); 
	initializeGenomicElement(g1, 0, 99999);	//this is the 'autosomal' genomic element  
	initializeGenomicElement(g2, 100000, 199999); //this is the 'X-linked' genomic element 
	initializeSex("A"); // turn on sex and model autosomes
	rates = c(0, 0.5, 0); ends = c(99999, 100000, 199999);
	initializeRecombinationRate(rates, ends); //force random segregation of g1 and g2 to simulate diff.
}
1 late() { // initialize the pop, with a Y marker for each male 
	sim.addSubpop("p1", 500);
	i = p1.individuals;
	i[i.sex == "M"].genome2.addNewMutation(m3, 0.0, 199999);
	i[i.sex == "F"].genome1.addNewMutation(m2, 0.0, 199999);
	mut = sim.mutationsOfType(m2);
	i[i.sex == "F"].genome2.addMutations(mut);
	i[i.sex == "M"].genome1.addMutations(mut);
	
	//now add adaptive mutations
	target = sample(i[i.sex=="M"].genome1, 5);
	target.addNewDrawnMutation(m4, 10000);
	target1 = sample(i[i.sex=="M"].genome1, 5);
	target1.addNewDrawnMutation(m5, 189999); //add X-linked mutation to 5 genomes
}

//fitnessEffect callbacks
fitnessEffect(){
	dominancecoef=0.75;

	fullcomb=1.6;//fitness of inidividual if autosomal mutation is homozygous and X-linked mutation also present
	heterozygouscomb=1.3+(0.3*dominancecoef);//fitness of inidividual if autosomal mutation is heterozygous and x-linked mutation also present
	fullsingle=1.3;//fitness of inidividual if only the x-linked mutation is present, OR if only the autosomal mutation is present and homozygous
	heterozygousauto=1.0+(0.3*dominancecoef);//fitness of individual if only the autosomal mutation is present, but heterozygous
	
	if (individual.sex == "M")
		if (individual.genome1.countOfMutationsOfType(m4) & individual.genome2.countOfMutationsOfType(m4) & individual.genome1.countOfMutationsOfType(m5))
			return fullcomb;//if homozygous Cw, and Fw is on X, full summed fitness benefit
		else if (individual.countOfMutationsOfType(m4)  & individual.genome1.countOfMutationsOfType(m5)) 
			return heterozygouscomb;//if heterozygous Cw, and Fw is on X, full benefit of Fw, half benefit of Cw
		else if (individual.genome1.countOfMutationsOfType(m5))
			return fullsingle;//if Fw is on X, full Fw benefit
		else if (individual.genome1.countOfMutationsOfType(m4) & individual.genome2.countOfMutationsOfType(m4))
			return fullsingle;//if Cw is homozygous, full Cw benefit
		else if (individual.countOfMutationsOfType(m4))
			return heterozygousauto;//if Cw is heterozygous, half Cw benefit
	else
		return 1.0;	
}

modifyChild() {
	numX = sum(child.genomes.containsMarkerMutation(m2, 199999)); // num of m2 (X)
	numY = sum(child.genomes.containsMarkerMutation(m3, 199999)); // num of m3 (Y)
	numCw = sum(child.genomes.containsMarkerMutation(m4, 10000)); // num of m4 (Cw)
	numFw = sum(child.genomes.containsMarkerMutation(m5, 189999)); // num of m5 (Fw)
	numFwg2 = sum(child.genome2.containsMarkerMutation(m5, 189999)); // num of m5 (Fw) on g2
	if (child.sex == "M" & numFwg2>1 ) return F;//Fw should never be on male genome2
	
	//no male should have 2 Fw
	if (child.sex == "M" & numFw > 1) return F;
	
	if (numY > 1) stop("### ERROR: too many Ys"); // females should have 0 Y's 
	if (child.sex == "F" & numY > 0) return F; // males should have 1 Y 
	if (child.sex == "M" & numY == 0) return F;
	if (child.sex == "F" & numX < 2) return F; // females should have 2 X 
	if (child.sex == "M" & numY > 1) return F;

	return T;
}

100 late(){
sim.simulationFinished();
outputFull();
}
```
I got essentially identical results using mutationEffect callbacks for each mutation separately -- see below. This can be substituted for the fitnessEffect callback above.
```

mutationEffect(m4) {
	
	domcoeff=0.75;dominance coefficient
	fitincabs = 0.3;
	fitinccomb = 0.265;
	
	if (individual.genome1.countOfMutationsOfType(m5) & individual.sex == "M")
		 if (homozygous)// in non-additive scenario, comment out homozygous if statement
			return 1.0 + fitinccomb;//if m4 homozygous, full combined fitness benefit
		else//in non-additive scenario, ignore m4 heterozygosity here
			return 1.0 + fitinccomb*domcoeff;//if m4 is heterozygous, multiply by dominance coefficient
	else if (individual.sex == "M")
		if (homozygous)
			return 1.0 + fitincabs;//if m4 homozygous, full individual fitness benefit
		else
			return 1.0 + fitincabs*domcoeff;//if m4 is heterozygous, multiply by dominance coefficient
	else
		return 1.0;
}
mutationEffect(m5) {
	domcoeff=0.75;//dominance coefficient	
	fitincabs = 0.3;//fit benefit when m4 absent
	fitinccomb = 0.265;//fit benefit when m5 present
	if (individual.countOfMutationsOfType(m4) & individual.sex == "M")
			return 1.0 + fitinccomb;//if m4 present, combined fitness benefit
	else if (individual.sex == "M")
		return 1.0 + fitincabs;
	else
		return 1.0;
}
``
