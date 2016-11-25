package pasim;

public class GeneticAlgorithm {
    private final int populationSize;
    private final double mutationRate;
    private final double crossoverRate;
    private final int elitismCount;
    private int patchDimX, patchDimY, portPosition;
    private Chromosome chrom;

    public GeneticAlgorithm(int populationSize, double mutationRate, double crossoverRate, int elitismCount) {
        this.populationSize = populationSize;
        this.mutationRate = mutationRate;
        this.crossoverRate = crossoverRate;
        this.elitismCount = elitismCount;
    }
    
    private int[][] extractMetalPatch(int[] chromosome){
        int[][] chromosomePatch = new int[patchDimX][patchDimY];
        for(int i=0;i<patchDimX;i++){
            for(int j=0;j<patchDimY;j++){
                chromosomePatch[i][j] = chromosome[patchDimY*i+j];
            }
        }
        return chromosomePatch;
    }
    
    private int[] extractPortPosition(int[] chromosome){
        int[] chromosomePort = new int[portPosition];
        
        for(int i=0;i<portPosition;i++){
            chromosomePort[i] = chromosome[patchDimX*patchDimY+i];
        }
        
        return chromosomePort;
    }

    public Population initPopulation(int patchDimX, int patchDimY, int portPosition) {
        this.patchDimX = patchDimX;
        this.patchDimY = patchDimY;
        this.portPosition = portPosition;
        chrom.layerCount =1;
        chrom.portCount = 2;
        chrom.xDim = patchDimX;
        chrom.yDim = patchDimY;
        chrom.portDim = portPosition;
	Population population = new Population(this.populationSize, this.patchDimX*this.patchDimY+this.portPosition);
	return population;
    }

    public double calcFitness(Individual individual, OptElementInstance pas) {
	// Calculate fitness
	double fitness = pas.calcFitness();

	// Store fitness
	individual.setFitness(fitness);

	return fitness;
    }

    public void evalPopulation(Population population) {
	double populationFitness = 0;
        int coreCount = Runtime.getRuntime().availableProcessors()-1;
        int groupSize;
        Individual[] individual = population.getIndividuals();
        
        for(int countIndividual = 0 ; countIndividual <= this.populationSize; countIndividual+=coreCount){
            
            if(this.populationSize - countIndividual >= coreCount){
                groupSize = coreCount;
            }
            else{
                groupSize = this.populationSize - countIndividual;
            }
            
            OptElementInstance[] pas = new OptElementInstance[groupSize];
            
            for (int index = 0 ; index < groupSize ; index++) {
                pas[index] = new OptElementInstance(extractMetalPatch(individual[countIndividual + index].getChromosome()),extractPortPosition(individual[countIndividual + index].getChromosome()));
                pas[index].t.setPriority(Thread.MAX_PRIORITY);
            }
            
            try{
                for (int index = 0 ; index < groupSize ; index++) {
                    if(pas[index].t.isAlive()){
                        pas[index].t.join();
                    }
                }
            }
            catch (InterruptedException e){
                e.printStackTrace();
            }
            
            for (int index = 0 ; index < groupSize ; index++) {
                populationFitness += calcFitness(individual[countIndividual + index],pas[index]);
            
                if(individual[countIndividual + index].getFitness()==1){
                    pas[index].saveData();
                }
                
                pas[index].saveSParam(countIndividual+index);
            }
        }
	population.setPopulationFitness(populationFitness);
    }

    public boolean isTerminationConditionMet(Population population) {
	for (Individual individual : population.getIndividuals()) {
            if (individual.getFitness() == 1) {
		return true;
            }
	}

	return false;
    }

    public Individual selectParent(Population population) {
	// Get individuals
	Individual individuals[] = population.getIndividuals();

	// Spin roulette wheel
	double populationFitness = population.getPopulationFitness();
	double rouletteWheelPosition = Math.random() * populationFitness;

	// Find parent
	double spinWheel = 0;
	for (Individual individual : individuals) {
            spinWheel += individual.getFitness();
            if (spinWheel >= rouletteWheelPosition) {
		return individual;
            }
	}
	return individuals[population.size() - 1];
    }

    public Population crossoverPopulation(Population population) {
	// Create new population
	Population newPopulation = new Population(population.size());

	// Loop over current population by fitness
	for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual parent1 = population.getFittest(populationIndex);

            // Apply crossover to this individual?
            if (this.crossoverRate > Math.random() && populationIndex >= this.elitismCount) {
		// Initialize offspring
		Individual offspring = new Individual(parent1.getChromosomeLength());
				
		// Find second parent
		Individual parent2 = selectParent(population);

		// Loop over genome
		for (int geneIndex = 0; geneIndex < parent1.getChromosomeLength(); geneIndex++) {
                    // Use half of parent1's genes and half of parent2's genes
                    if (0.5 > Math.random()) {
			offspring.setGene(geneIndex, parent1.getGene(geneIndex));
                    } else {
                    offspring.setGene(geneIndex, parent2.getGene(geneIndex));
                    }
                }

                offspring.makeSymmetric();
                // Add offspring to new population
                newPopulation.setIndividual(populationIndex, offspring);
            } 
            else {
		// Add individual to new population without applying crossover
		newPopulation.setIndividual(populationIndex, parent1);
            }
	}

	return newPopulation;
    }

    public Population mutatePopulation(Population population) {
	// Initialize new population
	Population newPopulation = new Population(this.populationSize);

	// Loop over current population by fitness
	for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual individual = population.getFittest(populationIndex);

            // Loop over individual's genes
            for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
		// Skip mutation if this is an elite individual
		if (populationIndex > this.elitismCount) {
                    // Does this gene need mutation?
                    if (this.mutationRate > Math.random()) {
			// Get new gene
			int newGene = 1;
			if (individual.getGene(geneIndex) == 1) {
                            newGene = 0;
			}
			// Mutate gene
                        individual.setGene(geneIndex, newGene);
                    }
                }
            }
            
            individual.makeSymmetric();

            // Add individual to population
            newPopulation.setIndividual(populationIndex, individual);
        }

        // Return mutated population
        return newPopulation;
    }
}
