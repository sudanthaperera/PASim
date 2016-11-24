package pasim;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;

public class Population {
    private Individual population[];
    private double populationFitness = -1;

    public Population(int populationSize) {
	this.population = new Individual[populationSize];
    }
    
    public void savePopulation(String fileName){
        double[][] saveState = new double[population.length][population[0].getChromosome().length];
        for(int i=0;i<population.length;i++){
            for(int j=0;j<population[0].getChromosome().length;j++){
                saveState[i][j] = population[i].getChromosome()[j];
            }
        }
        Common.save2DArray(saveState, fileName);
    }
    
    public void save(){
        double[][] fitnessArray = new double[this.population.length][this.population[0].getChromosomeLength()+1];
        for(int index1 = 0; index1<this.population.length;index1++){
            for(int index2=0; index2<this.population[index1].getChromosomeLength();index2++){
                fitnessArray[index1][index2] = this.population[index1].getGene(index2);
            }
            fitnessArray[index1][this.population[index1].getChromosomeLength()] = this.population[index1].getFitness();
        }
        Common.save2DArray(fitnessArray, "FitnessArray");
    }

    public Population(int populationSize, int chromosomeLength) {
	// Initialize the population as an array of individuals
	this.population = new Individual[populationSize];

        /*
	// Create each individual in turn
	for (int individualCount = 1; individualCount < populationSize; individualCount++) {
            // Create an individual, initializing its chromosome to the given
            // length
            Individual individual = new Individual(chromosomeLength);
            // Add individual to population
            this.population[individualCount] = individual;
	}
                */
        
        //Known good one
        int[] port = new int[5];
        port[0] = 1;
        port[1] = 0;
        port[2] = 1;
        port[3] = 1;
        port[4] = 0;
        
        for (int individualCount = 0; individualCount < populationSize; individualCount++) {
            this.population[individualCount] = new Individual(chromosomeLength,Common.genInt1DArray(65*65, 1),port);
        }
        /*
        port[0] = 0;
        port[1] = 1;
        port[2] = 1;
        port[3] = 1;
        port[4] = 0;
        
        this.population[1] = new Individual(chromosomeLength,Common.genInt1DArray(65*65, 1),port);
                */
    }

    public Individual[] getIndividuals() {
	return this.population;
    }

    public Individual getFittest(int offset) {
	// Order population by fitness
	Arrays.sort(this.population, new Comparator<Individual>() {
            @Override
            public int compare(Individual o1, Individual o2) {
                if (o1.getFitness() > o2.getFitness()) {
                    return -1;
                }
                else if (o1.getFitness() < o2.getFitness()){
                    return 1;
                }
                return 0;
            }
        });

	// Return the fittest individual
	return this.population[offset];
    }

    public void setPopulationFitness(double fitness) {
	this.populationFitness = fitness;
    }

    public double getPopulationFitness() {
	return this.populationFitness;
    }

    public int size() {
	return this.population.length;
    }

    public Individual setIndividual(int offset, Individual individual) {
	return population[offset] = individual;
    }

    public Individual getIndividual(int offset) {
	return population[offset];
    }
	
    public void shuffle() {
	Random rnd = new Random();
	for (int i = population.length - 1; i > 0; i--) {
            int index = rnd.nextInt(i + 1);
            Individual a = population[index];
            population[index] = population[i];
            population[i] = a;
	}
    }
}