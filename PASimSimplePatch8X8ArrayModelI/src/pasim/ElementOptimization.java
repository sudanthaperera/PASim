/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pasim;

/**
 *
 * @author brah3093
 */
public class ElementOptimization {
    public ElementOptimization(){
        double[] fitnessGeneration = new double[10000];
        double[] fitnessBest = new double[10000];
        long start = System.currentTimeMillis();
        
        System.out.println("Starting the Active Element Synthesizer");
	
	GeneticAlgorithm ga = new GeneticAlgorithm(21, 0.005, 0.95, 2);
        
        long startInitGeneration = System.currentTimeMillis();
	
	Population population = ga.initPopulation(65, 65, 5);
	
	ga.evalPopulation(population);
        
        long totalTimeForInitGeneration = System.currentTimeMillis() - startInitGeneration;
        System.out.printf("The time taken : %d ms ---> ",totalTimeForInitGeneration);
        System.out.printf("Fitness of the initial generation  : %.12f", population.getPopulationFitness() );
        System.out.printf("\tFitness of the of the fittest solution : %.12f\n", population.getFittest(0).getFitness());
        fitnessGeneration[0] = population.getPopulationFitness();
        fitnessBest[0] = population.getFittest(0).getFitness();
        Common.save1DArray(fitnessGeneration, "FitnessGeneration");
        Common.save1DArray(fitnessBest, "FitnessBest");
        population.getFittest(0).save(0);
	
	int generation = 1;
	long startGeneration,totalTimeForGeneration;
	while (ga.isTerminationConditionMet(population) == false) {
            population.savePopulation("Population");
            startGeneration = System.currentTimeMillis();
            
            population = ga.crossoverPopulation(population);
            
            population = ga.mutatePopulation(population);
            
            ga.evalPopulation(population);
            
            totalTimeForGeneration = System.currentTimeMillis() - startGeneration;
            System.out.printf("The time taken : %d ms ---> ",totalTimeForGeneration);
            System.out.printf("Fitness of the generation #%3d : %.12f", generation, population.getPopulationFitness());
            System.out.printf("\tFitness of the of the fittest solution : %.12f\n", population.getFittest(0).getFitness());
            fitnessGeneration[generation] = population.getPopulationFitness();
            fitnessBest[generation] = population.getFittest(0).getFitness();
            Common.save1DArray(fitnessGeneration, "FitnessGeneration");
            Common.save1DArray(fitnessBest, "FitnessBest");
            population.getFittest(0).save(generation);
            
            generation++;
	}
        
        System.out.println("Found solution in " + generation + " generations");
	System.out.println("Best solution: " + population.getFittest(0).toString());
        population.getFittest(0).save(7777777);
        
        System.out.println("Done!!");
        long totalTime = System.currentTimeMillis() - start;
        System.out.printf("The time taken : %d ms\n",totalTime);
    }
}
