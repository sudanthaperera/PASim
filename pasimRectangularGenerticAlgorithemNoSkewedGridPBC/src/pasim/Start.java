package pasim;

public class Start {
    public static void main(String[] args) {
        long start = System.currentTimeMillis();        
	// Create GA object
	GeneticAlgorithm ga = new GeneticAlgorithm(8, 0.001, 0.95, 2);
	// Initialize population
	Population population = ga.initPopulation(70, 70, 35);
	// Evaluate population
	ga.evalPopulation(population);
	// Keep track of current generation
	int generation = 1;
	//Start the evolution loop
	while (ga.isTerminationConditionMet(population) == false) {
            // Print fittest individual from population
            System.out.println("Best solution: " + population.getFittest(0).toString());
            // Apply crossover
            population = ga.crossoverPopulation(population);
            // Apply mutation
            population = ga.mutatePopulation(population);
            // Evaluate population
            ga.evalPopulation(population);
            // Increment the current generation
            generation++;
	}
        
        System.out.println("Found solution in " + generation + " generations");
	System.out.println("Best solution: " + population.getFittest(0).toString());
        
        System.out.println("Done!!");
        long totalTime = System.currentTimeMillis() - start;
        System.out.printf("The time taken : %d ms\n",totalTime);
    }
}
