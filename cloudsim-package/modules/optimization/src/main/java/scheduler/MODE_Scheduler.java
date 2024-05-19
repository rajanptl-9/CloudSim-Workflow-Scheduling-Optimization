package scheduler;

import java.util.Random;
import java.util.ArrayList;
import java.util.List;
import simulation.MOSimulation;
import simulation.Util;
import java.util.Comparator;


public class MODE_Scheduler extends MOTaskScheduler {
    
	
	 public MODE_Scheduler(MOSimulation sim) {
	        super(sim);
	        rng = sim.getRng();
	    }

	@Override
	public int[] schedule(int MAX_FES) {
		max_evaluations = MAX_FES;
				
		check_parameters();
		
		init();
		
		run_de();
		
		int[] mapping = Util.discretizeSol(bestSolution.getPosition());
//		for(int ele:mapping) {
//			System.out.print(ele + " ");
//		}
//		System.out.println();
		return mapping;
	}
	
	public void init() {
		eval = 0;
		for (int i = 0; i < size && eval < max_evaluations; i++) {
            population.add(new Solution(dim, sim.getNumOfVMs(), 0));
            population.get(i).init(rng);
            population.get(i).setCost(sim.predictCostValue(Util.discretizeSol(population.get(i).getPosition())));
            eval++;
        }
	}
	
	
	private boolean check_parameters() {
		
		dim = sim.getNumOfCloudlets();
		
		size = 40;
		
		population = new ArrayList<Solution>(size);
		
		type = TYPE_OF_ALGORITHM;
		
		return true;
	}
	
	public void run_de() {
		
//		System.out.println("Nos of VMs: " + sim.getNumOfVMs());
		
		nObj = 2;
		
//		System.out.println("Nos of Objectives: " + nObj);
		
//		System.out.println("Initialized bee population cost");
		
		while (eval < max_evaluations) {
	            // Update population using mutation and crossover
	        for (int i = 0; i < size; i++) {
	            double[] mutantVector = mutation(population, i);
	            double[] trialVector = crossover(population.get(i).getPosition(), mutantVector);

	            Solution trialSolution = new Solution(trialVector);
	            trialSolution.limitPosition();
	            trialSolution.setCost(sim.predictCostValue(Util.discretizeSol(trialSolution.getPosition())));
	            
	            // Replace target solution with trial solution if it dominates
	            if (trialSolution.dominates(population.get(i))) {
	                population.set(i, trialSolution);
	            }
	        }

	            // Update population based on crowding distance
	        updatePopulation(population);
	        eval++;
	   }		
		List<Solution> paretoOptimalSolutions = new ArrayList<>();
        for (Solution solution : population) {
            if (!isDominated(solution, population)) {
                paretoOptimalSolutions.add(solution);
            }
        }
        bestSolution = paretoOptimalSolutions.stream()
                .max(Comparator.comparingDouble(s -> s.crowdingDistance))
                .orElse(null);
	}
    
	private double[] crossover(double[] targetVector, double[] mutantVector) {
        double[] trialVector = new double[dim];
        int jRand = rng.nextInt(dim);
        for (int j = 0; j < dim; j++) {
            if (rng.nextDouble() <= CROSSOVER_PROBABILITY || j == jRand) {
                trialVector[j] = mutantVector[j];
            } else {
                trialVector[j] = targetVector[j];
            }
        }
        return trialVector;
    }
  
	private double[] mutation(List<Solution> population, int targetIndex) {
        int r1, r2, r3;
        do {
            r1 = rng.nextInt(size);
        } while (r1 != targetIndex);
        do {
            r2 = rng.nextInt(size);
        } while (r2 != targetIndex || r2 != r1);
        do {
            r3 = rng.nextInt(size);
        } while (r3 != targetIndex || r3 != r1 || r3 != r2);

        double[] mutantVector = new double[dim];
        for (int i = 0; i < dim; i++) {
            mutantVector[i] = population.get(r1).getPosition()[i]
                    + MUTATION_PROBABILITY * (population.get(r2).getPosition()[i] - population.get(r3).getPosition()[i]);
            if(mutantVector[i] < min_bound) {
            	mutantVector[i] = min_bound;
            }else if(mutantVector[i] > max_bound) {
            	mutantVector[i] = max_bound;
            }
        }
        return mutantVector;
    }
	    
	private static boolean isDominated(Solution solution, List<Solution> population) {
        for (Solution other : population) {
            if (!solution.equals(other) && other.dominates(solution)) {
                return true;
            }
        }
        return false;
    }
	
	// Update population using crowding distance
    private void updatePopulation(List<Solution> population) {
        // Evaluate objective values
        for (Solution solution : population) {
            solution.setCost(sim.predictCostValue(Util.discretizeSol(solution.getPosition())));
            solution.limitPosition();
            eval++;
        }

        // Calculate crowding distance for solutions in the population
        for (Solution solution : population) {
            solution.crowdingDistance = 0;
        }

        // Calculate crowding distance for each objective
        for (int objIndex = 0; objIndex < nObj; objIndex++) {
            final int index = objIndex;
            population.sort((s1, s2) -> Double.compare(s1.cost[index], s2.cost[index]));

            // Set extreme solutions' crowding distance to infinity
            population.get(0).crowdingDistance = Double.POSITIVE_INFINITY;
            population.get(size - 1).crowdingDistance = Double.POSITIVE_INFINITY;

            // Calculate crowding distance for the rest of the solutions
            for (int i = 1; i < size - 1; i++) {
                population.get(i).crowdingDistance +=
                        (population.get(i + 1).cost[objIndex] - population.get(i - 1).cost[objIndex]);
            }
        }
    }
	
	// population size
	private int size;
		
	private static int nObj;
		
    private int TYPE_OF_ALGORITHM = 0;    

    // Probability of crossover
    private static final double CROSSOVER_PROBABILITY = 0.8;
    
    // Probability of mutation
    private static final double MUTATION_PROBABILITY = 0.5;
	
	protected int dim;
	
	protected List<Solution> population;
	
	protected ArrayList<Solution> archive;
	
	protected Random rng;
	
	protected int max_evaluations;
	
	protected int min_bound = 0;
	protected int max_bound = sim.getNumOfVMs()-1;
	
	protected int type = 0;
	
	protected int eval;
	
	private static Solution bestSolution;

}
