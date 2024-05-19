package scheduler;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import simulation.MOSimulation;
import simulation.Util;
import java.util.Comparator;

public class NSGA_2 extends MOTaskScheduler {
    
    public NSGA_2(MOSimulation sim) {
        super(sim);
        rng = sim.getRng();
    }

    @Override
    public int[] schedule(int MAX_FES) {
        max_evaluations = MAX_FES;
                
        check_parameters();
        
        init();
        
        run_nsga2();
        
        int[] mapping = Util.discretizeSol(bestSolution.getPosition());
//		for(int ele:mapping) {
//			System.out.print(ele + " ");
//		}
//		System.out.println();
        return mapping;
    }
    

	private boolean check_parameters() {
		
		dim = sim.getNumOfCloudlets();
		
		size = 4;
		
		population = new ArrayList<Solution>(size);
		
		type = TYPE_OF_ALGORITHM;
		
		return true;
	}
    
    public void init() {
        eval = 0;
        for (int i = 0; i < size && eval < max_evaluations; i++) {
            Solution solution = new Solution(dim, max_bound, min_bound);
            solution.init(rng); // Initialize randomly within bounds
            solution.setCost(sim.predictCostValue(Util.discretizeSol(solution.getPosition())));
            population.add(solution);
            eval++;
        }
    }
    
    private void run_nsga2() {
        nObj = 2;
        
        while (eval < max_evaluations) {
            List<Solution> offspringPopulation = new ArrayList<>();
            for (int i = 0; i < size; i++) {
                Solution parent1 = selectParent();
                Solution parent2 = selectParent();
                double[] childPosition = crossover(parent1.getPosition(), parent2.getPosition());
                Solution child = new Solution(childPosition);
                mutate(child);
                child.setCost(sim.predictCostValue(Util.discretizeSol(child.getPosition())));
                offspringPopulation.add(child);
                eval++;
            }
            List<Solution> combinedPopulation = new ArrayList<>(population);
            combinedPopulation.addAll(offspringPopulation);
//            System.out.println("Population: " + combinedPopulation.size());
            crowdingDistanceAssignment(combinedPopulation);
            fronts = fastNonDominatedSort(combinedPopulation);
//            System.out.println("Population: " + combinedPopulation.size());
            population.clear();
            int index = 0;
            while (index < fronts.size() &&  population.size() + fronts.get(index).size() <= size) {
            	population.addAll(fronts.get(index));
                index++;
            }
            // Add individuals from the last front based on crowding distance until the population size is reached
            if (!fronts.isEmpty() && index < fronts.size()) {
                fronts.get(index).sort(Comparator.comparingDouble(s -> -s.crowdingDistance)); // Sort in descending order of crowding distance
                for (int i = 0; i < size - population.size(); i++) {
                    population.add(fronts.get(index).get(i));
                }
            }
        }
        bestSolution = population.stream().max(Comparator.comparingDouble(s -> s.crowdingDistance)).orElse(null);
    }
    
    private Solution selectParent() {
        int tournamentSize = 2; // You can adjust the tournament size as needed
        List<Solution> tournament = new ArrayList<>();
        for (int i = 0; i < tournamentSize; i++) {
            int randomIndex = rng.nextInt(population.size());
            tournament.add(population.get(randomIndex));
        }
        // Find the best individual in the tournament based on dominance rank and crowding distance
        Solution bestIndividual = tournament.stream()
                .min(Comparator.comparingInt(s -> s.rank))
                .orElse(null);
        return bestIndividual;
    }
    
    private double[] crossover(double[] parent1, double[] parent2) {
        double[] child = new double[parent1.length];
        double eta_c = 2.0; // Crossover distribution index, you may adjust it as needed
        double u, beta_q;
        double lb, ub; // Lower and upper bounds for variable values

        for (int i = 0; i < parent1.length; i++) {
            if (rng.nextDouble() <= 0.5) {
                if (Math.abs(parent1[i] - parent2[i]) > 1e-14) {
                    if (parent1[i] < parent2[i]) {
                        lb = parent1[i];
                        ub = parent2[i];
                    } else {
                        lb = parent2[i];
                        ub = parent1[i];
                    }
                    double rand = rng.nextDouble();
                    beta_q = 1.0 / (eta_c + 1.0);
                    double alpha = 2.0 - Math.pow((ub - lb) / (ub - lb), eta_c + 1.0);
//                    double beta = 2.0 - alpha;
                    if (rand <= 1.0 / alpha) {
                        u = Math.pow(rand * alpha, beta_q);
                    } else {
                        u = Math.pow(1.0 / (2.0 - rand * alpha), beta_q);
                    }
                    child[i] = 0.5 * ((1.0 + u) * parent1[i] + (1.0 - u) * parent2[i]);
                } else {
                    child[i] = parent1[i];
                }
            } else {
                child[i] = parent1[i];
            }
        }
        return child;
    }
    
	private void mutate(Solution solution) {
	    double eta_m = 20.0; // Mutation distribution index, you may adjust it as needed
	    double mutationProbability = 1.0 / solution.getPosition().length; // Probability of mutation per variable
	    double lb, ub; // Lower and upper bounds for variable values
	
	    for (int i = 0; i < solution.getPosition().length; i++) {
	        if (rng.nextDouble() <= mutationProbability) {
	            lb = solution.getPosition()[i] - (max_bound - min_bound) * 0.5;
	            ub = solution.getPosition()[i] + (max_bound - min_bound) * 0.5;
	            double delta1 = (solution.getPosition()[i] - lb) / (ub - lb);
	            double delta2 = (ub - solution.getPosition()[i]) / (ub - lb);
	            double rand = rng.nextDouble();
	            double mut_pow = 1.0 / (eta_m + 1.0);
	            double deltaq;
	            if (rand <= 0.5) {
	                double xy = 1.0 - delta1;
	                double val = 2.0 * rand + (1.0 - 2.0 * rand) * (Math.pow(xy, eta_m + 1.0));
	                deltaq = Math.pow(val, mut_pow) - 1.0;
	            } else {
	                double xy = 1.0 - delta2;
	                double val = 2.0 * (1.0 - rand) + 2.0 * (rand - 0.5) * (Math.pow(xy, eta_m + 1.0));
	                deltaq = 1.0 - Math.pow(val, mut_pow);
	            }
	            double x = solution.getPosition()[i] + deltaq * (ub - lb);
	            // Ensure the mutated value is within bounds
	            if (x < min_bound) {
	                x = min_bound;
	            } else if (x > max_bound) {
	                x = max_bound;
	            }
	            solution.getPosition()[i] = x;
	        }
	    }
	}

	private List<List<Solution>> fastNonDominatedSort(List<Solution> combinedPopulation) {
	    for (Solution p : combinedPopulation) {
	        p.dominationCount = 0;
	        p.dominatedSolutions = new ArrayList<>();
	        for (Solution q : combinedPopulation) {
	            if (p.dominates(q)) {
	                p.dominatedSolutions.add(q);
	            } else if (q.dominates(p)) {
	                p.dominationCount++;
	            }
	        }
	        if (p.dominationCount == 0) {
	            p.rank = 0; // Front 0
	        }
	    }
	    List<List<Solution>> newFronts = new ArrayList<>();
	    int currentRank = 0;
	    while (!combinedPopulation.isEmpty()) {
		    List<Solution> currentFront = new ArrayList<>();
	        for (Solution p : combinedPopulation) {
	            if (p.rank == currentRank) {
	                currentFront.add(p);
	            }
	        }
	        combinedPopulation.removeAll(currentFront);
	        for (Solution p : currentFront) {
	            for (Solution q : p.dominatedSolutions) {
	                q.dominationCount--;
	                if (q.dominationCount == 0) {
	                    q.rank = currentRank + 1;
	                    combinedPopulation.add(q);
	                }
	            }
	        }
	        newFronts.add(currentFront);
	        currentRank++;
	    }
	    return newFronts;
	}

    private static void crowdingDistanceAssignment(List<Solution> combinedPopulation) {
	    int nObjectives = combinedPopulation.get(0).cost.length;
	    
	    // Initialize crowding distance for each solution
	    for (Solution sol : combinedPopulation) {
	        sol.crowdingDistance = 0.0;
	    }	
	    
	    // Calculate crowding distance for each objective
	    for (int i = 0; i < nObjectives; i++) {
	        // Sort population based on the ith objective
	    	final int index = i;
	    	combinedPopulation.sort(Comparator.comparingDouble(sol -> sol.cost[index]));
	        
	        // Set boundary points' crowding distance to infinity
	    	combinedPopulation.get(0).crowdingDistance = Double.POSITIVE_INFINITY;
	    	combinedPopulation.get(combinedPopulation.size() - 1).crowdingDistance = Double.POSITIVE_INFINITY;
	        
	        // Calculate crowding distance for inner points
	        for (int j = 1; j < combinedPopulation.size() - 1; j++) {
	            Solution solPrev = combinedPopulation.get(j - 1);
	            Solution solNext = combinedPopulation.get(j + 1);
	            Solution solCurrent = combinedPopulation.get(j);
	            double distance = 0.0;
	            for (int k = 0; k < nObjectives; k++) {
	                distance += solNext.cost[k] - solPrev.cost[k];
	            }
	            solCurrent.crowdingDistance += distance;
	        }
	    }
	}
//
//    
    // Other existing methods remain unchanged
    
    // Additional fields and methods may be required for NSGA-II implementation

	// population size
	private int size;
		
	private int nObj;
		
    private int TYPE_OF_ALGORITHM = 0;    
	
	protected int dim;
	
	protected List<Solution> population;
	
	protected ArrayList<Solution> archive;
	
	protected List<List<Solution>> fronts = new ArrayList<>();
	
	protected Random rng;
	
	protected int max_evaluations;
	
	protected int min_bound = 0;
	protected int max_bound = sim.getNumOfVMs()-1;
	
	protected int type = 0;
	
	protected int eval;
	
	private static Solution bestSolution;
}
