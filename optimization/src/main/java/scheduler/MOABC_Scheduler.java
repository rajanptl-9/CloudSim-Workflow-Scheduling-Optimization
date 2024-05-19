package scheduler;

import java.util.Random;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import simulation.MOSimulation;
import simulation.Util;


public class MOABC_Scheduler extends MOTaskScheduler {
    
	
	 public MOABC_Scheduler(MOSimulation sim) {
	        super(sim);
	        rng = sim.getRng();
	    }

	@Override
	public int[] schedule(int MAX_FES) {
		max_evaluations = MAX_FES;
		
		limit = 5;
		
		check_parameters();
		
		init();
		
		run_abc();
		
		int[] mapping = Util.discretizeSol(best_bee.getPosition());
//		System.out.println(best_bee.getFitness());
		return mapping;
	}
	
	public void init() {
		eval = 0;
		for (int i = 0; i < size && eval < max_evaluations; i++) {
            colony.add(new MOABC_Bee(dim, sim.getNumOfVMs(), 0));
            colony.get(i).init(rng);
            colony.get(i).fitness = sim.predictFitnessValue(Util.discretizeSol(colony.get(i).getPosition()));
            eval++;
        }
	}
	
	
	private boolean check_parameters() {
		
		dim = sim.getNumOfCloudlets();
		
		size = nSolutions;
		
		colony = new ArrayList<MOABC_Bee>(size);
		
		type = TYPE_OF_ALGORITHM;
		
		return true;
	}
	
	public void run_abc() {
		
//		System.out.println("Initialized bee colony");
		
		for (int i = 0; i < colony.size() && eval < max_evaluations; i++) {
            colony.get(i).cost = sim.predictCostValue(Util.discretizeSol(colony.get(i).getPosition()));
            eval++;
            colony.get(i).best_cost = colony.get(i).cost;
        }
		
		nObj = colony.get(0).cost.length;
		
//		System.out.println("Initialized bee colony cost");
		
		while (eval < max_evaluations) {
			
			// Perform employed bee phase
            for (int i = 0; i < nEmployed && eval < max_evaluations; i++) {
                // Select a solution randomly
                int solutionIndex = new Random().nextInt(nSolutions);
                MOABC_Bee solution = colony.get(solutionIndex);

                // Generate a neighbor solution by perturbing the current solution
                MOABC_Bee neighborSolution = solution;
                for (int j = 0; j < dim; j++) {
                    neighborSolution.getPosition()[j] += (2 * new Random().nextDouble() - 1) * (max_bound - min_bound);
                    neighborSolution.getPosition()[j] = Math.max(min_bound, Math.min(max_bound, neighborSolution.getPosition()[j]));
                }

                // Evaluate objective functions for the neighbor solution
                double[] neighborObjectives = sim.predictCostValue(Util.discretizeSol(neighborSolution.getPosition()));
                eval++;
                // Update the colony if the neighbor solution dominates the current solution
                if (paretoDominance(neighborObjectives, solution.cost)) {
                    colony.get(solutionIndex).position = neighborSolution.getPosition();
                    colony.get(solutionIndex).cost = neighborObjectives;
                    colony.get(solutionIndex).trials = 0;
                }else {
                	colony.get(solutionIndex).trials += 1;
                }
            }
            
         // Perform onlooker bee phase
            for (int i = 0; i < nOnlookers && eval < max_evaluations; i++) {
                // Select a solution based on roulette wheel selection
                double[] probabilities = new double[size];
                for (int j = 0; j < size; j++) {
                    probabilities[j] = new Random().nextDouble();
                }
                double sumProbabilities = Arrays.stream(probabilities).sum();
                double randomProbability = new Random().nextDouble() * sumProbabilities;
                double cumulativeProbability = 0;
                int selectedSolutionIndex = 0;
                for (int j = 0; j < size; j++) {
                    cumulativeProbability += probabilities[j];
                    if (cumulativeProbability >= randomProbability) {
                        selectedSolutionIndex = j;
                        break;
                    }
                }
                MOABC_Bee selectedSolution = colony.get(selectedSolutionIndex);

                // Generate a neighbor solution by perturbing the selected solution
                MOABC_Bee neighborSolution = selectedSolution;
                for (int j = 0; j < dim; j++) {
                    neighborSolution.getPosition()[j] += (2 * new Random().nextDouble() - 1) * (max_bound - min_bound);
                    neighborSolution.getPosition()[j] = Math.max(min_bound, Math.min(max_bound, neighborSolution.getPosition()[j]));
                }

                // Evaluate objective functions for the neighbor solution
                double[] neighborObjectives = sim.predictCostValue(Util.discretizeSol(neighborSolution.getPosition()));
                eval++;
                // Update the colony if the neighbor solution dominates the selected solutsion
                if (paretoDominance(neighborObjectives, selectedSolution.cost)) {
                    colony.get(selectedSolutionIndex).position = neighborSolution.getPosition();
                    colony.get(selectedSolutionIndex).cost = neighborObjectives;
                    colony.get(selectedSolutionIndex).trials = 0;
                }else {
                	colony.get(selectedSolutionIndex).trials += 1;
                }
            }
            
         // Perform scout bee phase (re-initialization)
            for (int i = 0; i < size && eval < max_evaluations; i++) {
                // If a solution hasn't improved for a certain number of iterations, re-initialize it
                if (colony.get(i).trials > limit) { // Re-initialize 10% of the solutions
                    colony.get(i).init(rng);
                    colony.get(i).cost = sim.predictCostValue(Util.discretizeSol(colony.get(i).getPosition()));
                    eval++;
                }
            }
        }
				
		best_bee = findBestParticleUsingCrowdingDistance(colony);		
	}
	
	private boolean paretoDominance(double[] solution1, double[] solution2) {
		boolean xDominatesY = solution1[0] <= solution2[0] && solution1[1] >= solution2[1]; // Minimize first objective, maximize second objective
        return xDominatesY && (solution1[0] < solution2[0] || solution1[1] > solution2[1]);
	}
	
	private MOABC_Bee findBestParticleUsingCrowdingDistance(List<MOABC_Bee> archive) {
        // Calculate crowding distance for each particle
        for (MOABC_Bee particle : archive) {
            particle.crowdingDistance = 0; // Initialize crowding distance
        }
        for (int objIndex = 0; objIndex < nObj; objIndex++) {
            // Sort particles based on the objective value
        	final int finalObjIndex = objIndex;
            archive.sort((p1, p2) -> Double.compare(p1.cost[finalObjIndex], p2.cost[finalObjIndex]));

            // Set boundary crowding distance for the extreme particles
            archive.get(0).crowdingDistance = Double.POSITIVE_INFINITY;
            archive.get(archive.size() - 1).crowdingDistance = Double.POSITIVE_INFINITY;

            // Calculate crowding distance for other particles
            for (int i = 1; i < archive.size() - 1; i++) {
                double previousObjValue = archive.get(i - 1).cost[objIndex];
                double nextObjValue = archive.get(i + 1).cost[objIndex];
                archive.get(i).crowdingDistance += (nextObjValue - previousObjValue);
            }
        }

        // Find the particle with the highest crowding distance (best particle)
        MOABC_Bee bestParticle = archive.get(0);
        for (MOABC_Bee particle : archive) {
            if (particle.crowdingDistance > bestParticle.crowdingDistance) {
                bestParticle = particle;
            }
        }
        return bestParticle;
    }
	
	private int nSolutions = 50;
	
	private int nEmployed = 15;
	private int nOnlookers = 15;
	
	private static int nObj;
	
	private int limit = 10;
	
    private int TYPE_OF_ALGORITHM = 0;
	
	protected int dim;
	
	protected int size;
	
	protected ArrayList<MOABC_Bee> colony;
	
	protected Random rng;
	
	protected int max_evaluations;
	
	protected int min_bound = 0;
	protected int max_bound = sim.getNumOfVMs()-1;
	
	protected int type = 0;
	
	protected int eval;
	
	protected MOABC_Bee alpha, beta, gamma, best_bee; 

}
