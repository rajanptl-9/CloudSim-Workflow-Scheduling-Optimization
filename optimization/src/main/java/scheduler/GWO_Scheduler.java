package scheduler;

import java.util.Random;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Collections;
import simulation.Simulation;
import simulation.Util;


public class GWO_Scheduler extends TaskScheduler {

	private int WOLF_SIZE = 40;
    private int TYPE_OF_ALGORITHM = 0;
    
	
	public GWO_Scheduler(Simulation sim) {
		super(sim);
		rng = sim.getRng();
	}

	@Override
	public int[] schedule(int MAX_FES) {
		max_evaluations = MAX_FES;
		
		check_parameters();
		
		init();
		
		run_gwo();
		
		int[] mapping = Util.discretizeSol(best_wolf.getPosition());
//		System.out.println(best_wolf.getFitness());
		return mapping;
	}
	
	public void init() {
		eval = 0;
		for (int i = 0; i < size && eval < max_evaluations; i++) {
            wolves.add(new Wolf(dim, sim.getNumOfVMs(), 0));
            wolves.get(i).init(rng);
            wolves.get(i).fitness = sim.predictFitnessValue(Util.discretizeSol(wolves.get(i).getPosition()));
            eval++;
        }
	}
	
	
	private boolean check_parameters() {
		
		dim = sim.getNumOfCloudlets();
		
		size = WOLF_SIZE;
		
		wolves = new ArrayList<Wolf>(size);
		
		type = TYPE_OF_ALGORITHM;
		
		return true;
	}
	
	public void run_gwo() {
		Comparator<Wolf> fitnessCompare = new Comparator<Wolf>() {
			@Override
			public int compare(Wolf w1, Wolf w2) {
				return Double.compare(w1.getFitness(), w2.getFitness());
			}
		};
		
		Collections.sort(wolves, fitnessCompare);
		alpha = wolves.get(0);
		beta = wolves.get(1);
		gamma = wolves.get(2);
		for (int Iter = 0; Iter < max_evaluations  && eval < max_evaluations; Iter++) {
			
			double a = 2 * (1 - (double) Iter / max_evaluations);
			
            // Update each wolves member	
            for (int i = 0; i <  wolves.size() && eval < max_evaluations; i++) {
                // Update position using alpha, beta, and gamma
                wolves.get(i).updatePosition(sim,rng,a,alpha,beta,gamma);
                eval++;
            }
            Collections.sort(wolves, fitnessCompare);
    		alpha = wolves.get(0);
    		beta = wolves.get(1);
    		gamma = wolves.get(2);    		
        }
		best_wolf = alpha;		
	}
	
	protected int dim;
	
	protected int size;
	
	protected ArrayList<Wolf> wolves;
	
	protected Random rng;
	
	protected int max_evaluations;
	
	protected int min_bound = 0;
	protected int max_bound = sim.getNumOfVMs()-1;
	
	protected int type = 0;
	
	protected int eval;
	
	protected Wolf alpha, beta, gamma, best_wolf; 

}
