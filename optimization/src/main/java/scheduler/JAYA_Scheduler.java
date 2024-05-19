package scheduler;

import org.cloudbus.cloudsim.distributions.UniformDistr;
import java.util.ArrayList;
import simulation.Simulation;
import simulation.Util;


public class JAYA_Scheduler extends TaskScheduler {

	private int WOLF_SIZE = 40;
    private int TYPE_OF_ALGORITHM = 0;
    
	
	public JAYA_Scheduler(Simulation sim) {
		super(sim);
		rand = new UniformDistr(0, 1);
	} 

	@Override
	public int[] schedule(int MAX_FES) {
		max_evaluations = MAX_FES;
		
		check_parameters();
		
		init();
		
		run_jaya();
		
		int[] mapping = Util.discretizeSol(best_person.getPosition());
//		for(int i=0; i<mapping.length; i++) {
//			System.out.print(mapping[i] + " ");
//		}
//		System.out.println();
//		System.out.println(best_person.getFitness());
		return mapping;
	}
	
	public void init() {
		for (int i = 0; i < size && eval < max_evaluations; i++) {
            population.add(new Person(dim, sim.getNumOfVMs(), 0));
            population.get(i).init(rand);
        }
	}
	
	
	private boolean check_parameters() {
		
		dim = sim.getNumOfCloudlets();
		
		size = WOLF_SIZE;
		
		population = new ArrayList<Person>(size);
		
		type = TYPE_OF_ALGORITHM;
		
		return true;
	}
	
	public void run_jaya() {
		Person minObject = new Person(dim, sim.getNumOfVMs(), 0);
		Person maxObject = new Person(dim, sim.getNumOfVMs(), 0);
		
		minObject.fitness = Double.MAX_VALUE;
		maxObject.fitness = Double.MIN_VALUE;
		for (int i = 0; i < size && eval < max_evaluations; i++) {
            population.get(i).fitness = sim.predictFitnessValue(Util.discretizeSol(population.get(i).getPosition()));
            if(minObject.fitness > population.get(i).fitness) {
            	minObject = population.get(i);
            }
            if(maxObject.fitness < population.get(i).fitness) {
            	maxObject = population.get(i);
            }
            eval++;
        }
		
		while ( eval < max_evaluations ) {
			
            // Update each population member
            for (int i = 0; i <  population.size(); i++) {
            	
            	double[] Xnew = new double[dim];
                
                for (int j = 0; j < dim  && eval < max_evaluations; j++) {                    
                	double r1 = rand.sample();
                	double r2 = rand.sample();
                    
                    Xnew[j] = population.get(i).getPosition()[j] + (r1 * (minObject.getPosition()[j]-Math.abs(population.get(i).getPosition()[j]))) - (r2 * (maxObject.getPosition()[j]-Math.abs(population.get(i).getPosition()[j])));
                    if(Xnew[j]<min_bound) Xnew[j] = min_bound;
                    if(Xnew[j]>max_bound) Xnew[j] = max_bound;
                }
                

                // Fitness calculation of new solution
                double fnew = sim.predictFitnessValue(Util.discretizeSol(population.get(i).getPosition())); 
                eval++;
//                System.out.println("fitness: " + fnew);
                // Greedy selection
                if (fnew < population.get(i).getFitness()) {
                    population.get(i).position = Xnew;
                    population.get(i).fitness = fnew;
                    if(minObject.fitness > population.get(i).fitness) {
                    	minObject = population.get(i);
                    }
                    if(maxObject.fitness < population.get(i).fitness) {
                    	maxObject = population.get(i);
                    }
                }
            }
        }
		best_person = population.get(0);		
	}
	
	protected int dim;
	
	protected int size;
	
	protected ArrayList<Person> population;
	
	UniformDistr rand;
	
	protected int eval = 0;
	protected int max_evaluations;
	
	protected int min_bound = 0;
	protected int max_bound = sim.getNumOfVMs()-1;
	
	protected int type = 0;
		
	protected Person alpha, beta, gamma, best_person; 

}
