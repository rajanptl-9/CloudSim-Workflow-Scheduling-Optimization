package scheduler;

import java.util.ArrayList;
import simulation.Simulation;
import simulation.Util;
import java.util.Random;

public class TLBO_Scheduler extends TaskScheduler {

	private int WOLF_SIZE = 40;
    private int TYPE_OF_ALGORITHM = 0;
    
	
	public TLBO_Scheduler(Simulation sim) {
		super(sim);
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
		eval=0;
		for (int i = 0; i < size && eval < max_evaluations; i++) {
            students.add(new Student(dim, sim.getNumOfVMs(), 0));
            students.get(i).init(rand);
        }
	}	
	
	private boolean check_parameters() {
		
		dim = sim.getNumOfCloudlets();
		
		size = WOLF_SIZE;
		
		students = new ArrayList<Student>(size);
		
		type = TYPE_OF_ALGORITHM;
		
		return true;
	}
	
	public static int selectRandomInteger(int limit) {
        Random random = new Random();
        double u = random.nextDouble(); // Uniform random number between 0 and 1
        double p = -Math.log(1 - u) / 0.5; // Exponentially distributed random number
        int randomNumber = (int) Math.round(p); // Round to nearest integer

        // Ensure the selected number is within the range [0, limit]
        return Math.max(0, Math.min(limit, randomNumber));
    }
	
	public void run_jaya() {
		for (int i = 0; i < size && eval < max_evaluations; i++) {
            students.get(i).fitness = sim.predictFitnessValue(Util.discretizeSol(students.get(i).getPosition()));
            eval++;
        }
			
		teacher = students.stream().min((w1,w2) -> Double.compare(w1.getFitness(), w2.getFitness())).orElse(null);
		
		double[] Xmean = new double[dim];
		for (int i = 0; i <  dim; i++) {
			double sum = 0.0;
			for (int j = 0; j <  students.size(); j++) {
				sum = sum + students.get(j).getPosition()[i];
			}
			sum = sum / students.size();
			Xmean[i] = sum;
//			System.out.println(sum + " ");
		}
		
		while ( eval < max_evaluations ) {
			int tf = (int)Math.round(1+rand.nextDouble());
            // Update each students member
            for (int i = 0; i <  students.size() && eval < max_evaluations; i++) {
            	
            	double[] Xnewt = new double[dim];
                
                for (int j = 0; j < dim  && eval < max_evaluations; j++) {                    
                	double r = rand.nextDouble();
                    Xnewt[j] = students.get(i).getPosition()[j] + (r * (teacher.getPosition()[j]- (tf * Xmean[j])));

                    if(Xnewt[j]<min_bound) Xnewt[j] = min_bound;
                    if(Xnewt[j]>max_bound) Xnewt[j] = max_bound;
                }
                

                // Fitness calculation of new solution
                double fnewt = sim.predictFitnessValue(Util.discretizeSol(Xnewt)); 
                eval++;
                
                double[] Xp;
                int p = i;
                while(p != i) {
                	p = selectRandomInteger(24);
                }
                Xp = students.get(p).getPosition();
                double[] Xnewl = new double[dim];
                if(students.get(i).getFitness() < students.get(p).getFitness()) {
	                for (int j = 0; j < dim; j++) {
	                    double r = rand.nextDouble();
	                    Xnewl[j] = students.get(i).getPosition()[j] + ( r * ( students.get(i).getPosition()[j] - Xp[j] ) );	                    
	                    if(Xnewl[j]<min_bound) Xnewl[j] = min_bound;
	                    if(Xnewl[j]>max_bound) Xnewl[j] = max_bound;
	                }
                }else {
                	for (int j = 0; j < dim; j++) {                  		
	                    double r = rand.nextDouble();	                    
	                    Xnewl[j] = students.get(i).getPosition()[j] - ( r * ( students.get(i).getPosition()[j] - Xp[j] ) );	
	                    if(Xnewl[j]<min_bound) Xnewl[j] = min_bound;
	                    if(Xnewl[j]>max_bound) Xnewl[j] = max_bound;
	                }                	
                }
                

                // Fitness calculation of new solution
                double fnewl = sim.predictFitnessValue(Util.discretizeSol(Xnewl)); 
                eval++;
//                System.out.println("fitness: " + fnew);
                // Greedy selection
                if (fnewt < students.get(i).getFitness() && fnewt < fnewl) {
                    students.get(i).position = Xnewt;
                    students.get(i).fitness = fnewt;
                }else if(fnewl < students.get(i).getFitness()) {
                	students.get(i).position = Xnewl;
                    students.get(i).fitness = fnewl;
                }                
            }
    		
    		teacher = students.stream().min((w1,w2) -> Double.compare(w1.getFitness(), w2.getFitness())).orElse(null);
    		
    		Xmean = new double[dim];
    		for (int i = 0; i <  dim; i++) {
    			double sum = 0.0;
    			for (int j = 0; j <  students.size(); j++) {
    				sum = sum + students.get(j).getPosition()[i];
    			}
    			sum = sum / students.size();
    			Xmean[i] = sum;
//    			System.out.println(sum + " ");
    		}
        }
		best_person = teacher;
	}
	
	protected int dim;
	
	protected int size;
	
	protected ArrayList<Student> students;
	
	Random rand = new Random();
	
	protected int eval = 0;
	protected int max_evaluations;
	
	protected int min_bound = 0;
	protected int max_bound = sim.getNumOfVMs()-1;
	
	protected int type = 0;
		
	protected Student teacher, best_person; 

}
