package scheduler;

import java.util.Random;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Collections;
import simulation.Simulation;
import simulation.Util;

public class APO_Scheduler extends TaskScheduler {

	private int PROTOZOA_SIZE = 40;
    private int TYPE_OF_ALGORITHM = 0;
    
	
    public APO_Scheduler(Simulation sim) {
		super(sim);
		rng = sim.getRng();
	}

	@Override
	public int[] schedule(int MAX_FES) {
		max_evaluations = MAX_FES;
		
		check_parameters();
		
		init();
		
		run_apo();
		
		int[] mapping = Util.discretizeSol(best_protozoa.getPosition());
//		for(int i=0; i<mapping.length; i++) {
//			System.out.print(mapping[i] + " ");
//		}
//		System.out.println();
		return mapping;
	}
	
	private boolean check_parameters() {
		
		dim = sim.getNumOfCloudlets();
		
		size = PROTOZOA_SIZE;
		
		population = new ArrayList<Protozoa>(size);
		newPopulation = new ArrayList<Protozoa>(size);
		epn = new ArrayList<Protozoa>(size);
		
		
		type = TYPE_OF_ALGORITHM;
		
		return true;
	}
	
	public void init() {
		eval = 0;
		for (int i = 0; i < size && eval < max_evaluations; i++) {
            population.add(new Protozoa(dim, max_bound, 0));
            newPopulation.add(new Protozoa(dim, max_bound, 0));
            epn.add(new Protozoa(dim, max_bound, 0));
            population.get(i).setFitness(sim.predictFitnessValue(Util.discretizeSol(population.get(i).getPosition())));
            eval++;
        }
	}
	
	
	public void run_apo() {        
		Comparator<Protozoa> fitnessCompare = new Comparator<Protozoa>() {
			@Override
			public int compare(Protozoa w1, Protozoa w2) {
				return Double.compare(w1.getFitness(), w2.getFitness());
			}
		};
		while(eval < max_evaluations) {
			Collections.sort(population, fitnessCompare);
			pf = pfmax * Math.random();
			int len = (int) (Math.ceil(PROTOZOA_SIZE)*pf);
			int[] Drindex = randomIntegers(PROTOZOA_SIZE,len);
			
            // Update each population member	
            for (int i = 0; i <  population.size() && eval < max_evaluations; i++) {
            	final int target = i;
            	boolean isPresent = Arrays.stream(Drindex).anyMatch(element -> element == target);
            	if(isPresent) {
            		pdr = 0.5 * (1 + Math.cos((1 - i / (double) PROTOZOA_SIZE) * Math.PI));
            		if(pdr > Math.random() && i!=0) {
            			newPopulation.get(i).setPosition(new Protozoa(dim, max_bound, 0).getPosition());             			
            		}else {
            			int Flag = Math.random() > 0.5 ? 1 : -1; // +- (plus minus)
                        int[] Mr = new int[dim];
                        int[] randomDimIndices = getRandomDimIndices(dim);
                        for (int index : randomDimIndices) {
                            Mr[index] = 1;
                        }
            			double[] oldDim = population.get(i).getPosition();
             			for (int d = 0; d < dim; d++) {
                            newPopulation.get(i).position[d] = oldDim[d] + Flag * Math.random() * (min_bound + Math.random() * (max_bound - min_bound)) * Mr[d];
                        }
            		}
            	}else {
            		double f = Math.random() * (1 + Math.cos(eval / (double) max_evaluations * Math.PI)); // foraging factor
                    int[] Mf = new int[dim]; // Mf is a mapping vector in foraging
                    int[] randomDimIndices = getRandomDimIndices(dim * i / PROTOZOA_SIZE);
                    for (int index : randomDimIndices) {
                        Mf[index] = 1;
                    }
                    double pah = 0.5 * (1 + Math.cos(eval / (double) max_evaluations * Math.PI));
                    if (Math.random() < pah) { 
                        int j = (int) (Math.random() * PROTOZOA_SIZE); 
                        for (int k = 0; k < np; k++) { 
                            int km, kp;
                            if (i == 0) {
                                km = i;
                                kp = i + (int) (Math.random() * (PROTOZOA_SIZE - i));	
                            } else if (i == PROTOZOA_SIZE - 1) {
                                km = (int) (Math.random() * PROTOZOA_SIZE);
                                kp = i;
                            } else {
                                km = (int) (Math.random() * i);
                                kp = i + (int) (Math.random() * (PROTOZOA_SIZE - i));
                            }
                            double wa = Math.exp(-Math.abs(population.get(km).getFitness() / (population.get(kp).getFitness() + 1e-9)));
                            for (int d = 0; d < dim; d++) {
                                epn.get(k).position[d] = wa * (population.get(km).getPosition()[d] - population.get(kp).getPosition()[d]);
                            }
                        }                        
                        for (int d = 0; d < dim; d++) {
                        	double sum = 0;
                            for(int s=0; s<PROTOZOA_SIZE; s++) {
                            	sum += epn.get(s).getPosition()[d];
                            }
                            newPopulation.get(i).position[d] = population.get(i).getPosition()[d] + f * (population.get(j).getPosition()[d] - 										population.get(i).getPosition()[d] + (sum/np)) * Mf[d];
                        }
                    }else { 
                        for (int k = 0; k < np; k++) { 
                            int imk, ipk;
                            if (i == 0) {
                                imk = i;
                                ipk = i + k;
                            } else if (i == PROTOZOA_SIZE - 1) {
                                imk = PROTOZOA_SIZE - k - 1;
                                ipk = i;
                            } else {
                                imk = i - k;
                                ipk = i + k;
                            }
                            if (imk < 0) {
                                imk = 0;
                            } 
                            if (ipk >= PROTOZOA_SIZE) {
                                ipk = PROTOZOA_SIZE - 1;
                            }
                            double wh = Math.exp(-Math.abs(population.get(imk).getFitness() / (population.get(ipk).getFitness() + 1e-9)));
                            for (int d = 0; d < dim; d++) {
                                epn.get(k).position[d] = wh * (population.get(imk).getPosition()[d] - population.get(ipk).getPosition()[d]);
                            }
                        }
                        int Flag = Math.random() > 0.5 ? 1 : -1; // +- (plus minus)
                        double[] Xnear = new double[dim];
                        for (int d = 0; d < dim; d++) {
                            Xnear[d] = (1 + Flag * Math.random() * (1 - eval / (double) max_evaluations)) * population.get(i).getPosition()[d];
                        }
                        for (int d = 0; d < dim; d++) {
                        	double sum = 0;
                            for(int s=0; s<PROTOZOA_SIZE; s++) {
                            	sum += epn.get(s).getPosition()[d];
                            }
                            newPopulation.get(i).position[d] = newPopulation.get(i).getPosition()[d] + f * (Xnear[d] - population.get(i).getPosition()[d] + (1.0 / np * sum)) * Mf[d];
                        }
                    }
            	}
            }
            for (int j = 0; j < PROTOZOA_SIZE; j++) {
            	for (int d = 0; d < dim; d++) {
            		newPopulation.get(j).position[d] = Math.min(Math.max(min_bound, newPopulation.get(j).getPosition()[d]), max_bound);
            	}
                newPopulation.get(j).setFitness(sim.predictFitnessValue(Util.discretizeSol(newPopulation.get(j).getPosition()))) ;
                eval++;
            }
            for (int j = 0; j < PROTOZOA_SIZE; j++) {
            	if (population.get(j).getFitness() > newPopulation.get(j).getFitness()) {
            		population.get(j).setPosition(newPopulation.get(j).getPosition());
            		population.get(j).setFitness(newPopulation.get(j).getFitness());
            	}	
            }
		} 	
		Collections.sort(population, fitnessCompare);
		best_protozoa = population.get(0);
    }
	
	public int[] randomIntegers(int upperBound, int length) {
        int[] randomArray = new int[length];
        for (int i = 0; i < length; i++) {
            // Generate a random value between the lower and upper bounds (inclusive)
            randomArray[i] = rng.nextInt(upperBound + 1);
        }
        return randomArray;
    }

    private static int[] getRandomDimIndices(int count) {
        int[] result = new int[count];
        for (int i = 0; i < count; i++) {
            result[i] = (int) (Math.random() * result.length);
        }
        return result;
    } 
	
	protected int dim;
	
	protected int size;
	
	protected int np = 3;
	
	protected double pfmax = 0.5;
	
	protected double pf;
	
	protected ArrayList<Protozoa> population;

	protected ArrayList<Protozoa> newPopulation;
	
	protected ArrayList<Protozoa> epn;
	
	protected Random rng;
	
	protected double pdr;
	
	protected int max_evaluations;
	
	protected int min_bound = 0;
	protected int max_bound = sim.getNumOfVMs();
	
	protected int type = 0;
	
	protected int eval;
    
    /**
     * Grid 
     */
    protected List<GridDim> grid;
    
    /**
     * No of grids
     */
    protected int noGrid = 3;
	
	protected Protozoa best_protozoa; 

}
