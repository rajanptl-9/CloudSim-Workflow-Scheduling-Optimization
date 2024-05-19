package scheduler;

import java.util.Random;
import java.util.ArrayList;
import java.util.List;
import simulation.MOSimulation;
import simulation.Util;
import java.util.Arrays;


public class MOAPO_Scheduler extends MOTaskScheduler {    
	
	public MOAPO_Scheduler(MOSimulation sim) {
	    super(sim);
	    rng = sim.getRng();
	}

	@Override
	public int[] schedule(int MAX_FES) {
		max_evaluations = MAX_FES;
				
		check_parameters();
		
		init();
		
		run_moapo();
		
		int[] mapping = Util.discretizeSol(best_protozoa.getPosition());
//		for(int ele:mapping) {
//			System.out.print(ele + " ");
//		}
//		System.out.println();
		return mapping;
	}
	
	public void init() {
		eval = 0;
		for (int i = 0; i < size && eval < max_evaluations; i++) {
            population.add(new MOProtozoa(dim, max_bound, 0));
            newPopulation.add(new MOProtozoa(dim, max_bound, 0));
            epn.add(new MOProtozoa(dim, max_bound, 0));
            population.get(i).setCost(sim.predictCostValue(Util.discretizeSol(population.get(i).getPosition())));
            eval++;
        }
	}	
	
	private boolean check_parameters() {
		
		dim = sim.getNumOfCloudlets();
		
		size = 40;
		
		population = new ArrayList<MOProtozoa>(size);
		newPopulation = new ArrayList<MOProtozoa>(size);
		epn = new ArrayList<MOProtozoa>(size);
		
		type = TYPE_OF_ALGORITHM;
		
		return true;
	}
	
	public void run_moapo() {
		for (int i = 0; i < population.size() && eval < max_evaluations; i++) {
            population.get(i).cost = sim.predictCostValue(Util.discretizeSol(population.get(i).position));
            eval++;
            population.get(i).setCost(population.get(i).cost);
        }
		population = determineDomination(population);
		
		archive = new ArrayList<>();
        for (MOProtozoa particle : population) {
            if (!particle.isDominated) {
                archive.add(particle);
            }
        }
		while(eval < max_evaluations) {
			pf = pfmax * Math.random();
			int len = (int) (Math.ceil(PROTOZOA_SIZE)*pf);
			int[] Drindex = randomIntegers(PROTOZOA_SIZE,len);
			
            // Update each population member	
            for (int i = 0; i <  population.size() && eval < max_evaluations; i++) {
            	final int target = i;
            	boolean isPresent = Arrays.stream(Drindex).anyMatch(element -> element == target);
            	if(isPresent) {
            		pdr = 0.5 * (1 + Math.cos((1 - i / (double) PROTOZOA_SIZE) * Math.PI));
            		if(pdr > Math.random()) {
            			newPopulation.get(i).setPosition(new MOProtozoa(dim, max_bound, 0).getPosition());             			
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
                            Xnear[d] = (1 + Flag * Math.random() * (1 - (eval / (double) max_evaluations))) * population.get(i).getPosition()[d];
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
                newPopulation.get(j).setCost(sim.predictCostValue(Util.discretizeSol(newPopulation.get(j).getPosition()))) ;
                eval++;
            }            
            if (archive.size() > nRep) {
                int extra = archive.size() - nRep;
                for (int e = 0; e < extra; e++) {
                	archive = deleteOneRepositoryMember(archive, gamma);
                }
            }
            for (int j = 0; j < PROTOZOA_SIZE; j++) {
            	if (population.get(j).getFitness() > newPopulation.get(j).getFitness()) {
            		population.get(j).setPosition(newPopulation.get(j).getPosition());
            		population.get(j).setCost(newPopulation.get(j).getCost());
            	}	
            }
            for (MOProtozoa particle : population) {
            	archive.add(particle);
            }
            ArrayList<MOProtozoa> reposArray = new ArrayList<>(archive);
            reposArray = determineDomination(reposArray);
            archive = new ArrayList<>(reposArray);
            archive.removeIf(particle -> particle.isDominated);
		} 
		best_protozoa = findBestParticleUsingCrowdingDistance(archive);
    }
	
	protected boolean dominates(double[] x, double[] y) {
    	boolean xDominatesY = x[0] <= y[0] && x[1] >= y[1]; // Minimize first objective, maximize second objective
        return xDominatesY && (x[0] < y[0] || x[1] > y[1]);
    }	
	
	protected ArrayList<MOProtozoa> determineDomination(ArrayList<MOProtozoa> pop) {
        int popLen = pop.size();
        for (int i = 0; i < popLen; i++) {
            pop.get(i).isDominated = false;
        }
        for (int i = 0; i < popLen - 1; i++) {
            for (int j = i + 1; j < popLen; j++) {
                if (dominates(pop.get(i).cost, pop.get(j).cost)) {
                    pop.get(j).isDominated = true;
                }
                if (dominates(pop.get(j).cost, pop.get(i).cost)) {
                    pop.get(i).isDominated = true;
                }
            }
        }
        return pop;	
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
    
	private ArrayList<MOProtozoa> deleteOneRepositoryMember(ArrayList<MOProtozoa> rep, double gamma) {
        int[] gridIndices = new int[rep.size()];
        for (int i = 0; i < rep.size(); i++) {
            gridIndices[i] = rep.get(i).gridIndex;
        }
        int[] oCells = Arrays.stream(gridIndices).distinct().toArray();
        int[] N = new int[oCells.length];
        for (int k = 0; k < oCells.length; k++) {
            for (int gridIndex : gridIndices) {
                if (gridIndex == oCells[k]) {
                    N[k]++;
                }
            }
        }
        double[] p = new double[oCells.length];
        for (int i = 0; i < oCells.length; i++) {
            p[i] = Math.exp(gamma * N[i]);
        }
        double sumP = Arrays.stream(p).sum();
        for (int i = 0; i < p.length; i++) {
            p[i] /= sumP;
        }
        int sci = rouletteWheelSelection(p);
        int selectedCell = oCells[sci];
        List<MOProtozoa> selectedCellMembers = new ArrayList<>();
        for (MOProtozoa particle : rep) {
            if (particle.gridIndex == selectedCell) {
                selectedCellMembers.add(particle);
            }
        }
        int selectedMemberIndex = rng.nextInt(selectedCellMembers.size());
        rep.remove(selectedMemberIndex);
        return rep;
    }
	
	private int rouletteWheelSelection(double[] p) {
        double r = rng.nextDouble();
        double[] cumSum = new double[p.length];
        cumSum[0] = p[0];
        for (int i = 1; i < p.length; i++) {
            cumSum[i] = cumSum[i - 1] + p[i];
        }
        int count = 0;
        for (double cum : cumSum) {
            if (cum < r) {
                count++;
            }
        }
        return count;
    }
			
	private MOProtozoa findBestParticleUsingCrowdingDistance(List<MOProtozoa> archive) {
        // Calculate crowding distance for each particle
        for (MOProtozoa particle : archive) {
            particle.crowdingDistance = 0; // Initialize crowding distance
        }

        int nObj = archive.get(0).cost.length;
        for (int objIndex = 0; objIndex < nObj; objIndex++) {
            // Sort particles based on the objective value
        	final int finalObjIndex = objIndex;
//            if(objIndex == 0) {
//            	archive.sort((p1, p2) -> Double.compare(p1.cost[finalObjIndex], p2.cost[finalObjIndex]));
//            }else {
//            	archive.sort((p1, p2) -> Double.compare(-p1.cost[finalObjIndex], -p2.cost[finalObjIndex]));
//            }
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
        MOProtozoa bestParticle = archive.get(0);
        for (MOProtozoa particle : archive) {
            if (particle.crowdingDistance > bestParticle.crowdingDistance) {
                bestParticle = particle;
            }
        }
        return bestParticle;
    }

//	private static int nObj; 

	private int PROTOZOA_SIZE = 40;
    private int TYPE_OF_ALGORITHM = 0;
	
	protected ArrayList<MOProtozoa> archive;

	protected int dim;
	
	protected int size;
	
	protected int np = 3;
	
	protected double pfmax = 0.5;
	
	protected double pf;
	
	protected ArrayList<MOProtozoa> population;

	protected ArrayList<MOProtozoa> newPopulation;
	
	protected ArrayList<MOProtozoa> epn;
	
	protected int nRep = 10;
	
	protected Random rng;
	
	protected double pdr;
	
	protected int max_evaluations;
	
	protected int min_bound = 0;
	protected int max_bound = sim.getNumOfVMs();
	
	protected int type = 0;
	
	protected int eval;

    protected int gamma = 1;
	
	protected MOProtozoa best_protozoa; 

}

