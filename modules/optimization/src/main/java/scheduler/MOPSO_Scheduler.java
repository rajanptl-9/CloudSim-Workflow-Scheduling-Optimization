/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package scheduler;

import simulation.MOSimulation;
import simulation.Util;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Random;
/**
 *
 * @author osman From:https://sci2s.ugr.es/EAMHCO (Several advanced Particle
 * Swarm Optimization algorithms)
 */
public class MOPSO_Scheduler extends MOTaskScheduler {
    
    private int SWARM_SIZE = 40;
    private int TYPE_OF_ALGORITHM = 0;
    private int NEIGHBOR = SWARM_SIZE / 2; //[0, SWARM_SIZE]

    public MOPSO_Scheduler(MOSimulation sim) {
        super(sim);
        rng = sim.getRng();
    }

    @Override
    public int[] schedule(int MAX_FES) {
    	
        max_evaluations = MAX_FES;
        
        check_parameters();
        
        init();

        run_mopso();

        int[] mapping = Util.discretizeSol(best_particle.best_position);
        
//        for(int i=0; i<mapping.length; i++) {
//			System.out.print(mapping[i] + " ");
//		}
//		System.out.println(best_particle.cost[0] + " " + best_particle.cost[1]);
        
        return mapping;
    }

    /**
     * Swarm of particles
     */
    protected MOPSO_Particle[] swarm;
    /**
     * List of Archive
     */
    protected List<MOPSO_Particle> archive;
    /**
     * Number of evaluations spent
     */
    protected int eval;
    /**
     * Maximum evaluations
     */
    protected int max_evaluations;
    /**
     * Problem's dimension
     */
    protected int dim;

    /**
     * Swarm size
     */
    protected int nPop;
    /**
     * size of archive
     */
    int nRep = 8;
    /**
     * Best particle
     */
//    protected MOPSO_Particle best_particle;
    /**
     * Random generator
     */
    protected static Random rng = new Random();
    /**
     * Type of neighborhood
     */
    protected int type = 0;
    /**
     * Size of neighborhood
     */
    protected int neighbor;
    /**
     * Test function
     */
    /**
     * Seed of random generator
     */
    protected long seed;
    
    /** 
     * Ojective dimension
     */
    protected int nObj;
    
    /**
     * Grid 
     */
    protected List<GridDim> grid;
    
    /**
     * No of grids
     */
    protected int noGrid = 3;
    
    /**
     * best particle
     */
    protected MOPSO_Particle best_particle;
    
    /**
     * Some variables
     */
    protected double alpha = 0.1;
    protected int beta = 1;
    protected int gamma = 1;
//    double wDamping = 0.99;

    /**
     * Calculates the distance between two particles
     */
    private boolean check_parameters() {
        
        //Dimension
        dim = sim.getNumOfCloudlets();

        nPop = SWARM_SIZE;
        
        swarm = new MOPSO_Particle[nPop];
        
        //type of algorithm
        // 0 -> global
        // 1 -> social
        // 2 -> geographical
        type = TYPE_OF_ALGORITHM;
        
        if (type > 0) {
            neighbor = NEIGHBOR;
        }

        return true;
    }

    /**
     * Initializes the swarm
     */
    private void init() {
        eval = 0;

        for (int i = 0; i < nPop; i++) {
            swarm[i] = new MOPSO_Particle(dim, sim.getNumOfVMs()-1, 0);
            swarm[i].init(rng);
        }
    }
    
    /**
     * PSO Main function
     */
    public void run_mopso() {
//    	MOPSO_Archive archive = new MOPSO_Archive(nRep);
        /* Evaluates the initial population and find the first best*/
        
        for (int i = 0; i < swarm.length && eval < max_evaluations; i++) {
            swarm[i].cost = sim.predictCostValue(Util.discretizeSol(swarm[i].position));
            eval++;
            swarm[i].best_cost = swarm[i].cost;
            swarm[i].best_position = swarm[i].position.clone();
        }
        // finding the non dominating swarms
        swarm = determineDomination(swarm);
        
        // adding non dominating swarms to archive
        archive = new ArrayList<>();
        for (MOPSO_Particle particle : swarm) {
            if (!particle.isDominated) {
                archive.add(particle);
            }
        }
        
        nObj = archive.get(0).cost.length;
        
        grid = createGrid(archive,noGrid,alpha,nObj);        
        
        for (MOPSO_Particle repo : archive) {
            repo = findGridIndex(repo, grid);
        }
        
        PSO_Particle.update_intertia(eval, max_evaluations);
        
        /* PSO main boucle*/
        while(eval < max_evaluations) {
            for (int i = 0; i < nPop && eval < max_evaluations; i++) {
            	MOPSO_Particle leader = selectLeader(archive, beta);
                for (int j = 0; j < dim; j++) {
                    swarm[i].update_velocity(leader, rng);
                    swarm[i].update_position();
                }
                swarm[i].cost = sim.predictCostValue(Util.discretizeSol(swarm[i].position));
                eval++;
                if (dominates(swarm[i].cost, swarm[i].best_cost)) {
                    swarm[i].best_position = Arrays.copyOf(swarm[i].position, dim);
                    swarm[i].best_cost = Arrays.copyOf(swarm[i].cost, swarm[i].cost.length);
                } else {
                    if (rng.nextDouble() > 0.5) {
                        swarm[i].best_position = Arrays.copyOf(swarm[i].position, dim);
                        swarm[i].best_cost = Arrays.copyOf(swarm[i].cost, swarm[i].cost.length);
                    }
                }
//                swarm[i].mutate_particle();
            }

            for (MOPSO_Particle particle : swarm) {
            	archive.add(particle);
            }
            // repos = determineDomination(repos.toArray(new Particle[0]));
            MOPSO_Particle[] reposArray = archive.toArray(new MOPSO_Particle[0]);
            reposArray = determineDomination(reposArray);
            archive = new ArrayList<>(Arrays.asList(reposArray));
            archive.removeIf(particle -> particle.isDominated);
            
            grid = createGrid(archive, noGrid, alpha, nObj);
            
            for (MOPSO_Particle repo : archive) {
                repo = findGridIndex(repo, grid);
            }

            if (archive.size() > nRep) {
                int extra = archive.size() - nRep;
                for (int e = 0; e < extra; e++) {
                	archive = deleteOneRepositoryMember(archive, gamma);
                }
            }
            PSO_Particle.update_intertia(eval, max_evaluations);
        }
        best_particle = findBestParticleUsingCrowdingDistance(archive);
        
    }
    
    
    static boolean dominates(double[] x, double[] y) {
    	boolean xDominatesY = x[0] <= y[0] && x[1] >= y[1]; // Minimize first objective, maximize second objective
        return xDominatesY && (x[0] < y[0] || x[1] > y[1]);
    }	
    
    static MOPSO_Particle[] determineDomination(MOPSO_Particle[] pop) {
        int popLen = pop.length;
        for (int i = 0; i < popLen; i++) {
            pop[i].isDominated = false;
        }
        for (int i = 0; i < popLen - 1; i++) {
            for (int j = i + 1; j < popLen; j++) {
                if (dominates(pop[i].cost, pop[j].cost)) {
                    pop[j].isDominated = true;
                }
                if (dominates(pop[j].cost, pop[i].cost)) {
                    pop[i].isDominated = true;
                }
            }
        }
        return pop;	
    }
    
    static List<GridDim> createGrid(List<MOPSO_Particle> pop, int nGrid, double alpha, int nobj) {
        double[][] costs = pop.stream().map(p -> p.cost).toArray(double[][]::new);
        double[] cmin = Arrays.stream(costs).min((x, y) -> {
            for (int i = 0; i < x.length; i++) {
                if (x[i] != y[i]) {
                    return Double.compare(x[i], y[i]);
                }
            }
            return 0;
        }).orElse(new double[0]);

        double[] cmax = Arrays.stream(costs).max((x, y) -> {
            for (int i = 0; i < x.length; i++) {
                if (x[i] != y[i]) {
                    return Double.compare(x[i], y[i]);
                }
            }
            return 0;
        }).orElse(new double[0]);

        double[] deltaC = new double[cmin.length];
        for (int i = 0; i < cmin.length; i++) {
            deltaC[i] = cmax[i] - cmin[i];
            cmin[i] -= alpha * deltaC[i];
            cmax[i] += alpha * deltaC[i];
        }

        List<GridDim> grid = new ArrayList<>();
        for (int i = 0; i < nobj; i++) {
            GridDim gridDim = new GridDim();
            double[] dimValues = new double[nGrid + 1];
            for (int j = 0; j <= nGrid; j++) {
                dimValues[j] = cmin[i] + j * (cmax[i] - cmin[i]) / nGrid;
            }
            for (double dimValue : dimValues) {
                gridDim.lowerBounds.add(dimValue != Double.NEGATIVE_INFINITY ? dimValue : null);
            }
            for (double dimValue : dimValues) {
                gridDim.upperBounds.add(dimValue != Double.POSITIVE_INFINITY ? dimValue : null);
            }
            grid.add(gridDim);
        }
        return grid;
    }

    
    public double[] linspace(double Cmin, double Cmax, int noGrid)  {
    	double[] linspace = new double[noGrid];
        double step = (Cmax - Cmin) / (noGrid - 1);
        for (int i = 0; i < noGrid; i++) {
            linspace[i] = Cmin + step * i;
        }
        return linspace;
    }
    
    
    static MOPSO_Particle findGridIndex(MOPSO_Particle particle, List<GridDim> grid) {
        int nObj = particle.cost.length;
        particle.gridSubIndex = new int[nObj];
        for (int j = 0; j < nObj; j++) {
            int indexInDim = 0;
            for (Double upperBound : grid.get(j).upperBounds) {
                if (particle.cost[j] > upperBound) {
                    indexInDim++;
                } else {
                    break;
                }
            }
            particle.gridSubIndex[j] = indexInDim;
        }
        particle.gridIndex = particle.gridSubIndex[0];
        for (int j = 1; j < nObj; j++) {
            particle.gridIndex = particle.gridIndex * grid.size() + particle.gridSubIndex[j];
        }
        return particle;
    }
    
    
    static MOPSO_Particle selectLeader(List<MOPSO_Particle> rep, double beta) {
        double[] p = new double[rep.size()];
        int gridIndices[] = new int[rep.size()];
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
        for (int i = 0; i < oCells.length; i++) {
            p[i] = Math.exp(-beta * N[i]);
        }
        double sumP = Arrays.stream(p).sum();
        for (int i = 0; i < p.length; i++) {
            p[i] /= sumP;
        }
        int sci = rouletteWheelSelection(p);
        int selectedCell = oCells[sci];
        List<Integer> selectedCellMembers = new ArrayList<>();
        for (int gridIndex : gridIndices) {
            if (gridIndex == selectedCell) {
                selectedCellMembers.add(gridIndex);
            }
        }
        int selectedMemberIndex = rng.nextInt(selectedCellMembers.size());
        return rep.get(selectedMemberIndex);
    }

    
    static List<MOPSO_Particle> deleteOneRepositoryMember(List<MOPSO_Particle> rep, double gamma) {
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
        List<MOPSO_Particle> selectedCellMembers = new ArrayList<>();
        for (MOPSO_Particle particle : rep) {
            if (particle.gridIndex == selectedCell) {
                selectedCellMembers.add(particle);
            }
        }
        int selectedMemberIndex = rng.nextInt(selectedCellMembers.size());
        rep.remove(selectedMemberIndex);
        return rep;
    }

    
    static int rouletteWheelSelection(double[] p) {
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

    static MOPSO_Particle findBestParticleUsingCrowdingDistance(List<MOPSO_Particle> archive) {
        // Calculate crowding distance for each particle
        for (MOPSO_Particle particle : archive) {
            particle.crowdingDistance = 0; // Initialize crowding distance
        }

        int nObj = archive.get(0).cost.length;
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
        MOPSO_Particle bestParticle = archive.get(0);
        for (MOPSO_Particle particle : archive) {
            if (particle.crowdingDistance > bestParticle.crowdingDistance) {
                bestParticle = particle;
            }
        }
        return bestParticle;
    }
}
