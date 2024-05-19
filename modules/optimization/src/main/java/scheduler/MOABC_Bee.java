package scheduler;

import java.util.Random;

public class MOABC_Bee {
	
	/** Problem's dimension */
	protected int dimension;
	/** Actual position */
	protected double[] position;
	/** Problem bounds */
	protected static double max_bounds;
	protected static double min_bounds;
	/** Objective values */
	protected double[] cost;
	protected double[] best_cost;
	/** Fitness of actual position */
	protected double fitness;
	/** Fitness of best position */
	protected double best_fitness;
    /** Domination of particle*/
	protected boolean isDominated = false;
	/** crowding distance */
	protected double crowdingDistance = 0;
	/** trials */
	protected int trials;
	
	/** Default constructor empty */
	public MOABC_Bee(){}
	
	/** Particle constructor
	 * @param dim Problem's dimension
	 * @param max High bound for each dimension
	 * @param min Low bound for each dimension
	 */
	public MOABC_Bee(int dim, double max, double min){
		dimension = dim;
		position = new double[dimension];
		max_bounds = max;
		min_bounds = min;
		trials = 0;
	}

	/** Initializes a particle position and velocity
	 * @param rand Random generator.
	 */
	public void init(Random rand){
		for(int i = 0; i < dimension; i++){
			position[i] = rand.nextDouble()*(max_bounds - min_bounds) + min_bounds;
		}
	}
	
	public double[] getPosition() {
		return position;
	}
	
	public double[] getCost() {
		return cost;
	}
};

