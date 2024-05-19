package scheduler;

public class MOProtozoa {
	
	protected int dim;
	
	protected double[] position;
	
	protected double fitness;
	
	protected static double max_bound;
	protected static double min_bound;
	
	/** Objective values */
	protected double[] cost;
	protected double[] best_cost;
	
    /** Domination of particle*/
	protected boolean isDominated = false;
	/** sub grid index */
	protected int gridIndex;
	/** sub grid index */
	protected int[] gridSubIndex;
	
	/** crowding distance */
	protected double crowdingDistance = 0;
	
	public MOProtozoa() {}
	
	public MOProtozoa(int dim, int max, int min) {
		this.dim = dim;
		max_bound = max;
		min_bound = min;
		position = new double[dim];
		init();		
	}
	
	public void init() {
		for(int i=0; i<dim; i++) {
			position[i] = Math.random()*(max_bound - min_bound) + min_bound;
		}
	}
	
	public double getFitness() {
		return fitness;
	}
	
	public void setFitness(double fitness) {
		this.fitness = fitness;
	}
	
	public double[] getPosition() {
		return position;
	}
	
	public void setPosition(double[] position) {
		this.position = position;
	}   
	
	public double[] getCost() {
		return cost;
	}
	public void setCost(double[] cost) {
		this.cost = cost;
	}
    
    public double getCrowdingDistance() {
        return crowdingDistance;
    }
    
    public void setCrowdingDistance(double crowdingDistance) {
        this.crowdingDistance = crowdingDistance;
    }
    
    public void limitPosition() {
    	for(int i=0; i<this.getPosition().length; i++) {
    		if(this.getPosition()[i] < min_bound) {
    			this.getPosition()[i] = min_bound;
    		}else if(this.getPosition()[i] > max_bound) {
    			this.getPosition()[i] = max_bound;
    		}
    	}
    }
    
    public boolean dominates(MOProtozoa other) {
    	return (cost[0] <= other.cost[0] && cost[1] >= other.cost[1]) ||
                (cost[0] < other.cost[0] && cost[1] > other.cost[1]);
    }

}
