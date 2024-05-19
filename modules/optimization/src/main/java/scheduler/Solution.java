package scheduler;

import java.util.Random;
import java.util.List;

public class Solution {

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
	
	protected int rank;
	
	protected int dominationCount;
	
	public List<Solution> dominatedSolutions;
	
	// constructors
	public Solution(){}
	public Solution(double[] position) {
        this.position = position;
         // Initialize cost array
    }
	public Solution(int dim, double max, double min){
		dimension = dim;
		position = new double[dimension];
		max_bounds = max;
		min_bounds = min;
		trials = 0;
	}   
    
    public double getCrowdingDistance() {
        return crowdingDistance;
    }
    
    public void setCrowdingDistance(double crowdingDistance) {
        this.crowdingDistance = crowdingDistance;
    }
    
    public boolean dominates(Solution other) {
    	return (cost[0] <= other.cost[0] && cost[1] >= other.cost[1]) ||
                (cost[0] < other.cost[0] && cost[1] > other.cost[1]);
    }
    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Decision Variables: [");
        for (int i = 0; i < position.length; i++) {
            sb.append(position[i]);
            if (i < position.length - 1) {
                sb.append(", ");
            }
        }
        sb.append("], Objectives: [");
        for (int i = 0; i < cost.length; i++) {
            sb.append(cost[i]);
            if (i < cost.length - 1) {
                sb.append(", ");
            }
        }
        sb.append("], Crowding Distance: ").append(crowdingDistance);
        return sb.toString();
    }

    public boolean isBetterThan(Solution other) {
        return (cost[0] <= other.getCost()[0] && cost[1] >= other.getCost()[1]) || (cost[0] < other.getCost()[0] && cost[1] > other.getCost()[1]);
    }
    	
    	
    public void init(Random rand){
 		for(int i = 0; i < dimension; i++){
    		position[i] = rand.nextDouble()*(max_bounds - min_bounds) + min_bounds;
    	}
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
    
    public void limitPosition() {
    	for(int i=0; i<this.getPosition().length; i++) {
    		if(this.getPosition()[i] < min_bounds) {
    			this.getPosition()[i] = min_bounds;
    		}else if(this.getPosition()[i] > max_bounds) {
    			this.getPosition()[i] = max_bounds;
    		}
    	}
    }
}

 