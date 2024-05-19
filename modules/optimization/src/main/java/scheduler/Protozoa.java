package scheduler;

public class Protozoa {
	
	protected int dim;
	
	protected double[] position;
	
	protected double fitness;
	
	protected static double max_bound;
	protected static double min_bound;
	
	public Protozoa() {}
	
	public Protozoa(int dim, int max, int min) {
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

}
