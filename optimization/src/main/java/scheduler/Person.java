package scheduler;

import org.cloudbus.cloudsim.distributions.UniformDistr;


public class Person {
	
	protected int dim;
	
	protected double[] position;
	
	protected double fitness;
	
	protected static double max_bound;
	protected static double min_bound;
	
	public Person() {}
	
	public Person(int dim, int max, int min) {
		this.dim = dim;
		max_bound = max;
		min_bound = min;
		position = new double[dim];
		
	}
	
	public void init(UniformDistr rand) {
		for(int i=0; i<dim; i++) {
			position[i] = rand.sample()*(max_bound - min_bound) + min_bound;
		}
	}
	
	public double getFitness() {
		return fitness;
	}
	
	public double[] getPosition() {
		return position;
	}

}
