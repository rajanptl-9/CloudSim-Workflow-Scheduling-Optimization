package scheduler;

import java.util.Random;

public class Student {
	
	protected int dim;
	
	protected double[] position;
	
	protected double fitness;
	
	protected static double max_bound;
	protected static double min_bound;
	
	public Student() {}
	
	public Student(int dim, int max, int min) {
		this.dim = dim;
		max_bound = max;
		min_bound = min;
		position = new double[dim];
		
	}
	
	public void init(Random rand) {
		for(int i=0; i<dim; i++) {
			position[i] = rand.nextDouble()*(max_bound - min_bound) + min_bound;
		}	
	}
	
	public double getFitness() {
		return fitness;
	}
	
	public double[] getPosition() {
		return position;
	}

}
