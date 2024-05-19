package scheduler;

import java.util.Random;
import simulation.Simulation;
import simulation.Util;

public class Wolf {
	
	protected int dim;
	
	protected double[] position;
	
	protected double fitness;
	
	protected static double max_bound;
	protected static double min_bound;
	
	public Wolf() {}
	
	public Wolf(int dim, int max, int min) {
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
	
	public void updatePosition(Simulation sim,Random rand,double a,Wolf alpha, Wolf beta, Wolf gamma) {
		double[] Xnew = new double[dim];
        double A1 = a * (2 * rand.nextDouble() - 1);
        double A2 = a * (2 * rand.nextDouble() - 1);
        double A3 = a * (2 * rand.nextDouble() - 1);
        double C1 = 2 * rand.nextDouble();
        double C2 = 2 * rand.nextDouble();
        double C3 = 2 * rand.nextDouble();
        
        for (int j = 0; j < dim; j++) {                    

            double X1 = alpha.getPosition()[j] - (A1 * Math.abs(C1 * alpha.getPosition()[j] -position[j]));
            double X2 = beta.getPosition()[j] - (A2 * Math.abs(C2 * beta.getPosition()[j] - position[j]));
            double X3 = gamma.getPosition()[j] - (A3 * Math.abs(C3 * gamma.getPosition()[j] - position[j]));

            Xnew[j] = (X1 + X2 + X3) / 3.0;
            if(Xnew[j]<min_bound) Xnew[j] = min_bound;
            if(Xnew[j]>max_bound) Xnew[j] = max_bound;
        }
        

        // Fitness calculation of new solution
        double fnew = sim.predictFitnessValue(Util.discretizeSol(Xnew)); 
//        System.out.println("fitness: " + fnew);
        // Greedy selection
        if (fnew < fitness) {
            position = Xnew;
            fitness = fnew;
        }
	}

}
