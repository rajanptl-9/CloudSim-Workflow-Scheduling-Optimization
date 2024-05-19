/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package scheduler;

import simulation.MOSimulation;

/**
 *
 * @author osman
 */
public abstract class MOTaskScheduler {
    protected MOSimulation sim;
    
    public MOTaskScheduler(MOSimulation sim) {
        this.sim = sim;
    }
    
    /**
     * Maps given tasks to VMs.
     * @return task-to-VM mapping array.
     */
    public abstract int[] schedule(int MAX_FES);
}
