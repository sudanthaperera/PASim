/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pasim;

/**
 *
 * @author Sudantha
 */
public class Gaussian extends WaveForm{
    private double maximumFrequency;
    private double tau;
    private double t0;
    
    public Gaussian(int cellsPerWavelength, double amplitude, double phaseShift,int timeSteps){
        super(cellsPerWavelength, amplitude, phaseShift,timeSteps);
    }
    
    public void buildWave(int timeSteps){
        super.buildWave(timeSteps);
    }
    
    public double getMaxFreq(){
        return this.maximumFrequency;
    }

    public void setMaxFreq(double val){
        this.maximumFrequency = val;
    }
    
    public void setTau(double val){
        this.tau = val;
    }
    
    public void setT0(double val){
        this.t0 = val;
    }

    public double getT0(){
        return this.t0;
    }
    public double getTau(){
        return this.tau;
    }
}
