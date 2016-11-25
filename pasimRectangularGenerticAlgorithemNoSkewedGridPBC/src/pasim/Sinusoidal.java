/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pasim;

/**
 *
 * @author brah3093
 */
public class Sinusoidal extends WaveForm{
    private double Frequency;
    public Sinusoidal(int cellsPerWavelength, double amplitude, double phaseShift,int timeSteps){
        super(cellsPerWavelength, amplitude, phaseShift,timeSteps);
    }
    
    public double getFreq(){
        return this.Frequency;
    }
    
    public void setFreq(double Freq){
        this.Frequency = Freq;
    }
}
