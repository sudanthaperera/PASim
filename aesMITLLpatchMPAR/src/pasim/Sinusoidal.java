package pasim;

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
