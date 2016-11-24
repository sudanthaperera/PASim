package pasim;

public class CosineModulatedGaussian extends Gaussian{
    private double centerFrequency;
    public CosineModulatedGaussian(int cellsPerWavelength, double amplitude, double phaseShift,int timeSteps){
        super(cellsPerWavelength, amplitude, phaseShift,timeSteps);
    }
    
    public double getFreq(){
        return this.centerFrequency;
    }
    
    public void setFreq(double centerFrequency){
        this.centerFrequency = centerFrequency;
    }
}
