package pasim;

public class WaveForm {
    private int cellsPerWavelength;
    private double amplitude;
    private double phaseShift;
    private double[] waveform;
    
    public WaveForm(int cellsPerWavelength, double amplitude, double phaseShift, int timeSteps){
        this.cellsPerWavelength = cellsPerWavelength;
        this.amplitude = amplitude;
        this.phaseShift = phaseShift;
        this.waveform = new double[timeSteps];
    }
    
    public void buildWave(int timeSteps){
        waveform = new double[timeSteps];
    }
    
    public int getCellsPerWavelength(){
        return this.cellsPerWavelength;
    }
    
    public void setCellsPerWavelength(int cellsPerWavelength){
        this.cellsPerWavelength = cellsPerWavelength;
    }
    
    public void set(double value, int index){
        this.waveform[index] = value;
    }
    
    public double get(int index){
        return this.waveform[index];
    }
}
