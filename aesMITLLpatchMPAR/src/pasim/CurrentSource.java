package pasim;

public class CurrentSource extends Source {
    private Complex[] freqDomainValue;
    private double[] frequencies;
    private double[] waveform;
    private double resistance;
    private double magnitude;
    
    public CurrentSource(double Xmax,double Ymax,double Zmax,double Xmin,double Ymin,double Zmin, Cell c){
        super(Xmax,Ymax,Zmax,Xmin,Ymin,Zmin,c);
        this.resistance = 50.0;
        this.magnitude = 1.0;
    }
    
    public void SetDirection(int direction){
        this.direction = direction;
    }
    
    public void updateCurrentSourceHfiled(int timeIndex, HField H){
        //to be use
    }
}
