package pasim;

/**
 *
 * @author Sudantha
 */
public class Cell {
    private double deltaR;
    private double deltaA;
    private double deltaZ;
    private double deltaT;
    private double R;

    public Cell(){
        this.deltaR = 1e-3;
        this.deltaA = Math.PI/180000;
        this.deltaZ = 1e-3;
        this.deltaT = 1.0e-12;
        this.R = 0;
    }
    public Cell(double deltaR,double deltaA,double deltaZ,double R){
        this.deltaR = deltaR;
        this.deltaA = deltaA;
        this.deltaZ = deltaZ;
        this.deltaT = 1.0e-12;
        this.R = R;
    }
    
    public double getDeltaT(){
        return this.deltaT;
    }
    
    public void setDeltaT(double courantFactor,double speed,double R){
        this.deltaT = courantFactor/(speed*Math.sqrt((1/(Math.pow(deltaR, 2)))+(1/(Math.pow(R*deltaA, 2)))+(1/(Math.pow(deltaZ, 2)))));
    }
    
    public double getDeltaR(){
        return this.deltaR;
    }
    
    public double getDeltaA(){
        return this.deltaA;
    }
    
    public double getDeltaZ(){
        return this.deltaZ;
    }
    
    public double getR(){
        return this.R;
    }
    
    public void setDeltaR(double Value){
        this.deltaR = Value;
    }
    
    public void setDeltaA(double Value){
        this.deltaA = Value;
    }
        
    public void setDeltaZ(double Value){
        this.deltaZ = Value;
    }
    
    public double getMax(double R){
        return Math.max(deltaR, Math.max(R*deltaA, deltaZ));
    }
}
