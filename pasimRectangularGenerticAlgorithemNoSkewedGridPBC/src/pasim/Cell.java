package pasim;

/**
 *
 * @author Sudantha
 */
public class Cell {
    private double deltaX;
    private double deltaY;
    private double deltaZ;
    private double deltaT;

    public Cell(){
        this.deltaX = 1e-3;
        this.deltaY = 1e-3;
        this.deltaZ = 1e-3;        
    }
    public Cell(double deltaX,double deltaY,double deltaZ){
        this.deltaX = deltaX;
        this.deltaY = deltaY;
        this.deltaZ = deltaZ;
        this.deltaT = 1.0e-12;
    }
    
    public double getDeltaT(){
        return this.deltaT;
    }
    
    public void setDeltaT(double courantFactor,double speed){
        this.deltaT = courantFactor/(speed*Math.sqrt((1/(Math.pow(deltaX, 2)))+(1/(Math.pow(deltaY, 2)))+(1/(Math.pow(deltaZ, 2)))));
    }
    
    public double getDeltaX(){
        return this.deltaX;
    }
    
    public double getDeltaY(){
        return this.deltaY;
    }
    
    public double getDeltaZ(){
        return this.deltaZ;
    }
    
    public void setDeltaX(double Value){
        this.deltaX = Value;
    }
    
    public void setDeltaY(double Value){
        this.deltaY = Value;
    }
        
    public void setDeltaZ(double Value){
        this.deltaZ = Value;
    }
    
    public double getMax(){
        return Math.max(deltaX, Math.max(deltaY, deltaZ));
    }
}
