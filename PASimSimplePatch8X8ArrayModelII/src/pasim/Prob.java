package pasim;

public class Prob{
    protected double x;
    protected double y;
    protected double z;
    protected double Xmin;
    protected double Ymin;
    protected double Zmin;
    protected double Xmax;
    protected double Ymax;
    protected double Zmax;
    
    public Prob(Position probPos){
        this.x = probPos.x;
        this.y = probPos.y;
        this.z = probPos.z;
        
        this.Xmin = this.x;
        this.Ymin = this.y;
        this.Zmin = this.z;
        
        this.Xmax = this.x;
        this.Ymax = this.y;
        this.Zmax = this.z;
    }
}
