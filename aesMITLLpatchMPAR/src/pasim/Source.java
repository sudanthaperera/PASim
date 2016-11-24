package pasim;

public class Source{
    protected double Xmin;
    protected double Ymin;
    protected double Zmin;
    protected double Xmax;
    protected double Ymax;
    protected double Zmax;
    protected int is;
    protected int js;
    protected int ks;
    protected int ie;
    protected int je;
    protected int ke;
    protected int direction;
    
    protected Cell c;
    protected double dt;
    
    public Source(double Xmax,double Ymax,double Zmax,double Xmin,double Ymin,double Zmin,Cell c){
        this.c = c;
        this.dt = c.getDeltaT();
        this.Xmax = Xmax;
        this.Ymax = Ymax;
        this.Zmax = Zmax;
        this.Xmin = Xmin;
        this.Ymin = Ymin;
        this.Zmin = Zmin;        
    }
    
}
