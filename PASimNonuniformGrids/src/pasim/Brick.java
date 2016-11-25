package pasim;

/**
 *
 * @author Sudantha
 */
public class Brick extends Shape3D{
    private double Xmin;
    private double Ymin;
    private double Zmin;
    private double Xmax;
    private double Ymax;
    private double Zmax;
    
    public Brick(double Xmin, double Ymin, double Zmin, double Xmax, double Ymax, double Zmax, Material material){
        super(material);
        this.Xmin = Xmin;
        this.Ymin = Ymin;
        this.Zmin = Zmin;
        this.Xmax = Xmax;
        this.Ymax = Ymax;
        this.Zmax = Zmax;
    }
    
    public double getXmin(){
        return Xmin;
    }
    
    public double getYmin(){
        return Ymin;
    }

    public double getZmin(){
        return Zmin;
    }
    
    public double getXmax(){
        return Xmax;
    }
    
    public double getYmax(){
        return Ymax;
    }
    
    public double getZmax(){
        return Zmax;
    }

    public int getMaterialType(){
        return this.material.index();
    }
    
    public void setXmin(double value){
        Xmin = value;
    }
       
    public void setYmin(double value){
        Ymin = value;
    }
       
    public void setZmin(double value){
        Zmin = value;
    }
        
    public void setXmax(double value){
        Xmax = value;
    }
    
    public void setYmax(double value){
        Ymax = value;
    }
    
    public void setZmax(double value){
        Zmax = value;
    }
    
    public Material getMaterial(){
        return this.material;
    }
}
