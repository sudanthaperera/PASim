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
    
    public Brick(){
        super(new Material());
        this.Xmin = 0;
        this.Ymin = 0;
        this.Zmin = 0;
        this.Xmax = 0;
        this.Ymax = 0;
        this.Zmax = 0;
    }
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
    
    public void Xmin(double value){
        Xmin = value;
    }
       
    public void Ymin(double value){
        Ymin = value;
    }
       
    public void Zmin(double value){
        Zmin = value;
    }
        
    public void Xmax(double value){
        Xmax = value;
    }
    
    public void Ymax(double value){
        Ymax = value;
    }
    
    public void Zmax(double value){
        Zmax = value;
    }
    
    public void Material(Material m){
        this.material = m;
    }
    
    public Material getMaterial(){
        return this.material;
    }
}
