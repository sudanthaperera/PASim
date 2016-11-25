package pasim;

/**
 *
 * @author Sudantha
 */
public class Brick extends Shape3D{
    private double Rmin;
    private double Amin;
    private double Zmin;
    private double Rmax;
    private double Amax;
    private double Zmax;
    
    public Brick(double Rmin, double Amin, double Zmin, double Rmax, double Amax, double Zmax, Material material){
        super(material);
        this.Rmin = Rmin;
        this.Amin = Amin;
        this.Zmin = Zmin;
        this.Rmax = Rmax;
        this.Amax = Amax;
        this.Zmax = Zmax;
    }
    
    public double getRmin(){
        return Rmin;
    }
    
    public double getAmin(){
        return Amin;
    }

    public double getZmin(){
        return Zmin;
    }
    
    public double getRmax(){
        return Rmax;
    }
    
    public double getAmax(){
        return Amax;
    }
    
    public double getZmax(){
        return Zmax;
    }

    public int getMaterialType(){
        return this.material.index();
    }
    
    public void setRmin(double value){
        Rmin = value;
    }
       
    public void setAmin(double value){
        Amin = value;
    }
       
    public void setZmin(double value){
        Zmin = value;
    }
        
    public void setRmax(double value){
        Rmax = value;
    }
    
    public void setAmax(double value){
        Amax = value;
    }
    
    public void setZmax(double value){
        Zmax = value;
    }
    
    public Material getMaterial(){
        return this.material;
    }
}
