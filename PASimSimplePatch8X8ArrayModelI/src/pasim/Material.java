package pasim;

public class Material {
    private int index;
    private double epsR;
    private double muR;
    private double sigmaE;
    private double sigmaM;
    public static double operatingFrequency;
    
    public Material(){
        this.index=0;
        this.epsR = 1.0;
        this.muR = 1.0;
        this.sigmaE = 0;
        this.sigmaM = 1.0e-10;
    }
    
    public Material(int index, double epsR, double muR, double sigmaE, double sigmaM){
        this.index = index;
        
        if (epsR == 0) 
            this.epsR = 1.0;
        else 
            this.epsR = epsR;
        
        if (muR == 0)
            this.muR = 1.0;
        else
            this.muR = muR;
        
        if (sigmaE == 0)
            this.sigmaE = 0;
        else
            this.sigmaE = sigmaE;
        
        if (sigmaM == 0)
            this.sigmaM = 1.0e-10;
        else
            this.sigmaM = sigmaM;
    }
    public void tanDelta2Conductivity(double tanDelta, double OpeartingFreq){
        sigmaE = (2*Math.PI*OpeartingFreq)*(Constants.EPS0*epsR)*tanDelta;
    }
    
    public void index(int index){
        this.index = index;
    }
    
    public void epsR(double epsR){
        this.epsR = epsR;
    }
    
    public void muR(double muR){
        this.muR = muR;
    }
    
    public void sigmaE(double sigmaE){
        this.sigmaE = sigmaE;
    }
    
    public void sigmaM(double sigmaM){
        this.sigmaM = sigmaM;
    }    
    
    public int index(){
        return index;
    }
    
    public double epsR(){
        return epsR;
    }
    
    public double muR(){
        return muR;
    }
    
    public double sigmaE(){
        return sigmaE;
    }
    
    public double sigmaM(){
        return sigmaM;
    }
    
    public static Material Vacuum(int Index){
        Material mat = new Material();
        mat.index(Index);
        mat.epsR(1.0);
        mat.muR(1.0);
        mat.sigmaE(0);
        mat.sigmaM(1.0e-10);
        return mat;
    }
    
    public static Material Air(int Index){
        Material mat = new Material();
        mat.index(Index);
        mat.epsR(1.00058986);
        mat.muR(1.00000037);
        mat.sigmaE(5e-15);
        mat.sigmaM(1.0e-10);
        return mat;
    }
    
    public static Material PEC(int Index){
        Material mat = new Material();
        mat.index(Index);
        mat.epsR(1.0);
        mat.muR(1.0);
        mat.sigmaE(1.0e10);
        mat.sigmaM(1.0e-10);
        return mat;
    }
    
    public static Material PMC(int Index){
        Material mat = new Material();
        mat.index(Index);
        mat.epsR(1.0);
        mat.muR(1.0);
        mat.sigmaE(0);
        mat.sigmaM(1.0e10);
        return mat;
    }
    
    public static Material RogersRTduroid5880(int Index){
        Material mat = new Material();
        mat.index(Index);
        mat.epsR(2.2);
        mat.muR(1.0);
        mat.tanDelta2Conductivity(0.0009, operatingFrequency);
        mat.sigmaM(1.0e-10);
        return mat;
    }
    
    public static Material RogersRTduroid4350(int Index){
        Material mat = new Material();
        mat.index(Index);
        mat.epsR(3.48);
        mat.muR(1.0);
        mat.tanDelta2Conductivity(0.0037, operatingFrequency);
        mat.sigmaM(1.0e-10);
        return mat;
    }
    
    public static Material PerfectDielectric(int Index){
        Material mat = new Material();
        mat.index(Index);
        mat.epsR(2.2);
        mat.muR(1.0);
        mat.sigmaE(0);
        mat.sigmaM(1.0e-10);
        return mat;
    }
    
    public static Material Copper(int Index){
        Material mat = new Material();
        mat.index(Index);
        mat.epsR(1.0);
        mat.muR(0.999994);
        mat.sigmaE(5.96e7);
        mat.sigmaM(1.0e-10);
        return mat;
    }
}