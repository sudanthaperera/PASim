/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pasim;

/**
 *
 * @author Sudantha
 */
public class Material {
    private int index;
    private double epsR;
    private double muR;
    private double sigmaE;
    private double sigmaM;
    
    public Material(){
        
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
            this.sigmaM = 1.0e-20;
        else
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
}