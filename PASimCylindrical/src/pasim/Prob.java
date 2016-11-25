/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pasim;

/**
 *
 * @author brah3093
 */
public class Prob extends EMobject{
    protected double r;
    protected double a;
    protected double z;
    protected double Rmin;
    protected double Amin;
    protected double Zmin;
    protected double Rmax;
    protected double Amax;
    protected double Zmax;
    
    public Prob(Position probPos){
        this.r = probPos.r;
        this.a = probPos.a;
        this.z = probPos.z;
        
        this.Rmin = this.r;
        this.Amin = this.a;
        this.Zmin = this.z;
        
        this.Rmax = this.r;
        this.Amax = this.a;
        this.Zmax = this.z;
    }
}
