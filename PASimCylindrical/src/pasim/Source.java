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
public class Source extends EMobject{
    protected double Rmin;
    protected double Amin;
    protected double Zmin;
    protected double Rmax;
    protected double Amax;
    protected double Zmax;
    protected int is;
    protected int js;
    protected int ks;
    protected int ie;
    protected int je;
    protected int ke;
    protected int direction;
    
    public Source(double Rmax,double Amax,double Zmax,double Rmin,double Amin,double Zmin,Cell c){
        super(1,1,1,c);
        this.Rmax = Rmax;
        this.Amax = Amax;
        this.Zmax = Zmax;
        this.Rmin = Rmin;
        this.Amin = Amin;
        this.Zmin = Zmin;        
    }
    
}
