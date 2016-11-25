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
    
    public Source(double Xmax,double Ymax,double Zmax,double Xmin,double Ymin,double Zmin,Cell c){
        super(1,1,1,c);
        this.Xmax = Xmax;
        this.Ymax = Ymax;
        this.Zmax = Zmax;
        this.Xmin = Xmin;
        this.Ymin = Ymin;
        this.Zmin = Zmin;        
    }
    
}
