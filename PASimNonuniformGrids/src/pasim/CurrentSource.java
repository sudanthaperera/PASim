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
public class CurrentSource extends Source {
    private Complex[] freqDomainValue;
    private double[] frequencies;
    private double[] waveform;
    private double resistance;
    private double magnitude;
    
    public static int currentSourceCount = 0;
    
    public CurrentSource(double Xmax,double Ymax,double Zmax,double Xmin,double Ymin,double Zmin, Cell c){
        super(Xmax,Ymax,Zmax,Xmin,Ymin,Zmin,c);
        this.resistance = 50.0;
        this.magnitude = 1.0;
    }
    
    public void SetDirection(int direction){
        this.direction = direction;
    }
    
    public void updateCurrentSourceHfiled(int timeIndex){
        //to be use
    }
}
