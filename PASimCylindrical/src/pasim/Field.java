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
public class Field extends EMobject{
    protected double dt;
    
    public Field(int nr,int na,int nz, Cell c){
        super(nr,na,nz,c);
        this.dt = c.getDeltaT();
    }

    public static void setHR(int rIndex, int aIndex, int zIndex, double val){
        Hr[rIndex][aIndex][zIndex] = val;
    }

    public static void setHA(int rIndex, int aIndex, int zIndex, double val){
        Ha[rIndex][aIndex][zIndex] = val;
    }
    
    public static void setHZ(int rIndex, int aIndex, int zIndex, double val){
        Hz[rIndex][aIndex][zIndex] = val;
    }
    
    public static double getHR(int rIndex, int aIndex, int zIndex){
        return Hr[rIndex][aIndex][zIndex];
    }

    public static double getHA(int rIndex, int aIndex, int zIndex){
        return Ha[rIndex][aIndex][zIndex];
    }
    
    public static double getHZ(int rIndex, int aIndex, int zIndex){
        return Hz[rIndex][aIndex][zIndex];
    }
    
    
    public static void setER(int rIndex, int aIndex, int zIndex, double val){
        Er[rIndex][aIndex][zIndex] = val;
    }

    public static void setEA(int rIndex, int aIndex, int zIndex, double val){
        Ea[rIndex][aIndex][zIndex] = val;
    }
    
    public static void setEZ(int rIndex, int aIndex, int zIndex, double val){
        Ez[rIndex][aIndex][zIndex] = val;
    }
    
    public static double getER(int rIndex, int aIndex, int zIndex){
        return Er[rIndex][aIndex][zIndex];
    }

    public static double getEA(int rIndex, int aIndex, int zIndex){
        return Ea[rIndex][aIndex][zIndex];
    }
    
    public static double getEZ(int rIndex, int aIndex, int zIndex){
        return Ez[rIndex][aIndex][zIndex];
    }
}
