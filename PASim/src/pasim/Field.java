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
    protected Boundary b;
    
    public Field(int nx,int ny,int nz, Cell c){
        super(nx,ny,nz,c);
        this.dt = c.getDeltaT();
    }
    
    public void setBoundary(Boundary b){
        this.b=b;
    }

    public void setHX(int xIndex, int yIndex, int zIndex, double val){
        Hx[xIndex][yIndex][zIndex] = val;
    }

    public void setHY(int xIndex, int yIndex, int zIndex, double val){
        Hy[xIndex][yIndex][zIndex] = val;
    }
    
    public void setHZ(int xIndex, int yIndex, int zIndex, double val){
        Hz[xIndex][yIndex][zIndex] = val;
    }
    
    public double getHX(int xIndex, int yIndex, int zIndex){
        return Hx[xIndex][yIndex][zIndex];
    }

    public double getHY(int xIndex, int yIndex, int zIndex){
        return Hy[xIndex][yIndex][zIndex];
    }
    
    public double getHZ(int xIndex, int yIndex, int zIndex){
        return Hz[xIndex][yIndex][zIndex];
    }
    
    
    public void setEX(int xIndex, int yIndex, int zIndex, double val){
        Ex[xIndex][yIndex][zIndex] = val;
    }

    public void setEY(int xIndex, int yIndex, int zIndex, double val){
        Ey[xIndex][yIndex][zIndex] = val;
    }
    
    public void setEZ(int xIndex, int yIndex, int zIndex, double val){
        Ez[xIndex][yIndex][zIndex] = val;
    }
    
    public double getEX(int xIndex, int yIndex, int zIndex){
        return Ex[xIndex][yIndex][zIndex];
    }

    public double getEY(int xIndex, int yIndex, int zIndex){
        return Ey[xIndex][yIndex][zIndex];
    }
    
    public double getEZ(int xIndex, int yIndex, int zIndex){
        return Ez[xIndex][yIndex][zIndex];
    }
}
