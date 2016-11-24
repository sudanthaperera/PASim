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
public class EMobject {
    protected int nx, ny, nz;
    protected double[][][] Hx,Hy,Hz;
    protected double[][][] Ex,Ey,Ez;
    protected double[][][] Cexe,Cexhz,Cexhy,Ceye,Ceyhx,Ceyhz,Ceze,Cezhy,Cezhx;
    protected double[][][] Chxh,Chxez,Chxey,Chyh,Chyex,Chyez,Chzh,Chzey,Chzex;
    protected Cell c;
    
        
    public void saveAllCoefArrays(){
        Common.save3DArray(Cexe,"Cexe");
        Common.save3DArray(Cexhz,"Cexhz");
        Common.save3DArray(Cexhy,"Cexhy");
        Common.save3DArray(Ceye,"Ceye");
        Common.save3DArray(Ceyhx,"Ceyhx");
        Common.save3DArray(Ceyhz,"Ceyhz");
        Common.save3DArray(Ceze,"Ceze");
        Common.save3DArray(Cezhy,"Cezhy");
        Common.save3DArray(Cezhx,"Cezhx");
        
        Common.save3DArray(Chxh,"Chxh");
        Common.save3DArray(Chxez,"Chxez");
        Common.save3DArray(Chxey,"Chxey");
        Common.save3DArray(Chyh,"Chyh");
        Common.save3DArray(Chyex,"Chyex");
        Common.save3DArray(Chyez,"Chyez");
        Common.save3DArray(Chzh,"Chzh");
        Common.save3DArray(Chzey,"Chzey");
        Common.save3DArray(Chzex,"Chzex");
    }
    
    public EMobject(int nx,int ny,int nz,Cell c){
        this.nx = nx;
        this.ny = ny;
        this.nz = nz;
        this.c = c;
    }
    
    public EMobject(){
        
    }
   
    public double getCexe(int x,int y, int z){
        return this.Cexe[x][y][z];
    }
    
    public double getCexhz(int x,int y, int z){
        return this.Cexhz[x][y][z];
    }
    
    public double getCexhy(int x,int y, int z){
        return this.Cexhy[x][y][z];
    }
    
    public double getCeye(int x,int y, int z){
        return this.Ceye[x][y][z];
    }
    
    public double getCeyhx(int x,int y, int z){
        return this.Ceyhx[x][y][z];
    }
    
    public double getCeyhz(int x,int y, int z){
        return this.Ceyhz[x][y][z];
    }
    
    public double getCeze(int x,int y, int z){
        return this.Ceze[x][y][z];
    }
    
    public double getCezhy(int x,int y, int z){
        return this.Cezhy[x][y][z];
    }
    
    public double getCezhx(int x,int y, int z){
        return this.Cezhx[x][y][z];
    }
    
    public double getChxh(int x,int y, int z){
        return this.Chxh[x][y][z];
    }
    
    public double getChxez(int x,int y, int z){
        return this.Chxez[x][y][z];
    }
    
    public double getChxey(int x,int y, int z){
        return this.Chxey[x][y][z];
    }
    
    public double getChyh(int x,int y, int z){
        return this.Chyh[x][y][z];
    }
    
    public double getChyex(int x,int y, int z){
        return this.Chyex[x][y][z];
    }
    
    public double getChyez(int x,int y, int z){
        return this.Chyez[x][y][z];
    }
    
    public double getChzh(int x,int y, int z){
        return this.Chzh[x][y][z];
    }
    
    public double getChzey(int x,int y, int z){
        return this.Chzey[x][y][z];
    }
    
    public double getChzex(int x,int y, int z){
        return this.Chzex[x][y][z];
    }
    
    /////
    
    public void setCexe(int x,int y, int z, double value){
        this.Cexe[x][y][z] = value;
    }
    
    public void setCexhz(int x,int y, int z, double value){
        this.Cexhz[x][y][z] = value;
    }
    
    public void setCexhy(int x,int y, int z, double value){
        this.Cexhy[x][y][z] = value;
    }
    
    public void setCeye(int x,int y, int z, double value){
        this.Ceye[x][y][z] = value;
    }
    
    public void setCeyhx(int x,int y, int z, double value){
        this.Ceyhx[x][y][z] = value;
    }
    
    public void setCeyhz(int x,int y, int z, double value){
        this.Ceyhz[x][y][z] = value;
    }
    
    public void setCeze(int x,int y, int z, double value){
        this.Ceze[x][y][z] = value;
    }
    
    public void setCezhy(int x,int y, int z, double value){
        this.Cezhy[x][y][z] = value;
    }
    
    public void setCezhx(int x,int y, int z, double value){
        this.Cezhx[x][y][z] = value;
    }
    
    public void setChxh(int x,int y, int z, double value){
        this.Chxh[x][y][z] = value;
    }
    
    public void setChxez(int x,int y, int z, double value){
        this.Chxez[x][y][z] = value;
    }
    
    public void setChxey(int x,int y, int z, double value){
        this.Chxey[x][y][z] = value;
    }
    
    public void setChyh(int x,int y, int z, double value){
        this.Chyh[x][y][z] = value;
    }
    
    public void setChyex(int x,int y, int z, double value){
        this.Chyex[x][y][z] = value;
    }
    
    public void setChyez(int x,int y, int z, double value){
        this.Chyez[x][y][z] = value;
    }
    
    public void setChzh(int x,int y, int z, double value){
        this.Chzh[x][y][z] = value;
    }
    
    public void setChzey(int x,int y, int z, double value){
        this.Chzey[x][y][z] = value;
    }
    
    public void setChzex(int x,int y, int z, double value){
        this.Chzex[x][y][z] = value;
    }
}
