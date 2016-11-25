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
    protected static int nr, na, nz;
    protected static double[][][] Hr,Ha,Hz;
    protected static double[][][] Er,Ea,Ez;
    protected static double[][][] Cere,Cerhz,Cerha,Ceae,Ceahr,Ceahz,Ceze,Cezha,Cezhr;
    protected static double[][][] Chrh,Chrez,Chrea,Chah,Chaer,Chaez,Chzh,Chzea,Chzer;
    protected static Cell c;
    
        
    public static void saveAllCoefArrays(){
        Common.save3DArray(Cere,"Cere");
        Common.save3DArray(Cerhz,"Cerhz");
        Common.save3DArray(Cerha,"Cerha");
        Common.save3DArray(Ceae,"Ceae");
        Common.save3DArray(Ceahr,"Ceahr");
        Common.save3DArray(Ceahz,"Ceahz");
        Common.save3DArray(Ceze,"Ceze");
        Common.save3DArray(Cezha,"Cezha");
        Common.save3DArray(Cezhr,"Cezhr");
        
        Common.save3DArray(Chrh,"Chrh");
        Common.save3DArray(Chrez,"Chrez");
        Common.save3DArray(Chrea,"Chrea");
        Common.save3DArray(Chah,"Chah");
        Common.save3DArray(Chaer,"Chaer");
        Common.save3DArray(Chaez,"Chaez");
        Common.save3DArray(Chzh,"Chzh");
        Common.save3DArray(Chzea,"Chzea");
        Common.save3DArray(Chzer,"Chzer");
    }
    
    public EMobject(int nr,int na,int nz,Cell c){
        this.nr = nr;
        this.na = na;
        this.nz = nz;
        this.c = c;
    }
    
    public EMobject(){
        
    }
   
    public double getCere(int r,int a, int z){
        return this.Cere[r][a][z];
    }
    
    public double getCerhz(int r,int a, int z){
        return this.Cerhz[r][a][z];
    }
    
    public double getCerha(int r,int a, int z){
        return this.Cerha[r][a][z];
    }
    
    public double getCeae(int r,int a, int z){
        return this.Ceae[r][a][z];
    }
    
    public double getCeahr(int r,int a, int z){
        return this.Ceahr[r][a][z];
    }
    
    public double getCeahz(int r,int a, int z){
        return this.Ceahz[r][a][z];
    }
    
    public double getCeze(int r,int a, int z){
        return this.Ceze[r][a][z];
    }
    
    public double getCezha(int r,int a, int z){
        return this.Cezha[r][a][z];
    }
    
    public double getCezhr(int r,int a, int z){
        return this.Cezhr[r][a][z];
    }
    
    public double getChrh(int r,int a, int z){
        return this.Chrh[r][a][z];
    }
    
    public double getChrez(int r,int a, int z){
        return this.Chrez[r][a][z];
    }
    
    public double getChrea(int r,int a, int z){
        return this.Chrea[r][a][z];
    }
    
    public double getChah(int r,int a, int z){
        return this.Chah[r][a][z];
    }
    
    public double getChaer(int r,int a, int z){
        return this.Chaer[r][a][z];
    }
    
    public double getChaez(int r,int a, int z){
        return this.Chaez[r][a][z];
    }
    
    public double getChzh(int r,int a, int z){
        return this.Chzh[r][a][z];
    }
    
    public double getChzea(int r,int a, int z){
        return this.Chzea[r][a][z];
    }
    
    public double getChzer(int r,int a, int z){
        return this.Chzer[r][a][z];
    }
    
    /////
    
    public void setCere(int r,int a, int z, double value){
        this.Cere[r][a][z] = value;
    }
    
    public void setCerhz(int r,int a, int z, double value){
        this.Cerhz[r][a][z] = value;
    }
    
    public void setCerha(int r,int a, int z, double value){
        this.Cerha[r][a][z] = value;
    }
    
    public void setCeae(int r,int a, int z, double value){
        this.Ceae[r][a][z] = value;
    }
    
    public void setCeahr(int r,int a, int z, double value){
        this.Ceahr[r][a][z] = value;
    }
    
    public void setCeahz(int r,int a, int z, double value){
        this.Ceahz[r][a][z] = value;
    }
    
    public void setCeze(int r,int a, int z, double value){
        this.Ceze[r][a][z] = value;
    }
    
    public void setCezha(int r,int a, int z, double value){
        this.Cezha[r][a][z] = value;
    }
    
    public void setCezhr(int r,int a, int z, double value){
        this.Cezhr[r][a][z] = value;
    }
    
    public void setChrh(int r,int a, int z, double value){
        this.Chrh[r][a][z] = value;
    }
    
    public void setChrez(int r,int a, int z, double value){
        this.Chrez[r][a][z] = value;
    }
    
    public void setChrea(int r,int a, int z, double value){
        this.Chrea[r][a][z] = value;
    }
    
    public void setChah(int r,int a, int z, double value){
        this.Chah[r][a][z] = value;
    }
    
    public void setChaer(int r,int a, int z, double value){
        this.Chaer[r][a][z] = value;
    }
    
    public void setChaez(int r,int a, int z, double value){
        this.Chaez[r][a][z] = value;
    }
    
    public void setChzh(int r,int a, int z, double value){
        this.Chzh[r][a][z] = value;
    }
    
    public void setChzea(int r,int a, int z, double value){
        this.Chzea[r][a][z] = value;
    }
    
    public void setChzer(int r,int a, int z, double value){
        this.Chzer[r][a][z] = value;
    }
}
