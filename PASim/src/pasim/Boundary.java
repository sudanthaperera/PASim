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
public class Boundary {
    private boolean cpmlXn;
    private boolean pbcX;
    private int airBufferXn;
    private int cpmlCellsXn;

    private boolean cpmlXp;
    private int airBufferXp;
    private int cpmlCellsXp;

    private boolean cpmlYn;
    private boolean pbcY;
    private int airBufferYn;
    private int cpmlCellsYn;

    private boolean cpmlYp;
    private int airBufferYp;
    private int cpmlCellsYp;

    private boolean cpmlZn;
    private boolean pbcZ;
    private int airBufferZn;
    private int cpmlCellsZn;

    private boolean cpmlZp;
    private int airBufferZp;
    private int cpmlCellsZp;

    private int order; 
    private double sigmaFactor;
    private double kappaMax;
    private double alphaMin;
    private double alphaMax;
    
    public Boundary(){
        this.setAirBuffer(10,10,10);
        this.setCPML(true, true, true);
        this.setCpmlCellNumber(8, 8, 8);
        this.setCpmlParam(3, 1.3, 7.0, 0.0, 0.05);
    }
    
    public void setAirBuffer(int airBufferX, int airBufferY, int airBufferZ){
        this.airBufferXn = airBufferX;
        this.airBufferXp = airBufferX;
        this.airBufferYn = airBufferY;
        this.airBufferYp = airBufferY;
        this.airBufferZn = airBufferZ;
        this.airBufferZp = airBufferZ;
    }
    
    public void setAirBufferX(int airBufferX){
        this.airBufferXn = airBufferX;
        this.airBufferXp = airBufferX;
    }
    
    public void setAirBufferY(int airBufferY){
        this.airBufferYn = airBufferY;
        this.airBufferYp = airBufferY;
    }
    
    public void setAirBufferZ(int airBufferZ){
        this.airBufferZn = airBufferZ;
        this.airBufferZp = airBufferZ;
    }
    
    public void setCpmlCellNumber(int cpmlCellsX, int cpmlCellsY, int cpmlCellsZ){
        this.cpmlCellsXn = cpmlCellsX;
        this.cpmlCellsXp = cpmlCellsX;
        this.cpmlCellsYn = cpmlCellsY;
        this.cpmlCellsYp = cpmlCellsY;
        this.cpmlCellsZn = cpmlCellsZ;
        this.cpmlCellsZp = cpmlCellsZ;
    }
    
    public void setCpmlCellNumberX(int cpmlCellsX){
        this.cpmlCellsXn = cpmlCellsX;
        this.cpmlCellsXp = cpmlCellsX;
    }
    
    public void setCpmlCellNumberY(int cpmlCellsY){
        this.cpmlCellsYn = cpmlCellsY;
        this.cpmlCellsYp = cpmlCellsY;
    }
    
    public void setCpmlCellNumberZ(int cpmlCellsZ){
        this.cpmlCellsZn = cpmlCellsZ;
        this.cpmlCellsZp = cpmlCellsZ;
    }
    
    public void setCpmlParam(int order,double sigmaFactor,double kappaMax,double alphaMin,double alphaMax){
        this.order = order;
        this.sigmaFactor = sigmaFactor;
        this.kappaMax = kappaMax;
        this.alphaMin = alphaMin;
        this.alphaMax = alphaMax;
    }
    
    public void setOrder(int value){
        this.order = value;
    }
    
    public void setSigmaFactor(double value){
        this.sigmaFactor = value;
    }
    
    public void setKappaMax(double value){
        this.kappaMax = value;
    }
    
    public void setAlphaMax(double value){
        this.alphaMax = value;
    }
    
    public void setAlphaMin(double value){
        this.alphaMin = value;
    }
    
    public void setPBC(boolean pbcX, boolean pbcY, boolean pbcZ){
        this.pbcX = pbcX;
        this.pbcY = pbcY;
        this.pbcZ = pbcZ;
    }
    
    public void setPBCx(boolean pbcX){
        this.pbcX = pbcX;
    }
    
    public void setPBCy(boolean pbcY){
        this.pbcY = pbcY;
    }
    
    public void setPBCz(boolean pbcZ){
        this.pbcZ = pbcZ;
    }
    
    public boolean getPBCX(){
        return this.pbcX;
    }
    
    public boolean getPBCY(){
        return this.pbcY;
    }
    
    public boolean getPBCZ(){
        return this.pbcZ;
    }
    
    public void setCPML(boolean cpmlX, boolean cpmlY, boolean cpmlZ){
        this.cpmlXn = cpmlX;
        this.cpmlXp = cpmlX;
        this.cpmlYn = cpmlY;
        this.cpmlYp = cpmlY;
        this.cpmlZn = cpmlZ;
        this.cpmlZp = cpmlZ;
    }
    
    public void setCPMLx(boolean value){
        this.cpmlXn = value;
        this.cpmlXp = value;
    }
    
    public void setCPMLy(boolean value){
        this.cpmlYn = value;
        this.cpmlYp = value;
    }
    
    public void setCPMLz(boolean value){
        this.cpmlZn = value;
        this.cpmlZp = value;
    }
    
    public int getAirBufferCellCountXn(){
        return this.airBufferXn;
    }

    public int getAirBufferCellCountXp(){
        return this.airBufferXp;
    }
    
    public int getAirBufferCellCountYn(){
        return this.airBufferYn;
    }
    
    public int getAirBufferCellCountYp(){
        return this.airBufferYp;
    }
    
    public int getAirBufferCellCountZn(){
        return this.airBufferZn;
    }
    
    public int getAirBufferCellCountZp(){
        return this.airBufferZp;
    }
    
    public boolean getCPMLxn(){
        return this.cpmlXn;
    }
    
    public boolean getCPMLxp(){
        return this.cpmlXp;
    }

    public boolean getCPMLyn(){
        return this.cpmlYn;
    }
    
    public boolean getCPMLyp(){
        return this.cpmlYp;
    }
    
    public boolean getCPMLzn(){
        return this.cpmlZn;
    }
    
    public boolean getCPMLzp(){
        return this.cpmlZn;
    }
    
    public boolean isAnySideCPML(){
        return cpmlXn|cpmlXp|cpmlYn|cpmlYp|cpmlZn|cpmlZp;
    }
    
    public int getCellCountXn(){
        return this.cpmlCellsXn;
    }

    public int getCellCountXp(){
        return this.cpmlCellsXp;
    }
    
    public int getCellCountYn(){
        return this.cpmlCellsYn;
    }
    
    public int getCellCountYp(){
        return this.cpmlCellsYp;
    }
    
    public int getCellCountZn(){
        return this.cpmlCellsZn;
    }
    
    public int getCellCountZp(){
        return this.cpmlCellsZp;
    }
    
    public int getOrder(){
        return order;
    }
    
    public double getSigmaFactor(){
        return sigmaFactor;
    }
    
    public double getKappaMax(){
        return kappaMax;
    }
    

    public double getAlphaMin(){
        return alphaMin;
    }
    
    public double getAlphaMax(){
        return alphaMax;
    }
}