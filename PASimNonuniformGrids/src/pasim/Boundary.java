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
    private int airBufferXn;
    private int cpmlCellsXn;

    private boolean cpmlXp;
    private int airBufferXp;
    private int cpmlCellsXp;

    private boolean cpmlYn;
    private int airBufferYn;
    private int cpmlCellsYn;

    private boolean cpmlYp;
    private int airBufferYp;
    private int cpmlCellsYp;

    private boolean cpmlZn;
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
        this.setAirBuffer(10,10,10,10,10,10);
        this.setBoundaryType(true, true, true, true, true, true);
        this.setCpmlCellNumber(8, 8, 8, 8, 8, 8);
        this.setCpmlParam(3, 1.3, 7.0, 0.0, 0.05);
    }
    
    public void setAirBuffer(int airBufferXn, int airBufferXp, int airBufferYn, int airBufferYp, int airBufferZn, int airBufferZp){
        this.airBufferXn = airBufferXn;
        this.airBufferXp = airBufferXp;
        this.airBufferYn = airBufferYn;
        this.airBufferYp = airBufferYp;
        this.airBufferZn = airBufferZn;
        this.airBufferZp = airBufferZp;
    }
    
    public void setCpmlCellNumber(int cpmlCellsXn, int cpmlCellsXp, int cpmlCellsYn, int cpmlCellsYp, int cpmlCellsZn, int cpmlCellsZp){
        this.cpmlCellsXn = cpmlCellsXn;
        this.cpmlCellsXp = cpmlCellsXp;
        this.cpmlCellsYn = cpmlCellsYn;
        this.cpmlCellsYp = cpmlCellsYp;
        this.cpmlCellsZn = cpmlCellsZn;
        this.cpmlCellsZp = cpmlCellsZp;
    }
    
    public void setCpmlParam(int order,double sigmaFactor,double kappaMax,double alphaMin,double alphaMax){
        this.order = order; 
        this.sigmaFactor = sigmaFactor;
        this.kappaMax = kappaMax;
        this.alphaMin = alphaMin;
        this.alphaMax = alphaMax;
    }
    
    public void setBoundaryType(boolean cpmlXn, boolean cpmlXp, boolean cpmlYn, boolean cpmlYp, boolean cpmlZn, boolean cpmlZp){
        this.cpmlXn = cpmlXn;
        this.cpmlXp = cpmlXp;
        this.cpmlYn = cpmlYn;
        this.cpmlYp = cpmlYp;
        this.cpmlZn = cpmlZn;
        this.cpmlZp = cpmlZp;
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