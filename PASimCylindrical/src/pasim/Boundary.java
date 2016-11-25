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
    private boolean cpmlRn;
    private int airBufferRn;
    private int cpmlCellsRn;

    private boolean cpmlRp;
    private int airBufferRp;
    private int cpmlCellsRp;

    private boolean cpmlAn;
    private int airBufferAn;
    private int cpmlCellsAn;

    private boolean cpmlAp;
    private int airBufferAp;
    private int cpmlCellsAp;

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
    
    public void setAirBuffer(int airBufferRn, int airBufferRp, int airBufferAn, int airBufferAp, int airBufferZn, int airBufferZp){
        this.airBufferRn = airBufferRn;
        this.airBufferRp = airBufferRp;
        this.airBufferAn = airBufferAn;
        this.airBufferAp = airBufferAp;
        this.airBufferZn = airBufferZn;
        this.airBufferZp = airBufferZp;
    }
    
    public void setCpmlCellNumber(int cpmlCellsRn, int cpmlCellsRp, int cpmlCellsAn, int cpmlCellsAp, int cpmlCellsZn, int cpmlCellsZp){
        this.cpmlCellsRn = cpmlCellsRn;
        this.cpmlCellsRp = cpmlCellsRp;
        this.cpmlCellsAn = cpmlCellsAn;
        this.cpmlCellsAp = cpmlCellsAp;
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
    
    public void setBoundaryType(boolean cpmlRn, boolean cpmlRp, boolean cpmlAn, boolean cpmlAp, boolean cpmlZn, boolean cpmlZp){
        this.cpmlRn = cpmlRn;
        this.cpmlRp = cpmlRp;
        this.cpmlAn = cpmlAn;
        this.cpmlAp = cpmlAp;
        this.cpmlZn = cpmlZn;
        this.cpmlZp = cpmlZp;
    }
    
    public int getAirBufferCellCountRn(){
        return this.airBufferRn;
    }

    public int getAirBufferCellCountRp(){
        return this.airBufferRp;
    }
    
    public int getAirBufferCellCountAn(){
        return this.airBufferAn;
    }
    
    public int getAirBufferCellCountAp(){
        return this.airBufferAp;
    }
    
    public int getAirBufferCellCountZn(){
        return this.airBufferZn;
    }
    
    public int getAirBufferCellCountZp(){
        return this.airBufferZp;
    }
    
    public boolean getCPMLrn(){
        return this.cpmlRn;
    }
    
    public boolean getCPMLrp(){
        return this.cpmlRp;
    }

    public boolean getCPMLan(){
        return this.cpmlAn;
    }
    
    public boolean getCPMLap(){
        return this.cpmlAp;
    }
    
    public boolean getCPMLzn(){
        return this.cpmlZn;
    }
    
    public boolean getCPMLzp(){
        return this.cpmlZn;
    }
    
    public boolean isAnySideCPML(){
        return cpmlRn|cpmlRp|cpmlAn|cpmlAp|cpmlZn|cpmlZp;
    }
    
    public int getCellCountRn(){
        return this.cpmlCellsRn;
    }

    public int getCellCountRp(){
        return this.cpmlCellsRp;
    }
    
    public int getCellCountAn(){
        return this.cpmlCellsAn;
    }
    
    public int getCellCountAp(){
        return this.cpmlCellsAp;
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