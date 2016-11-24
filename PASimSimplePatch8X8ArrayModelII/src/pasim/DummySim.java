package pasim;

import pasim.*;

public class DummySim{
    private DummyFarField ff;
    private Boundary b;
    private Cell c;
    private ProblemSpace ps;
    private double speed;
    private double courantFactor;
    private int brickCount;
    private Brick[] bricks;
    private Material[] m;
    private double operatingFrequency;
    private int arraySizeX,arraySizeY;
    private InfiniteByInfinite IByI;
    private InfiniteByFinite IByF;
    private FiniteByInfinite FByI;
    private FiniteByFinite FByF;
    
    public DummySim(){
        initializeProblemSpace();
        defineGeometry();
        defineOutput();
        initMaterialGrid();
        initFDTDparams();
        initParam();
        ps=null;
        b=null;
    }
    
    private void simCompArray(){
        IByI = new InfiniteByInfinite(c);
        IByF = new InfiniteByFinite(c);
        
        IByI.getFF().printSize("IByI");
        IByF.getFF().printSize("IByF");
        
        ff.printSize("Dummy");
        
        IByI.start();
        IByF.start();
        
        try {
            IByI.t.join();
            IByF.t.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        
        FByI = new FiniteByInfinite(c);
        FByI.getFF().printSize("FByI");
        FByI.start();
        
        try {
            FByI.t.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        
        FByF = new FiniteByFinite(c);
        FByF.getFF().printSize("FByF");
        FByF.start();
        
        try {
            FByF.t.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        
        /*
        try {
            IByI.t.join();
            IByF.t.join();
            FByI.t.join();
            FByF.t.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        */
        
        this.assignFFdata();
    }
    
    private void assignFFdata(){
        ff.setCjxyp(this.cjxypBuild());
        ff.setCjxzp(this.cjxzpBuild());
        ff.setCjyxp(this.cjyxpBuild());
        ff.setCjyzp(this.cjyzpBuild());
        ff.setCjzxp(this.cjzxpBuild());
        ff.setCjzyp(this.cjzypBuild());
        ff.setCjxyn(this.cjxynBuild());
        ff.setCjxzn(this.cjxznBuild());
        ff.setCjyxn(this.cjyxnBuild());
        ff.setCjyzn(this.cjyznBuild());
        ff.setCjzxn(this.cjzxnBuild());
        ff.setCjzyn(this.cjzynBuild());
        ff.setCmxyp(this.cmxypBuild());
        ff.setCmxzp(this.cmxzpBuild());
        ff.setCmyxp(this.cmyxpBuild());
        ff.setCmyzp(this.cmyzpBuild());
        ff.setCmzxp(this.cmzxpBuild());
        ff.setCmzyp(this.cmzypBuild());
        ff.setCmxyn(this.cmxynBuild());
        ff.setCmxzn(this.cmxznBuild());
        ff.setCmyxn(this.cmyxnBuild());
        ff.setCmyzn(this.cmyznBuild());
        ff.setCmzxn(this.cmzxnBuild());
        ff.setCmzyn(this.cmzynBuild());
    }
    
    private void initializeProblemSpace(){
        courantFactor =0.9;
        operatingFrequency = 2.8e9; // Hz
        arraySizeX = 10;
        arraySizeY = 10;
        c = new Cell(1e-3,1e-3,0.5e-3);
        createAllBoundaries();
        createAllMaterials();
    }
    
    private double elementDimension(double ratio){
        return ratio*wavelengthInFreeSpace();
    }
    
    private double wavelengthInFreeSpace(){
        return Constants.C/operatingFrequency;
    }
    
    private void createAllBoundaries(){
        b = new Boundary();
        b.setAirBuffer(10,10,10);
        b.setCPML(true, true, true);
        b.setPBC(false, false, false);
        b.setCpmlCellNumber(8, 8, 8);
        b.setCpmlParam(3, 1.3, 7.0, 0.0, 0.05);
    }
    
    private void createAllMaterials(){
        Material.operatingFrequency = this.operatingFrequency;
        m = new Material[1];
        m[0] = Material.Air(0);
    }
    
    private void defineGeometry(){
        brickCount =1;
        bricks = new Brick[brickCount];
        bricks[0] = CreateSimulationBox();
    }
    
    private Brick CreateSimulationBox(){
        double dimX = c.getDeltaX()*((int)Math.round(((elementDimension(0.5)+c.getDeltaX())/2)/c.getDeltaX()));
        double dimY = c.getDeltaY()*((int)Math.round(((elementDimension(0.5)+c.getDeltaY())/2)/c.getDeltaY()));
        double dimZ = c.getDeltaZ()*((int)Math.round(((elementDimension(0.5)+c.getDeltaZ())/2)/c.getDeltaZ()));

        Brick b = new Brick();
        b.Material(this.m[0]);
        b.Ymin(-dimY - (27.0e-3)*arraySizeY);
        b.Xmin(-dimX - (27.0e-3)*arraySizeX);
        b.Zmin(-dimZ -1.5e-3);
        b.Ymax(dimY + (27.0e-3)*arraySizeY);
        b.Xmax(dimX + (27.0e-3)*arraySizeX);
        b.Zmax(dimZ + 0);
        return b;
    }
    
    private void defineOutput(){
        ff = new DummyFarField(2.7e9,2.9e9,3,c);
        ff.setCellCountOuterBoundary(13);                
    }
    
    private void initMaterialGrid(){   
        ps = new ProblemSpace(this.brickCount,this.bricks,this.c,this.b);
        ps.initMaterialGrid(true);
    }
    
    private void initFDTDparams(){
        speed = 1/Math.sqrt(Constants.MU0*Constants.EPS0);
        c.setDeltaT(courantFactor,speed);      
    }
    
    private void initParam(){
        ff.initFarFieldArrays(ps,b);
    }
    
    public void run(){
        simCompArray();
        this.postProcessing();
        this.saveData();
    }
    
    ////////////////////////////////////////////////////////////////////////////      
    public Complex[][][][] cjxypBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Iperiodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*Iperiodic + (this.FByF.getFF().ui()-this.FByF.getFF().li()));
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCjxypYmaxEnd();
        Complex[][][][] JYmin = this.getCjxypYminEnd();
        Complex[][][][] JY = this.getCjxypPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    public Complex[][][][] cjxynBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Iperiodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*Iperiodic + (this.FByF.getFF().ui()-this.FByF.getFF().li()));
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCjxynYmaxEnd();
        Complex[][][][] JYmin = this.getCjxynYminEnd();
        Complex[][][][] JY = this.getCjxynPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    public Complex[][][][] cjzypBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Iperiodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*Iperiodic + (this.FByF.getFF().ui()-this.FByF.getFF().li()));
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCjzypYmaxEnd();
        Complex[][][][] JYmin = this.getCjzypYminEnd();
        Complex[][][][] JY = this.getCjzypPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    public Complex[][][][] cjzynBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Iperiodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*Iperiodic + (this.FByF.getFF().ui()-this.FByF.getFF().li()));
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCjzynYmaxEnd();
        Complex[][][][] JYmin = this.getCjzynYminEnd();
        Complex[][][][] JY = this.getCjzynPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    
    public Complex[][][][] cjyxpBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Jperiodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*Jperiodic + (this.FByF.getFF().uj()-this.FByF.getFF().lj()));
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int I = 1;
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCjyxpYmaxEnd();
        Complex[][][][] JYmin = this.getCjyxpYminEnd();
        Complex[][][][] JY = this.getCjyxpPeriodic();
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    public Complex[][][][] cjyxnBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Jperiodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*Jperiodic + (this.FByF.getFF().uj()-this.FByF.getFF().lj()));
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int I = 1;
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCjyxnYmaxEnd();
        Complex[][][][] JYmin = this.getCjyxnYminEnd();
        Complex[][][][] JY = this.getCjyxnPeriodic();
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    public Complex[][][][] cjzxpBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Jperiodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*Jperiodic + (this.FByF.getFF().uj()-this.FByF.getFF().lj()));
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int I = 1;
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCjzxpYmaxEnd();
        Complex[][][][] JYmin = this.getCjzxpYminEnd();
        Complex[][][][] JY = this.getCjzxpPeriodic();
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    
    public Complex[][][][] cjzxnBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Jperiodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*Jperiodic + (this.FByF.getFF().uj()-this.FByF.getFF().lj()));
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int I = 1;
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCjzxnYmaxEnd();
        Complex[][][][] JYmin = this.getCjzxnYminEnd();
        Complex[][][][] JY = this.getCjzxnPeriodic();
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    
    public Complex[][][][] cjxzpBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        
        int IedgePeriodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int Iperiodic = this.IByI.getFF().ui() - this.IByI.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*IedgePeriodic + 2*Iends);
        
        int JedgePeriodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int Jperiodic = this.IByI.getFF().uj() - this.IByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*JedgePeriodic + 2*Jends);
        
        int K = 1;
        
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JXminYmin = this.getCjxzpXminYminCorner();
        Complex[][][][] JXmaxYmin = this.getCjxzpXmaxYminCorner();
        Complex[][][][] JXminYmax = this.getCjxzpXminYmaxCorner();
        Complex[][][][] JXmaxYmax = this.getCjxzpXmaxYmaxCorner();
        Complex[][][][] JYminEdge = this.getCjxzpYminEdge();
        Complex[][][][] JYmaxEdge = this.getCjxzpYmaxEdge();
        Complex[][][][] JXminEdge = this.getCjxzpXminEdge();
        Complex[][][][] JXmaxEdge = this.getCjxzpXmaxEdge();
        Complex[][][][] JInfByInf = this.getCjxzpPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JXminYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JXminYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JXmaxYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][J-Jends+j][k][f] = JXmaxYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JYminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][J-Jends+j][k][f] = JYmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JXminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[I-Iends+i][Jends + yCount*Jperiodic + j][k][f] = JXmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeY-arrayFbyFSize) ; yCount++){
                            for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                                val[Iends + xCount*Iperiodic + i][Jends + yCount*Jperiodic + j][k][f] = JInfByInf[i][j][k][f];
                            }
                        }
                    }
                }
            }
        }
         
        return val;
    }
    
    public Complex[][][][] cjxznBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        
        int IedgePeriodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int Iperiodic = this.IByI.getFF().ui() - this.IByI.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*IedgePeriodic + 2*Iends);
        
        int JedgePeriodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int Jperiodic = this.IByI.getFF().uj() - this.IByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*JedgePeriodic + 2*Jends);
        
        int K = 1;
        
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JXminYmin = this.getCjxznXminYminCorner();
        Complex[][][][] JXmaxYmin = this.getCjxznXmaxYminCorner();
        Complex[][][][] JXminYmax = this.getCjxznXminYmaxCorner();
        Complex[][][][] JXmaxYmax = this.getCjxznXmaxYmaxCorner();
        Complex[][][][] JYminEdge = this.getCjxznYminEdge();
        Complex[][][][] JYmaxEdge = this.getCjxznYmaxEdge();
        Complex[][][][] JXminEdge = this.getCjxznXminEdge();
        Complex[][][][] JXmaxEdge = this.getCjxznXmaxEdge();
        Complex[][][][] JInfByInf = this.getCjxznPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JXminYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JXminYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JXmaxYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][J-Jends+j][k][f] = JXmaxYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JYminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][J-Jends+j][k][f] = JYmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JXminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[I-Iends+i][Jends + yCount*Jperiodic + j][k][f] = JXmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeY-arrayFbyFSize) ; yCount++){
                            for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                                val[Iends + xCount*Iperiodic + i][Jends + yCount*Jperiodic + j][k][f] = JInfByInf[i][j][k][f];
                            }
                        }
                    }
                }
            }
        }
         
        return val;
    }
    
    public Complex[][][][] cjyzpBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        
        int IedgePeriodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int Iperiodic = this.IByI.getFF().ui() - this.IByI.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*IedgePeriodic + 2*Iends);
        
        int JedgePeriodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int Jperiodic = this.IByI.getFF().uj() - this.IByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*JedgePeriodic + 2*Jends);
        
        int K = 1;
        
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JXminYmin = this.getCjyzpXminYminCorner();
        Complex[][][][] JXmaxYmin = this.getCjyzpXmaxYminCorner();
        Complex[][][][] JXminYmax = this.getCjyzpXminYmaxCorner();
        Complex[][][][] JXmaxYmax = this.getCjyzpXmaxYmaxCorner();
        Complex[][][][] JYminEdge = this.getCjyzpYminEdge();
        Complex[][][][] JYmaxEdge = this.getCjyzpYmaxEdge();
        Complex[][][][] JXminEdge = this.getCjyzpXminEdge();
        Complex[][][][] JXmaxEdge = this.getCjyzpXmaxEdge();
        Complex[][][][] JInfByInf = this.getCjyzpPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JXminYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JXminYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JXmaxYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][J-Jends+j][k][f] = JXmaxYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JYminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][J-Jends+j][k][f] = JYmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JXminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[I-Iends+i][Jends + yCount*Jperiodic + j][k][f] = JXmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeY-arrayFbyFSize) ; yCount++){
                            for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                                val[Iends + xCount*Iperiodic + i][Jends + yCount*Jperiodic + j][k][f] = JInfByInf[i][j][k][f];
                            }
                        }
                    }
                }
            }
        }
         
        return val;
    }
    
    public Complex[][][][] cjyznBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        
        int IedgePeriodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int Iperiodic = this.IByI.getFF().ui() - this.IByI.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*IedgePeriodic + 2*Iends);
        
        int JedgePeriodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int Jperiodic = this.IByI.getFF().uj() - this.IByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*JedgePeriodic + 2*Jends);
        
        int K = 1;
        
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JXminYmin = this.getCjyznXminYminCorner();
        Complex[][][][] JXmaxYmin = this.getCjyznXmaxYminCorner();
        Complex[][][][] JXminYmax = this.getCjyznXminYmaxCorner();
        Complex[][][][] JXmaxYmax = this.getCjyznXmaxYmaxCorner();
        Complex[][][][] JYminEdge = this.getCjyznYminEdge();
        Complex[][][][] JYmaxEdge = this.getCjyznYmaxEdge();
        Complex[][][][] JXminEdge = this.getCjyznXminEdge();
        Complex[][][][] JXmaxEdge = this.getCjyznXmaxEdge();
        Complex[][][][] JInfByInf = this.getCjyznPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JXminYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JXminYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JXmaxYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][J-Jends+j][k][f] = JXmaxYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JYminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][J-Jends+j][k][f] = JYmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JXminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[I-Iends+i][Jends + yCount*Jperiodic + j][k][f] = JXmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeY-arrayFbyFSize) ; yCount++){
                            for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                                val[Iends + xCount*Iperiodic + i][Jends + yCount*Jperiodic + j][k][f] = JInfByInf[i][j][k][f];
                            }
                        }
                    }
                }
            }
        }
         
        return val;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////      
    public Complex[][][][] cmxypBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Iperiodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*Iperiodic + (this.FByF.getFF().ui()-this.FByF.getFF().li()));
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCmxypYmaxEnd();
        Complex[][][][] JYmin = this.getCmxypYminEnd();
        Complex[][][][] JY = this.getCmxypPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    public Complex[][][][] cmxynBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Iperiodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*Iperiodic + (this.FByF.getFF().ui()-this.FByF.getFF().li()));
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCmxynYmaxEnd();
        Complex[][][][] JYmin = this.getCmxynYminEnd();
        Complex[][][][] JY = this.getCmxynPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    
    public Complex[][][][] cmzypBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Iperiodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*Iperiodic + (this.FByF.getFF().ui()-this.FByF.getFF().li()));
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCmzypYmaxEnd();
        Complex[][][][] JYmin = this.getCmzypYminEnd();
        Complex[][][][] JY = this.getCmzypPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    public Complex[][][][] cmzynBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Iperiodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*Iperiodic + (this.FByF.getFF().ui()-this.FByF.getFF().li()));
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCmzynYmaxEnd();
        Complex[][][][] JYmin = this.getCmzynYminEnd();
        Complex[][][][] JY = this.getCmzynPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    
    public Complex[][][][] cmyxpBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Jperiodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*Jperiodic + (this.FByF.getFF().uj()-this.FByF.getFF().lj()));
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int I = 1;
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCmyxpYmaxEnd();
        Complex[][][][] JYmin = this.getCmyxpYminEnd();
        Complex[][][][] JY = this.getCmyxpPeriodic();
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    public Complex[][][][] cmyxnBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Jperiodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*Jperiodic + (this.FByF.getFF().uj()-this.FByF.getFF().lj()));
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int I = 1;
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCmyxnYmaxEnd();
        Complex[][][][] JYmin = this.getCmyxnYminEnd();
        Complex[][][][] JY = this.getCmyxnPeriodic();
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    public Complex[][][][] cmzxpBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Jperiodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*Jperiodic + (this.FByF.getFF().uj()-this.FByF.getFF().lj()));
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int I = 1;
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCmzxpYmaxEnd();
        Complex[][][][] JYmin = this.getCmzxpYminEnd();
        Complex[][][][] JY = this.getCmzxpPeriodic();
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    
    public Complex[][][][] cmzxnBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        int Jperiodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*Jperiodic + (this.FByF.getFF().uj()-this.FByF.getFF().lj()));
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int I = 1;
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JYmax = this.getCmzxnYmaxEnd();
        Complex[][][][] JYmin = this.getCmzxnYminEnd();
        Complex[][][][] JY = this.getCmzxnPeriodic();
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JY[i][j][k][f];
                        }
                    }
                }
            }
        }   
         
        return val;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    
    public Complex[][][][] cmxzpBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        
        int IedgePeriodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int Iperiodic = this.IByI.getFF().ui() - this.IByI.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*IedgePeriodic + 2*Iends);
        
        int JedgePeriodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int Jperiodic = this.IByI.getFF().uj() - this.IByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*JedgePeriodic + 2*Jends);
        
        int K = 1;
        
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JXminYmin = this.getCmxzpXminYminCorner();
        Complex[][][][] JXmaxYmin = this.getCmxzpXmaxYminCorner();
        Complex[][][][] JXminYmax = this.getCmxzpXminYmaxCorner();
        Complex[][][][] JXmaxYmax = this.getCmxzpXmaxYmaxCorner();
        Complex[][][][] JYminEdge = this.getCmxzpYminEdge();
        Complex[][][][] JYmaxEdge = this.getCmxzpYmaxEdge();
        Complex[][][][] JXminEdge = this.getCmxzpXminEdge();
        Complex[][][][] JXmaxEdge = this.getCmxzpXmaxEdge();
        Complex[][][][] JInfByInf = this.getCmxzpPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JXminYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JXminYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JXmaxYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][J-Jends+j][k][f] = JXmaxYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JYminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][J-Jends+j][k][f] = JYmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JXminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[I-Iends+i][Jends + yCount*Jperiodic + j][k][f] = JXmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeY-arrayFbyFSize) ; yCount++){
                            for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                                val[Iends + xCount*Iperiodic + i][Jends + yCount*Jperiodic + j][k][f] = JInfByInf[i][j][k][f];
                            }
                        }
                    }
                }
            }
        }
         
        return val;
    }
    
    public Complex[][][][] cmxznBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        
        int IedgePeriodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int Iperiodic = this.IByI.getFF().ui() - this.IByI.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*IedgePeriodic + 2*Iends);
        
        int JedgePeriodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int Jperiodic = this.IByI.getFF().uj() - this.IByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*JedgePeriodic + 2*Jends);
        
        int K = 1;
        
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JXminYmin = this.getCmxznXminYminCorner();
        Complex[][][][] JXmaxYmin = this.getCmxznXmaxYminCorner();
        Complex[][][][] JXminYmax = this.getCmxznXminYmaxCorner();
        Complex[][][][] JXmaxYmax = this.getCmxznXmaxYmaxCorner();
        Complex[][][][] JYminEdge = this.getCmxznYminEdge();
        Complex[][][][] JYmaxEdge = this.getCmxznYmaxEdge();
        Complex[][][][] JXminEdge = this.getCmxznXminEdge();
        Complex[][][][] JXmaxEdge = this.getCmxznXmaxEdge();
        Complex[][][][] JInfByInf = this.getCmxznPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JXminYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JXminYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JXmaxYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][J-Jends+j][k][f] = JXmaxYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JYminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][J-Jends+j][k][f] = JYmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JXminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[I-Iends+i][Jends + yCount*Jperiodic + j][k][f] = JXmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeY-arrayFbyFSize) ; yCount++){
                            for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                                val[Iends + xCount*Iperiodic + i][Jends + yCount*Jperiodic + j][k][f] = JInfByInf[i][j][k][f];
                            }
                        }
                    }
                }
            }
        }
         
        return val;
    }
    
    public Complex[][][][] cmyzpBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        
        int IedgePeriodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int Iperiodic = this.IByI.getFF().ui() - this.IByI.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*IedgePeriodic + 2*Iends);
        
        int JedgePeriodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int Jperiodic = this.IByI.getFF().uj() - this.IByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*JedgePeriodic + 2*Jends);
        
        int K = 1;
        
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JXminYmin = this.getCmyzpXminYminCorner();
        Complex[][][][] JXmaxYmin = this.getCmyzpXmaxYminCorner();
        Complex[][][][] JXminYmax = this.getCmyzpXminYmaxCorner();
        Complex[][][][] JXmaxYmax = this.getCmyzpXmaxYmaxCorner();
        Complex[][][][] JYminEdge = this.getCmyzpYminEdge();
        Complex[][][][] JYmaxEdge = this.getCmyzpYmaxEdge();
        Complex[][][][] JXminEdge = this.getCmyzpXminEdge();
        Complex[][][][] JXmaxEdge = this.getCmyzpXmaxEdge();
        Complex[][][][] JInfByInf = this.getCmyzpPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JXminYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JXminYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JXmaxYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][J-Jends+j][k][f] = JXmaxYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JYminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][J-Jends+j][k][f] = JYmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JXminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[I-Jends+i][Jends + yCount*Jperiodic + j][k][f] = JXmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeY-arrayFbyFSize) ; yCount++){
                            for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                                val[Iends + xCount*Iperiodic + i][Jends + yCount*Jperiodic + j][k][f] = JInfByInf[i][j][k][f];
                            }
                        }
                    }
                }
            }
        }
         
        return val;
    }
    
    public Complex[][][][] cmyznBuild(){
        int arrayFbyFSize = this.FByF.arraySideCount();
        
        int IedgePeriodic = this.IByF.getFF().ui() - this.IByF.getFF().li();
        int Iends = (this.FByF.getFF().ui()-this.FByF.getFF().li())/2;
        int Iperiodic = this.IByI.getFF().ui() - this.IByI.getFF().li();
        int I = ((this.arraySizeX-arrayFbyFSize)*IedgePeriodic + 2*Iends);
        
        int JedgePeriodic = this.FByI.getFF().uj() - this.FByI.getFF().lj();
        int Jends = (this.FByF.getFF().uj()-this.FByF.getFF().lj())/2;
        int Jperiodic = this.IByI.getFF().uj() - this.IByI.getFF().lj();
        int J = ((this.arraySizeY-arrayFbyFSize)*JedgePeriodic + 2*Jends);
        
        int K = 1;
        
        int F = this.IByF.getFF().ffFreqCount();
        
        Complex[][][][] val = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] JXminYmin = this.getCmyznXminYminCorner();
        Complex[][][][] JXmaxYmin = this.getCmyznXmaxYminCorner();
        Complex[][][][] JXminYmax = this.getCmyznXminYmaxCorner();
        Complex[][][][] JXmaxYmax = this.getCmyznXmaxYmaxCorner();
        Complex[][][][] JYminEdge = this.getCmyznYminEdge();
        Complex[][][][] JYmaxEdge = this.getCmyznYmaxEdge();
        Complex[][][][] JXminEdge = this.getCmyznXminEdge();
        Complex[][][][] JXmaxEdge = this.getCmyznXmaxEdge();
        Complex[][][][] JInfByInf = this.getCmyznPeriodic();
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][j][k][f] = JXminYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[i][J-Jends+j][k][f] = JXminYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][j][k][f] = JXmaxYmin[i][j][k][f];
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        val[I-Iends+i][J-Jends+j][k][f] = JXmaxYmax[i][j][k][f];
                    }
                }
            }
        }        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][j][k][f] = JYminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jends ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                            val[Iends + xCount*Iperiodic + i][J-Jends+j][k][f] = JYmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[i][Jends + yCount*Jperiodic + j][k][f] = JXminEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iends ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeX-arrayFbyFSize) ; yCount++){
                            val[I-Iends+i][Jends + yCount*Jperiodic + j][k][f] = JXmaxEdge[i][j][k][f];
                        }
                    }
                }
            }
        }
        
        for(int i=0 ; i<Iperiodic ; i++){
            for(int j=0 ; j<Jperiodic ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        for(int yCount=0 ; yCount<(this.arraySizeY-arrayFbyFSize) ; yCount++){
                            for(int xCount=0 ; xCount<(this.arraySizeX-arrayFbyFSize) ; xCount++){
                                val[Iends + xCount*Iperiodic + i][Jends + yCount*Jperiodic + j][k][f] = JInfByInf[i][j][k][f];
                            }
                        }
                    }
                }
            }
        }
         
        return val;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    /*
    Y direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCjxypYmaxEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjxyp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cjxyp at y=Ymin edge
    */
    private Complex[][][][] getCjxypYminEnd(){        
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjxyp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cjxyn at y=Ymax edge
    */
    private Complex[][][][] getCjxynYmaxEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;       
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjxyn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cjxyn at y=Ymin edge
    */
    private Complex[][][][] getCjxynYminEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjxyn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cmxyp at y=Ymax edge
    */
    private Complex[][][][] getCmxypYmaxEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmxyp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cmxyp at y=Ymin edge
    */
    private Complex[][][][] getCmxypYminEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmxyp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cmxyn at y=Ymax edge
    */
    private Complex[][][][] getCmxynYmaxEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;       
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmxyn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cmxyn at y=Ymin edge
    */
    private Complex[][][][] getCmxynYminEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;       
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmxyn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cjxyp
    */
    private Complex[][][][] getCjxypPeriodic(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = 1;
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cjxyp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cjxyn
    */
    private Complex[][][][] getCjxynPeriodic(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = 1;       
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cjxyn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
        /*
    Periodic edge element Cmxyp
    */
    private Complex[][][][] getCmxypPeriodic(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = 1;
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cmxyp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cmxyn
    */
    private Complex[][][][] getCmxynPeriodic(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = 1;       
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cmxyn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    /*
    X direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCjzypYmaxEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjzyp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    X direction cjxyp at y=Ymin edge
    */
    private Complex[][][][] getCjzypYminEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjzyp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    X direction cjxyn at y=Ymax edge
    */
    private Complex[][][][] getCjzynYmaxEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;       
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjzyn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    X direction cjxyn at y=Ymin edge
    */
    private Complex[][][][] getCjzynYminEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjzyn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    X direction cmxyp at y=Ymax edge
    */
    private Complex[][][][] getCmzypYmaxEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmzyp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmzypYminEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmzyp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmzynYmaxEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;       
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmzyn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmzynYminEnd(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = 1;       
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmzyn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cjxyp
    */
    private Complex[][][][] getCjzypPeriodic(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = 1;
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cjzyp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cjxyn
    */
    private Complex[][][][] getCjzynPeriodic(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = 1;       
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cjzyn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
        /*
    Periodic edge element Cjxyp
    */
    private Complex[][][][] getCmzypPeriodic(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = 1;
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cmzyp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cjxyn
    */
    private Complex[][][][] getCmzynPeriodic(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = 1;       
        int K = this.IByF.getFF().uk() - this.IByF.getFF().lk();
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cmzyn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    /*
    X direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCjyxpYmaxEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjyxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    X direction cjxyp at y=Ymin edge
    */
    private Complex[][][][] getCjyxpYminEnd(){        
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjyxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    X direction cjxyn at y=Ymax edge
    */
    private Complex[][][][] getCjyxnYmaxEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjyxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    X direction cjxyn at y=Ymin edge
    */
    private Complex[][][][] getCjyxnYminEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjyxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    X direction cmxyp at y=Ymax edge
    */
    private Complex[][][][] getCmyxpYmaxEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmyxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyxpYminEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmyxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyxnYmaxEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmyxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyxnYminEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmyxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cjxyp
    */
    private Complex[][][][] getCjyxpPeriodic(){
        int I = 1;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cmyxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cjxyn
    */
    private Complex[][][][] getCjyxnPeriodic(){
        int I = 1;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjyxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
        /*
    Periodic edge element Cjxyp
    */
    private Complex[][][][] getCmyxpPeriodic(){
        int I = 1;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjyxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cjxyn
    */
    private Complex[][][][] getCmyxnPeriodic(){
        int I = 1;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjyxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    ////////////////////////////////////////////////////////////////////////////
    /*
    X direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCjzxpYmaxEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjzxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    X direction cjxyp at y=Ymin edge
    */
    private Complex[][][][] getCjzxpYminEnd(){        
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjzxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    X direction cjxyn at y=Ymax edge
    */
    private Complex[][][][] getCjzxnYmaxEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjzxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    X direction cjxyn at y=Ymin edge
    */
    private Complex[][][][] getCjzxnYminEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjzxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    X direction cmxyp at y=Ymax edge
    */
    private Complex[][][][] getCmzxpYmaxEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmzxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmzxpYminEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmzxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmzxnYmaxEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmzxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmzxnYminEnd(){
        int I = 1;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = this.FByF.getFF().uk() - this.FByF.getFF().lk();
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmzxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cjxyp
    */
    private Complex[][][][] getCjzxpPeriodic(){
        int I = 1;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cmzxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cjxyn
    */
    private Complex[][][][] getCjzxnPeriodic(){
        int I = 1;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjzxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
        /*
    Periodic edge element Cjxyp
    */
    private Complex[][][][] getCmzxpPeriodic(){
        int I = 1;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjzxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cjxyn
    */
    private Complex[][][][] getCmzxnPeriodic(){
        int I = 1;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = this.FByI.getFF().uk() - this.FByI.getFF().lk();
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjzxp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    /*
    Y direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCjxzpXmaxYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjxzp();
        for(int i=0 ; i<I ; i++){
            for(int j=0 ; j<J ; j++){
                for(int k=0 ; k<K ; k++){
                    for(int f=0 ; f<F ; f++){
                        cut[i][j][k][f] = oragin[I+i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjxzpXminYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjxzpXmaxYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjxzpXminYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjxznXmaxYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjxznXminYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjxznXmaxYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjxznXminYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCmxzpXmaxYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmxzpXminYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmxzpXmaxYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmxzpXminYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
        private Complex[][][][] getCmxznXmaxYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmxznXminYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmxznXmaxYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmxznXminYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cjxzp
    */
    private Complex[][][][] getCjxzpPeriodic(){
        return this.IByI.getFF().cjxzp();
    }
    
    /*
    Periodic edge element Cjxzn
    */
    private Complex[][][][] getCjxznPeriodic(){
        return this.IByI.getFF().cjxzn();
    }
    
    /*
    Periodic edge element Cmxzp
    */
    private Complex[][][][] getCmxzpPeriodic(){
        return this.IByI.getFF().cmxzp();
    }
    
    /*
    Periodic edge element Cmxzn
    */
    private Complex[][][][] getCmxznPeriodic(){
        return this.IByI.getFF().cmxzn();
    }
    
    /////////////////////////////
    
    /*
    Y direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCjxzpXmaxEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int Il = ((this.FByI.getFF().ui() - this.FByI.getFF().li()) - I);
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[Il+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjxzpXminEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjxzpYmaxEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int Jl = ((this.IByF.getFF().uj() - this.IByF.getFF().lj()) - J);
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cjxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][Jl+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjxzpYminEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cjxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCjxznXmaxEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int Il = ((this.FByI.getFF().ui() - this.FByI.getFF().li()) - I);
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[Il+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjxznXminEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjxznYmaxEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int Jl = ((this.IByF.getFF().uj() - this.IByF.getFF().lj()) - J);
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cjxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][Jl+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjxznYminEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cjxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCmxzpXmaxEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int Il = ((this.FByI.getFF().ui() - this.FByI.getFF().li()) - I);
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cmxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[Il+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmxzpXminEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cmxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmxzpYmaxEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int Jl = ((this.IByF.getFF().uj() - this.IByF.getFF().lj()) - J);
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cmxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][Jl+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmxzpYminEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cmxzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCmxznXmaxEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int Il = ((this.FByI.getFF().ui() - this.FByI.getFF().li()) - I);
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cmxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[Il+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmxznXminEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cmxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmxznYmaxEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int Jl = ((this.IByF.getFF().uj() - this.IByF.getFF().lj()) - J);
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cmxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][Jl+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmxznYminEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cmxzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    ////////////////////////////////////////////////////////////////////////////

    /*
    Y direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCjyzpXmaxYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjyzpXminYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjyzpXmaxYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjyzpXminYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjyznXmaxYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjyznXminYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjyznXmaxYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjyznXminYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cjyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCmyzpXmaxYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyzpXminYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyzpXmaxYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyzpXminYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyznXmaxYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyznXminYmaxCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][J+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyznXmaxYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[I+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyznXminYminCorner(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;;
        int K = 1;
        int F = this.FByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByF.getFF().cmyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Periodic edge element Cjxzp
    */
    private Complex[][][][] getCjyzpPeriodic(){
        return this.IByI.getFF().cjyzp();
    }
    
    /*
    Periodic edge element Cjxzn
    */
    private Complex[][][][] getCjyznPeriodic(){
        return this.IByI.getFF().cjyzn();
    }
    
    /*
    Periodic edge element Cmxzp
    */
    private Complex[][][][] getCmyzpPeriodic(){
        return this.IByI.getFF().cmyzp();
    }
    
    /*
    Periodic edge element Cmxzn
    */
    private Complex[][][][] getCmyznPeriodic(){
        return this.IByI.getFF().cmyzn();
    }
    
    //////////////////////
    
   /*
    Y direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCjyzpXmaxEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int Il = ((this.FByI.getFF().ui() - this.FByI.getFF().li()) - I);
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[Il+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjyzpXminEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjyzpYmaxEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int Jl = ((this.IByF.getFF().uj() - this.IByF.getFF().lj()) - J);
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cjyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][Jl+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjyzpYminEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cjyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCjyznXmaxEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int Il = ((this.FByI.getFF().ui() - this.FByI.getFF().li()) - I);
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[Il+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjyznXminEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cjyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjyznYmaxEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int Jl = ((this.IByF.getFF().uj() - this.IByF.getFF().lj()) - J);
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cjyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][Jl+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCjyznYminEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cjyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCmyzpXmaxEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int Il = ((this.FByI.getFF().ui() - this.FByI.getFF().li()) - I);
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cmyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[Il+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyzpXminEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cmyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyzpYmaxEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int Jl = ((this.IByF.getFF().uj() - this.IByF.getFF().lj()) - J);
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cmyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][Jl+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyzpYminEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cmyzp();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    /*
    Y direction cjxyp at y=Ymax edge
    */
    private Complex[][][][] getCmyznXmaxEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int Il = ((this.FByI.getFF().ui() - this.FByI.getFF().li()) - I);
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cmyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[Il+i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyznXminEdge(){
        int I = (this.FByF.getFF().ui() - this.FByF.getFF().li())/2;
        int J = (this.FByI.getFF().uj() - this.FByI.getFF().lj());
        int K = 1;
        int F = this.FByI.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.FByI.getFF().cmyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyznYmaxEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int Jl = ((this.IByF.getFF().uj()-this.IByF.getFF().lj()) - J);
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cmyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][Jl+j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    
    private Complex[][][][] getCmyznYminEdge(){
        int I = (this.IByF.getFF().ui() - this.IByF.getFF().li());
        int J = (this.FByF.getFF().uj() - this.FByF.getFF().lj())/2;
        int K = 1;
        int F = this.IByF.getFF().ffFreqCount();
        Complex[][][][] cut = Common.genComplex4DArray(I , J , K , F , new Complex(0.0,0.0));
        Complex[][][][] oragin = this.IByF.getFF().cmyzn();
        for(int i = 0 ; i<I ; i++){
            for(int j = 0 ; j<J ; j++){
                for(int k = 0 ; k<K ; k++){
                    for(int f=0;f<F;f++){
                        cut[i][j][k][f] = oragin[i][j][k][f];
                    }
                }
            }
        }
        
        return cut;
    }
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    private void postProcessing(){
        ff.calTotalRadiatedPower();
        ff.setupTheta(-180, 179, 360);
        ff.setupPhi(0, 0, 360);
        ff.calFarField();
    }

    public void saveData(){
        ff.saveTotalRadiatedPower();
        ff.saveFarFieldThetaPhiPattern();
        this.radXcut();
        this.radYcut();
        this.radDcut();
        this.rad3D();
    }
    
    public void radXcut(){
        ff.setupTheta(-180, 179, 360);
        ff.setupPhi(0, 0, 360);
        ff.calFarField();
        ff.saveFarFieldThetaPhiPattern("FarFieldPhiPatternX","FarFieldThetaPatternX","FarFIeldPhiX","FarFieldThetaX");
    }
    
    public void radYcut(){
        ff.setupTheta(-180, 179, 360);
        ff.setupPhi(45, 45, 360);
        ff.calFarField();
        ff.saveFarFieldThetaPhiPattern("FarFieldPhiPatternD","FarFieldThetaPatternD","FarFIeldPhiD","FarFieldThetaD");
    }
    
    public void radDcut(){
        ff.setupTheta(-180, 179, 360);
        ff.setupPhi(90, 90, 360);
        ff.calFarField();
        ff.saveFarFieldThetaPhiPattern("FarFieldPhiPatternY","FarFieldThetaPatternY","FarFIeldPhiY","FarFieldThetaY");
    }

    public void rad3D(){
        ff.setupTheta(0, 180, 181);
        ff.setupPhi(0, 360, 361);
        ff.calFarField3D();
        ff.saveFarFieldThetaPhi3DPattern("DPhi3D","DTheta3D","Phi3D","Theta3D");
    }
}