package pasim;


/**
 *
 * @author Sudantha
 */
public class Pasim {
    private int timeSteps;
    private double courantFactor;
    private int brickCount, eFieldProbCount, hFieldProbCount, voltageProbCount, currentProbCount, portCount;
    private int gaussianWaveformCount, sinusoidalWaveformCount,cosineModulatedGaussianWaveformCount,normalizedDerivativeGaussianWaveformCount;
    private double speed;
    private double eta0;
    private double[] time;
    private double operatingFrequency, farfieldFrequencyMin, farfieldFrequencyMax;
    
    private HField H;
    private EField E;
    private ProblemSpace ps;
    private PMLmedia pml;
    private Cell c;
    private Boundary b;
    private Material[] m;
    private Brick[] bricks;
    private VoltageSource[] vs;
    private CurrentSource[] cs;
    private Gaussian[] w;
    private Sinusoidal[] sw;
    private CosineModulatedGaussian[] cgw;
    private NormalizedDerivativeGaussian[] ndgw;
    private Position[] currentProbPos, voltageProbPos, eFieldProbPos, hFieldProbPos;
    private VoltageProb[] voltageProbObject;
    private CurrentProb[] currentProbObject;
    private HFieldProb[] hProbObject;
    private EFieldProb[] eProbObject;
    private Port ports[];
    private FarField ff;
    private FrequencyDomain freqDomain;
    private PeriodicBoundary pbc;
    private boolean isActiveElem, isIsoElem;
    
    public Pasim(){
        initializeProblemSpace();
        defineGeometry();
        defineSources();
        defineWaveform();
        defineOutput();
        initMaterialGrid();
        initFDTDparams();
        initFDTDarrays();
        initSource();
        initWaveForm();
        genCoefficients();
        initParam();
    }
    
    private void initializeProblemSpace(){
        timeSteps = 50000;
        courantFactor =0.9;
        operatingFrequency = 2.8e9; // Hz
        farfieldFrequencyMin = 1.800e9; // Hz
        farfieldFrequencyMax = 3.800e9; // Hz
        VoltageSource.voltageSourceCount = 2;
        CurrentSource.currentSourceCount = 0;
        this.isActiveElem = false;
        this.isIsoElem = true;
        double delta = 0.315e-3;
        c = new Cell(delta,delta,delta);
        createAllBoundaries();
        createAllMaterials();
    }
    
    private void createAllBoundaries(){
        if(this.isActiveElem)
            createAllBoundariesForActiveELem();
        if(this.isIsoElem)
            createAllBoundariesForIsoElem();
    }
    
    
    private void createAllBoundariesForIsoElem(){
        b = new Boundary();
        b.setAirBuffer(10,10,10);
        b.setCPML(true, true, true);
        b.setPBC(false, false, false);
        b.setCpmlCellNumber(8, 8, 8);
        b.setCpmlParam(3, 1.3, 7.0, 0.0, 0.05);
    }
    
    private void createAllBoundariesForActiveELem(){
        b = new Boundary();
        b.setAirBuffer(0,0,10);
        b.setCPML(false, false, true);
        b.setPBC(true, true, false);
        b.setCpmlCellNumber(0, 0, 8);
        b.setCpmlParam(3, 1.3, 7.0, 0.0, 0.05);
    }
    
    private void createAllMaterials(){
        Material.operatingFrequency = this.operatingFrequency;
        m = new Material[7];
        m[0] = Material.Vacuum(0);
        m[1] = Material.Air(1);
        m[2] = Material.PEC(2);
        m[3] = Material.PMC(3);
        m[4] = Material.RogersRTduroid5880(4);
        m[5] = Material.Copper(5);
        m[6] = Material.PerfectDielectric(6);
    }
    
    private void defineGeometry(){
        brickCount =4;
        bricks = new Brick[brickCount];
        
        bricks[0] = CreateSimulationBox();
        
        bricks[1] = new Brick();
        bricks[1].Material(this.m[4]);
        bricks[1].Xmin(-27.09e-3);
        bricks[1].Ymin(-27.09e-3);
        bricks[1].Zmin(-1.575e-3);
        bricks[1].Xmax(27.09e-3);
        bricks[1].Ymax(27.09e-3);
        bricks[1].Zmax(0);
        
        bricks[2] = new Brick();
        bricks[2].Material(this.m[5]);
        bricks[2].Xmin(-27.09e-3);
        bricks[2].Ymin(-27.09e-3);
        bricks[2].Zmin(-3.575e-3);
        bricks[2].Xmax(27.09e-3);
        bricks[2].Ymax(27.09e-3);
        bricks[2].Zmax(-1.575e-3);
        
        bricks[3] = new Brick();
        bricks[3].Material(this.m[5]);
        bricks[3].Xmin(-17.325e-3);
        bricks[3].Ymin(-17.325e-3);
        bricks[3].Zmin(0);
        bricks[3].Xmax(17.325e-3);
        bricks[3].Ymax(17.325e-3);
        bricks[3].Zmax(0);
    }
    
    private Brick CreateSimulationBox(){
        double frequencyMin = 2.65e9;
        double airGap = (3e8/frequencyMin)/(2*Math.PI);
        Brick b = new Brick();
        b.Material(this.m[1]);
        
        if(this.isActiveElem){
            b.Xmin(-27.09e-3);
            b.Ymin(-27.09e-3);
            b.Xmax(27.09e-3);
            b.Ymax(27.09e-3);
        }
        
        if(this.isIsoElem){
            b.Xmin(-27.09e-3 - airGap);
            b.Ymin(-27.09e-3 - airGap);
            b.Xmax(27.09e-3 + airGap);
            b.Ymax(27.09e-3 + airGap);
        }
        b.Zmin(-3.575e-3 - airGap);
        b.Zmax(0 + airGap);
        return b;
    }
    
    private void defineSources(){
        vs = new VoltageSource[2];
        
        vs[0] = new VoltageSource(0, 6.825e-3, -1.575e-3, 0, 6.825e-3, 0, c);
        vs[0].setDirection(Constants.ZP);
        vs[0].setResistance(50.0);
        vs[0].setMagnitude(1.0);
        
        vs[1] = new VoltageSource(6.825e-3, 0, -1.575e-3, 6.825e-3, 0, 0, c);
        vs[1].setDirection(Constants.ZP);
        vs[1].setResistance(50.0);
        vs[1].setMagnitude(0.0);

        /*
	if (VoltageSource.voltageSourceCount > 0){
            vs = new VoltageSource[VoltageSource.voltageSourceCount];
            for (int i = 0; i<VoltageSource.voltageSourceCount;i++){
                vs[i] = new VoltageSource(0, 5.25e-3, -1.575e-3, 0, 5.25e-3, 0, c);
                vs[i].setDirection(Constants.ZP);
                vs[i].setResistance(50.0);
                vs[i].setMagnitude(1.0);
            }
	}
                */

	if (CurrentSource.currentSourceCount > 0){
            cs = new CurrentSource[CurrentSource.currentSourceCount];
	}
    }
    
    private void defineWaveform(){
        gaussianWaveformCount = 0;
        sinusoidalWaveformCount = 0;
        cosineModulatedGaussianWaveformCount = 2;
        normalizedDerivativeGaussianWaveformCount = 0;
        
        if (gaussianWaveformCount != 0){
            w = new Gaussian[gaussianWaveformCount];
        
            w[0] = new Gaussian(20,1.0,0.0,timeSteps);
            w[1] = new Gaussian(20,1.0,0.0,timeSteps);
        }
        
        if (sinusoidalWaveformCount != 0){
            sw = new Sinusoidal[sinusoidalWaveformCount];
            
            sw[0] = new Sinusoidal(20,1.0,0.0,timeSteps);
            sw[1] = new Sinusoidal(20,1.0,0.0,timeSteps);
        }
        
        if (cosineModulatedGaussianWaveformCount != 0){
            cgw = new CosineModulatedGaussian[cosineModulatedGaussianWaveformCount];
            
            cgw[0] = new CosineModulatedGaussian(20,1.0,0.0,timeSteps);
            cgw[1] = new CosineModulatedGaussian(20,1.0,0.0,timeSteps);
        }
        
        if (normalizedDerivativeGaussianWaveformCount != 0){
            ndgw = new NormalizedDerivativeGaussian[normalizedDerivativeGaussianWaveformCount];
            
            ndgw[0] = new NormalizedDerivativeGaussian(20,1.0,0.0,timeSteps);
            ndgw[1] = new NormalizedDerivativeGaussian(20,1.0,0.0,timeSteps);
        }
    }
    
    private void defineOutput(){
        eFieldProbCount  = 1;
        hFieldProbCount  = 0;
        voltageProbCount  = 2;
        currentProbCount  = 2;
        portCount = 2;
        
        currentProbPos = new Position[2];
        voltageProbPos = new Position[2];
        eFieldProbPos = new Position[1];
        hFieldProbPos = new Position[1];
        
        currentProbPos[0] = new Position(0, 6.825e-3, -1e-3);
        voltageProbPos[0] = new Position(0, 6.825e-3, 0);
        currentProbPos[1] = new Position(6.825e-3, 0, -1e-3);
        voltageProbPos[1] = new Position(6.825e-3, 0, 0);
        eFieldProbPos[0] = new Position(20e-3, 20e-3, 0e-3);
        hFieldProbPos[0] = new Position(20e-3, 20e-3, 0e-3);
        
        voltageProbObject = new VoltageProb[voltageProbCount];
        currentProbObject = new CurrentProb[currentProbCount];
        ports = new Port[portCount];
        freqDomain = new FrequencyDomain(2e9,3.5e9,1501);
        ff = new FarField(2.7e9,2.9e9,3);
        ff.setCellCountOuterBoundary(13);
        
        if(voltageProbCount > 0){
            for (int i=0;i<voltageProbCount;i++){
                voltageProbObject[i] = new VoltageProb(voltageProbPos[i],Constants.ZP,-1.575e-3);
            }
        }
        
        if(currentProbCount > 0){
            for (int i=0;i<currentProbCount;i++){
                currentProbObject[i] = new CurrentProb(currentProbPos[i],Constants.ZP);
            }
        }
        
        if(portCount > 0){
            for (int i=0;i<portCount;i++){
                ports[i] = new Port(50,true);
            }
        }                
    }
    
    private void initMaterialGrid(){   
        ps = new ProblemSpace(this.brickCount,this.bricks,this.c,this.b,this.m);
	System.out.printf("Number of cells in X direction = %d\n" , ps.getNX());
	System.out.printf("Number of cells in Y direction = %d\n" , ps.getNY());
	System.out.printf("Number of cells in Z direction = %d\n" , ps.getNZ());
	System.out.printf("Number of cells = %d \n" , ps.getNX()*ps.getNY()*ps.getNZ());
        System.out.println("Initializing FDTD grid...");
        ps.initMaterialGrid();
    }
    
    private void initFDTDparams(){
        speed = 1/Math.sqrt(Constants.MU0*Constants.EPS0);
        eta0 = Math.sqrt(Constants.MU0/Constants.EPS0);
        c.setDeltaT(courantFactor,speed);
        time = new double[this.timeSteps];
        
	for(int index = 0; index < this.timeSteps; index++){
            time[index] = (0.5 + index)*c.getDeltaT();
	}        
    }

    private void initFDTDarrays(){
        H = new HField(this.ps.getNX(),this.ps.getNY(),this.ps.getNZ(),c);
        E = new EField(this.ps.getNX(),this.ps.getNY(),this.ps.getNZ(),c);
        
        H.setAllX(0.0);
        H.setAllY(0.0);
        H.setAllZ(0.0);
        
        E.setAllX(0.0);
        E.setAllY(0.0);
        E.setAllZ(0.0);
    }
    
    private void initSource(){
        double insValue;
	for (int ind = 0; ind < this.gaussianWaveformCount ; ind++){
            if (w[ind].getCellsPerWavelength() == 0){
		w[ind].setCellsPerWavelength(Constants.CELLS_PER_WAVELENGTH);
            }

            w[ind].setMaxFreq(speed/(w[ind].getCellsPerWavelength()*c.getMax()));
            w[ind].setTau((w[ind].getCellsPerWavelength()*c.getMax()/(2*speed)));
            w[ind].setT0(4.5*w[ind].getTau());

            for (int time_ind = 0; time_ind < this.timeSteps; time_ind++){
                insValue = Math.exp(-((time[time_ind] - w[ind].getT0())/w[ind].getTau())* ((time[time_ind] - w[ind].getT0())/w[ind].getTau()));
		w[ind].set(insValue, time_ind);
            }
	}
        
        for (int ind = 0; ind < this.sinusoidalWaveformCount ; ind++){
            if (sw[ind].getCellsPerWavelength() == 0){
		sw[ind].setCellsPerWavelength(Constants.CELLS_PER_WAVELENGTH);
            }
            sw[ind].setFreq(2.7e9);
            for (int time_ind = 0; time_ind < this.timeSteps; time_ind++){
                if(time_ind>20&&time_ind<1000){
                    insValue = Math.sin(2*(Math.PI)*sw[ind].getFreq()*(time_ind+1)*c.getDeltaT());
                }
                else{
                    insValue = 0;
                }
		sw[ind].set(insValue, time_ind);
            }
        }
        
        for (int ind = 0; ind < this.cosineModulatedGaussianWaveformCount ; ind++){
            if (cgw[ind].getCellsPerWavelength() == 0){
		cgw[ind].setCellsPerWavelength(Constants.CELLS_PER_WAVELENGTH);
            }


            cgw[ind].setMaxFreq(speed/(cgw[ind].getCellsPerWavelength()*c.getMax()));
            cgw[ind].setTau(0.966/(this.farfieldFrequencyMax-this.farfieldFrequencyMin));
            cgw[ind].setT0(4.5*cgw[ind].getTau());
            cgw[ind].setFreq(this.operatingFrequency);
            
            for (int time_ind = 0; time_ind < this.timeSteps; time_ind++){
                insValue = (Math.exp(-((time[time_ind] - cgw[ind].getT0())/cgw[ind].getTau())* ((time[time_ind] - cgw[ind].getT0())/cgw[ind].getTau())))*(Math.cos(2*(Math.PI)*cgw[ind].getFreq()*(time[time_ind] - cgw[ind].getT0())));
		cgw[ind].set(insValue, time_ind);
            }
	}
        
        for (int ind = 0; ind < this.normalizedDerivativeGaussianWaveformCount ; ind++){
            if (ndgw[ind].getCellsPerWavelength() == 0){
		ndgw[ind].setCellsPerWavelength(Constants.CELLS_PER_WAVELENGTH);
            }

            ndgw[ind].setMaxFreq(speed/(ndgw[ind].getCellsPerWavelength()*c.getMax()));
            ndgw[ind].setTau((ndgw[ind].getCellsPerWavelength()*c.getMax()/(2*speed)));
            ndgw[ind].setT0(4.5*ndgw[ind].getTau());

            for (int time_ind = 0; time_ind < this.timeSteps; time_ind++){
                insValue = -(Math.sqrt(2*Math.E)*(time[time_ind] - ndgw[ind].getT0())/ndgw[ind].getTau())*Math.exp(-((time[time_ind] - ndgw[ind].getT0())/ndgw[ind].getTau())* ((time[time_ind] - ndgw[ind].getT0())/ndgw[ind].getTau()));
		ndgw[ind].set(insValue, time_ind);
            }
	}
    }
    
    private void initWaveForm(){
	int is,js,ks,ie,je,ke;
	for(int i = 0; i < VoltageSource.voltageSourceCount; i++){
            is = (int)Math.round((this.vs[i].getXmin() - ps.getXmin())/c.getDeltaX());
            js = (int)Math.round((this.vs[i].getYmin() - ps.getYmin())/c.getDeltaY());
            ks = (int)Math.round((this.vs[i].getZmin() - ps.getZmin())/c.getDeltaZ());
            ie = (int)Math.round((this.vs[i].getXmax() - ps.getXmin())/c.getDeltaX());
            je = (int)Math.round((this.vs[i].getYmax() - ps.getYmin())/c.getDeltaY());
            ke = (int)Math.round((this.vs[i].getZmax() - ps.getZmin())/c.getDeltaZ());
    
            this.vs[i].setIS(is);
            this.vs[i].setJS(js);
            this.vs[i].setKS(ks);
            this.vs[i].setIE(ie);
            this.vs[i].setJE(je);
            this.vs[i].setKE(ke);
		
            int n_fields = 0;
            double r_magnitude_factor = 0.0;
            double v_magnitude_factor = 0.0;
	
            if(this.vs[i].getDirection() == Constants.XN){
		n_fields = ie - is;
		r_magnitude_factor = (double)((1 + je - js) * (1 + ke - ks))/(ie - is); 
		v_magnitude_factor = (-1*vs[i].getMagnitude())/n_fields;
            }
            else if(vs[i].getDirection() == Constants.XP){
		n_fields = ie - is;
		r_magnitude_factor = (double)((1 + je - js) * (1 + ke - ks))/(ie - is); 
		v_magnitude_factor = (1*vs[i].getMagnitude())/n_fields;
            }
            else if(vs[i].getDirection() == Constants.YN){
		n_fields = je - js;
		r_magnitude_factor = (double)((1 + ie - is) * (1 + ke - ks))/(je - js); 
		v_magnitude_factor = (-1*vs[i].getMagnitude())/n_fields;
            }
            else if(vs[i].getDirection() == Constants.YP){
		n_fields = je - js;
		r_magnitude_factor = (double)((1 + ie - is) * (1 + ke - ks))/(je - js); 
		v_magnitude_factor = (1*vs[i].getMagnitude())/n_fields;
            }
            else if(vs[i].getDirection() == Constants.ZN){
		n_fields = ke - ks;
		r_magnitude_factor = (double)((1 + ie - is) * (1 + je - js))/(ke - ks); 
		v_magnitude_factor = (-1*vs[i].getMagnitude())/n_fields;
            }
            else if(vs[i].getDirection() == Constants.ZP){
		n_fields = ke - ks;
		r_magnitude_factor = (double)((1 + ie - is) * (1 + je - js))/(ke - ks); 
		v_magnitude_factor = (1*vs[i].getMagnitude())/n_fields;
            }

            vs[i].setResistancePerComponent(r_magnitude_factor * vs[i].getResistance());
            vs[i].initVoltagePerEfield(timeSteps);
            vs[i].initWaveform(timeSteps);

            if (this.gaussianWaveformCount != 0){
                for (int j=0;j<this.timeSteps;j++){
                    vs[i].setVoltagePerEfield(j, (v_magnitude_factor)*(w[0].get(j)));
                    vs[i].setWaveform(j, (v_magnitude_factor)*(w[0].get(j)*(n_fields)));
                }
            }
            
            if (this.sinusoidalWaveformCount != 0){
                for (int j=0;j<this.timeSteps;j++){
                    vs[i].setVoltagePerEfield(j, (v_magnitude_factor)*(sw[0].get(j)));
                    vs[i].setWaveform(j, (v_magnitude_factor)*(sw[0].get(j)*(n_fields)));
                }
            }
            
            if (this.cosineModulatedGaussianWaveformCount != 0){
                for (int j=0;j<this.timeSteps;j++){
                    vs[i].setVoltagePerEfield(j, (v_magnitude_factor)*(cgw[0].get(j)));
                    vs[i].setWaveform(j, (v_magnitude_factor)*(cgw[0].get(j)*(n_fields)));
                }
            }
            
            if (this.normalizedDerivativeGaussianWaveformCount != 0){
                for (int j=0;j<this.timeSteps;j++){
                    vs[i].setVoltagePerEfield(j, (v_magnitude_factor)*(ndgw[0].get(j)));
                    vs[i].setWaveform(j, (v_magnitude_factor)*(ndgw[0].get(j)*(n_fields)));
                }
            }
            
            vs[i].setTime(time);
            vs[i].setFreq(this.freqDomain.getFreqArray());
            vs[i].calFreqDomainValues();
            vs[i].saveWaveform("Waveform".concat(String.valueOf(i)));
	}
    }
    
    private void genCoefficients(){
        E.updatingCoefficients(ps);
        H.updatingCoefficients(ps);
        vs[0].updatingCoefficients(ps);
        vs[1].updatingCoefficients(ps);
    }
    
    private void initParam(){
        pml = new PMLmedia(ps,b,c);
        pbc = new PeriodicBoundary();
        pml.initAllCPMLboundary();
	
	if (eFieldProbCount > 0){
            eProbObject = new EFieldProb[eFieldProbCount];
		
            for (int ind = 0; ind < eFieldProbCount; ind++){  
                this.eProbObject[ind] = new EFieldProb(eFieldProbPos[ind],ps,c,timeSteps); 
            }
        }
        
        if (hFieldProbCount > 0){
            hProbObject = new HFieldProb[hFieldProbCount];
		
            for (int ind = 0; ind < hFieldProbCount; ind++){  
                this.hProbObject[ind] = new HFieldProb(hFieldProbPos[ind],ps,c,timeSteps); 
            } 
        }
        
        if (voltageProbCount > 0){
            for (int i=0; i < voltageProbCount; i++){
                this.voltageProbObject[i].initVoltageProb(ps, c, timeSteps);
            }       
        }

        if (currentProbCount > 0){
            for (int i=0; i < currentProbCount; i++){
                this.currentProbObject[i].initCurrentProb(ps, c, timeSteps);
            }       
        }
        ff.initFarFieldArrays(ps);
    }
    
    private void run(){
        System.out.println("Starting the time marching loop...");
        long startTimeMs, taskTimeMs;
        H.setBoundary(b);
        E.setBoundary(b);
        
        for (int timeIndex = 0; timeIndex < timeSteps; timeIndex++){
            startTimeMs = System.currentTimeMillis( );
            
            H.updateHField();
            
            pml.applyCPML2Hfield();
            
            for (int i=0; i<CurrentSource.currentSourceCount; i++){
                this.cs[i].updateCurrentSourceHfiled(timeIndex);
            }
            
            for (int i=0; i<hFieldProbCount; i++){
                this.hProbObject[i].CaptureHField(timeIndex);
            }
            
            for (int i=0; i<currentProbCount; i++){
                this.currentProbObject[i].CaptureCurrent(H,timeIndex);
            }
            
            
            E.updateEField();
            
            if(this.isActiveElem)
                pbc.updatePBCinXY();
            
            pml.applyCPML2Efield();
            
            for (int i=0; i<this.vs.length; i++){
                this.vs[i].updateVoltageSourceEfiled(timeIndex);
            }
            
            for (int i=0; i<eFieldProbCount; i++){
                this.eProbObject[i].CaptureEField(timeIndex);
            }
            
            for (int i=0; i<voltageProbCount; i++){
                this.voltageProbObject[i].CaptureVoltage(E,timeIndex);
            }
            
            ff.CalculateJandM(timeIndex);
            
            if(timeIndex%(timeSteps/100) == 0){
                taskTimeMs = System.currentTimeMillis( ) - startTimeMs;
                System.out.printf("**** %d %% Done! ****  %d Seconds Remaining\n" , (100*timeIndex/timeSteps) , taskTimeMs*(timeSteps-timeIndex)/1000 );
            }
        }
    }
    
    private void postProcessing(){
        postProcessingSparam();
        postProcessingNearField();
        postProcessingFarField();
    }
    
    private void postProcessingSparam(){        
        voltageProbObject[0].time2freq(freqDomain.getFreqArray());
        currentProbObject[0].time2freq(freqDomain.getFreqArray());
        voltageProbObject[1].time2freq(freqDomain.getFreqArray());
        currentProbObject[1].time2freq(freqDomain.getFreqArray());
        ports[0].setIProb(this.currentProbObject[0]);
        ports[1].setIProb(this.currentProbObject[1]);
        ports[0].setVProb(this.voltageProbObject[0]);
        ports[1].setVProb(this.voltageProbObject[1]);
        ports[0].calAB(freqDomain.getFreqArray());
        ports[1].calAB(freqDomain.getFreqArray());
        ports[0].calSparam(ports);
        ports[1].calSparam(ports);
        saveDataSparam();
    }
    
    private void postProcessingFarField(){
        ff.calTotalRadiatedPower();
        saveDataPower();
        
        ff.setupTheta(-180, 179, 360);
        ff.setupPhi(0, 0, 360);
        ff.calFarField2D();
        saveDataFarField("0Degree");
        
        ff.setupTheta(-180, 179, 360);
        ff.setupPhi(90, 90, 360);
        ff.calFarField2D();
        saveDataFarField("90Degree");
        
        ff.setupTheta(-180, 179, 360);
        ff.setupPhi(45, 45, 360);
        ff.calFarField2D();
        saveDataFarField("45Degree");
        
        this.postProcessing3DPattern();
    }
    
    public void postProcessing3DPattern(){
        ff.setupTheta(-90, 90, 181);
        ff.setupPhi(0, 180, 181);
        ff.calFarField3D();
        ff.saveFarFieldThetaPhi3DPattern("DPhi3D","DTheta3D","Phi3D","Theta3D");
    }
    
    private void postProcessingNearField(){
        saneDataNearField();
    }
    
    private void saneDataNearField(){
        eProbObject[0].saveEfiled();
    }
    
    private void saveDataSparam(){        
        voltageProbObject[0].saveTimeDomainValue();
        voltageProbObject[0].saveFreqDomainValue();
        voltageProbObject[0].saveTime();
        voltageProbObject[0].saveFreq();
        
        currentProbObject[0].saveTimeDomainValue();
        currentProbObject[0].saveFreqDomainValue();
        currentProbObject[0].saveTime();
        currentProbObject[0].saveFreq();
        
        voltageProbObject[1].saveTimeDomainValue();
        voltageProbObject[1].saveFreqDomainValue();
        voltageProbObject[1].saveTime();
        voltageProbObject[1].saveFreq();
        
        currentProbObject[1].saveTimeDomainValue();
        currentProbObject[1].saveFreqDomainValue();
        currentProbObject[1].saveTime();
        currentProbObject[1].saveFreq();
        
        ports[0].saveSParam("H");
        ports[1].saveSParam("V");
    }

    private void saveDataFarField(String cut){
        ff.saveFarFieldThetaPhiPattern(cut);
    }
    
    private void saveDataPower(){
        ff.saveTotalRadiatedPower();
    }
    
    public static void main(String[] args) {
        long start = System.currentTimeMillis();
        Pasim pas = new Pasim();
        pas.run();
        pas.postProcessing();
        System.out.println("Done!!");
        long totalTime = System.currentTimeMillis() - start;
        System.out.printf("The time taken : %d ms\n",totalTime);
    }
}