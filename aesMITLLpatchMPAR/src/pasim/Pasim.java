package pasim;

public class Pasim implements Runnable {
    private int timeSteps;
    private double courantFactor;
    private int brickCount, eFieldProbCount, hFieldProbCount, voltageProbCount, currentProbCount, portCount, voltageSourceCount,currentSourceCount;
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
    public Thread t;
    private double[][] geometryToSave;
    
    public Pasim(int[][] patchMatrix){
        initializeProblemSpace();
        defineGeometry(patchMatrix);
        defineSources();
        defineWaveform();
        defineOutput(protPosition);
        initMaterialGrid();
        initFDTDparams();
        initFDTDarrays();
        initSource();
        initWaveForm();
        genCoefficients();
        initParam();
        t = new Thread(this);
        t.start();
    }
    
    public void buildTheGeometryToTest(){
        geometryToSave = new double[ps.getNX()+1][ps.getNY()+1];
        for(int i=0;i<ps.getNZ();i++){
            Common.sumMetrix(geometryToSave, E.getEFieldArray(i));
        }
        Common.save2DArray(geometryToSave, "Test");
    }
    
    private void initializeProblemSpace(){
        timeSteps = 20000;
        courantFactor =0.9;
        operatingFrequency = 2.8e9; // Hz
        farfieldFrequencyMin = 2.5e9; // Hz
        farfieldFrequencyMax = 3.1e9; // Hz
        voltageSourceCount = 1;
        currentSourceCount = 0;
        double delta = 0.525e-3;
        c = new Cell(delta,delta,delta);
        createAllBoundaries();
        createAllMaterials();
    }
    
    private double elementDimension(double ratio){
        return ratio*wavelengthInFreeSpace();
    }
    
    private double wavelengthInFreeSpace(){
        return Constants.C/operatingFrequency;
    }
    
    private double wavelengthInMaterial(double epsR){
        return wavelengthInFreeSpace()/Math.sqrt(epsR);
    }
    
    private double resonantLength(double epsR){
        return 0.49*wavelengthInMaterial(epsR);
    }
    
    private void createAllBoundaries(){
        b = new Boundary();
        b.setAirBuffer(0,0,10);
        b.setCPML(false, false, true);
        b.setPBC(true, true, false);
        b.setCpmlCellNumber(0, 0, 8);
        b.setCpmlParam(3, 1.3, 7.0, 0.0, 0.05);
    }
    
    private void createAllMaterials(){
        Material.operatingFrequency = this.operatingFrequency;
        m = new Material[8];
        m[0] = Material.Vacuum(0);
        m[1] = Material.Air(1);
        m[2] = Material.PEC(2);
        m[3] = Material.PMC(3);
        m[4] = Material.RogersRTduroid5880(4);
        m[5] = Material.Copper(5);
        m[6] = Material.PerfectDielectric(6);
        m[7] = Material.RogersRTduroid4350(7);
    }
    
    private void defineGeometry(int[][] patchMatrix){
        int X = patchMatrix.length;
        int Y = patchMatrix[0].length;
        double cellX = c.getDeltaX();
        double cellY = c.getDeltaY();
        
        brickCount =3;
        for (int patchGeneX = 0; patchGeneX < X; patchGeneX++) {
            for (int patchGeneY = 0; patchGeneY < Y; patchGeneY++) {
                if(patchMatrix[patchGeneX][patchGeneY] == 1){
                     brickCount++;
                }
            }
        }
        
        bricks = new Brick[brickCount];
        
        double dim = elementDimension(0.5)/2;
        
        bricks[0] = CreateSimulationBox();
        
        bricks[1] = new Brick();
        bricks[1].Material(this.m[4]);
        bricks[1].Xmin(-dim);
        bricks[1].Ymin(-dim);
        bricks[1].Zmin(-1.575e-3);
        bricks[1].Xmax(dim);
        bricks[1].Ymax(dim);
        bricks[1].Zmax(0);
        
        bricks[2] = new Brick();
        bricks[2].Material(this.m[5]);
        bricks[2].Xmin(-dim);
        bricks[2].Ymin(-dim);
        bricks[2].Zmin(-1.575e-3);
        bricks[2].Xmax(dim);
        bricks[2].Ymax(dim);
        bricks[2].Zmax(-1.575e-3);
        
        
        int indexBrick = 2;
        for (int patchGeneX = 0; patchGeneX < X; patchGeneX++) {
            for (int patchGeneY = 0; patchGeneY < Y; patchGeneY++) {
                if(patchMatrix[patchGeneX][patchGeneY] == 1){
                    indexBrick++;
                    bricks[indexBrick] = new Brick();
                    bricks[indexBrick].Material(this.m[5]);
                    bricks[indexBrick].Xmin((-X/2 + patchGeneX)*cellX);
                    bricks[indexBrick].Ymin((-Y/2 + patchGeneY)*cellY);
                    bricks[indexBrick].Zmin(0);
                    bricks[indexBrick].Xmax((-X/2 + (patchGeneX + 1))*cellX);
                    bricks[indexBrick].Ymax((-Y/2 + (patchGeneY + 1))*cellY);
                    bricks[indexBrick].Zmax(0);
                }
            }
        }
    }
    
    private Brick CreateSimulationBox(){
        //double airGap = (3e8/this.farfieldFrequencyMin)/(2*Math.PI);
        double dim = (elementDimension(0.5)+c.getDeltaX())/2;
        double airGap = c.getDeltaZ();
        Brick b = new Brick();
        b.Material(this.m[1]);
        b.Xmin(-dim - c.getDeltaX());
        b.Ymin(-dim - c.getDeltaY());
        b.Zmin(-1.575e-3 - airGap);
        b.Xmax(dim + c.getDeltaX());
        b.Ymax(dim + c.getDeltaY());
        b.Zmax(0 + airGap);
        return b;
    }
    
    private void defineSources(){
        double leng = 0;
        for(int i=0;i<protPosition.length;i++){
            leng = leng + protPosition[i]*Math.pow(2, i);
        }
        
	if (voltageSourceCount > 0){
            vs = new VoltageSource[voltageSourceCount];
            for (int i = 0; i<voltageSourceCount;i++){
                vs[i] = new VoltageSource(0, leng*c.getDeltaY(), -1.575e-3, 0, leng*c.getDeltaY(), 0, c);
                vs[i].setDirection(Constants.ZP);
                vs[i].setResistance(50.0);
                vs[i].setMagnitude(1.0);
            }
	}

	if (currentSourceCount > 0){
            cs = new CurrentSource[currentSourceCount];
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
        eFieldProbCount  = 0;
        hFieldProbCount  = 0;
        voltageProbCount  = 2;
        currentProbCount  = 2;
        portCount = 2;
        double dim = elementDimension(0.5)/2;
        
        currentProbPos = new Position[2];
        voltageProbPos = new Position[2];
        eFieldProbPos = new Position[1];
        hFieldProbPos = new Position[1];
        
        double leng = 0;
        for(int i=0;i<protPosition.length;i++){
            leng = leng + protPosition[i]*Math.pow(2, i);
        }
        
        currentProbPos[0] = new Position(0, leng*c.getDeltaY(), -1e-3);
        voltageProbPos[0] = new Position(0, leng*c.getDeltaY(), 0);
        currentProbPos[1] = new Position(leng*c.getDeltaX(), 0, -1e-3);
        voltageProbPos[1] = new Position(leng*c.getDeltaX(), 0, 0);
        eFieldProbPos[0] = new Position(dim, dim, 1e-3);
        hFieldProbPos[0] = new Position(dim, dim, 1e-3);
        
        voltageProbObject = new VoltageProb[voltageProbCount];
        currentProbObject = new CurrentProb[currentProbCount];
        ports = new Port[portCount];
        freqDomain = new FrequencyDomain(2e9,3.5e9,1501);
        ff = new FarField(2.5e9,3.1e9,7,c);
        ff.setCellCountOuterBoundary(13);
        
        if(voltageProbObject.length > 0){
            for (int i=0;i<voltageProbObject.length;i++){
                voltageProbObject[i] = new VoltageProb(voltageProbPos[i],Constants.ZP,-1.575e-3);
            }
        }
        
        if(currentProbObject.length > 0){
            for (int i=0;i<currentProbObject.length;i++){
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
            cgw[ind].setFreq(2.8e9);
            
            for (int time_ind = 0; time_ind < this.timeSteps; time_ind++){
                insValue = (Math.exp(-((time[time_ind] - cgw[ind].getT0())/cgw[ind].getTau())* ((time[time_ind] - cgw[ind].getT0())/cgw[ind].getTau())))*(Math.sin(2*(Math.PI)*cgw[ind].getFreq()*(time[time_ind] - cgw[ind].getT0())));
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
	for(int i = 0; i < voltageSourceCount; i++){
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
            //vs[i].saveWaveform("Waveform".concat(String.valueOf(i)));
	}
    }
    
    private void genCoefficients(){
        E.updatingCoefficients(ps);
        H.updatingCoefficients(ps);
        vs[0].updatingCoefficients(ps,E,H);
    }
    
    private void initParam(){
        pml = new PMLmedia(ps,b,c);
        pbc = new PeriodicBoundary(ps);
        pml.initAllCPMLboundary(E,H);
	
	if (eFieldProbCount > 0){
            eProbObject = new EFieldProb[eFieldProbCount];
		
            for (int ind = 0; ind < eFieldProbCount; ind++){  
                this.eProbObject[ind] = new EFieldProb(eFieldProbPos[ind],ps,c,timeSteps); 
            }
        }
        
        if (hFieldProbCount > 0){
            hProbObject = new HFieldProb[hFieldProbCount];
		
            for (int ind = 0; ind < hProbObject.length; ind++){  
                this.hProbObject[ind] = new HFieldProb(hFieldProbPos[ind],ps,c,timeSteps); 
            } 
        }
        
        if (voltageProbObject.length > 0){
            for (int i=0; i < voltageProbObject.length; i++){
                this.voltageProbObject[i].initVoltageProb(ps, c, timeSteps);
            }       
        }

        if (currentProbObject.length > 0){
            for (int i=0; i < currentProbObject.length; i++){
                this.currentProbObject[i].initCurrentProb(ps, c, timeSteps);
            }       
        }
        ff.initFarFieldArrays(ps,b);
    }
    
    private void updateHfieldFromCurrentSources(int timeIndex, HField H){
        for (int i=0; i<currentSourceCount; i++){
            this.cs[i].updateCurrentSourceHfiled(timeIndex, H);
        }
    }
    
    private void samplingHfield(int timeIndex, HField H){
        for (int i=0; i<hFieldProbCount; i++){
            this.hProbObject[i].CaptureHField(timeIndex,H);
        }
    }
    
    private void samplingCurrentAtPorts(int timeIndex, HField H){
        for (int i=0; i<currentProbObject.length; i++){
            this.currentProbObject[i].CaptureCurrent(H,timeIndex);
        }
    }
    
    private void updateEfieldFromVoltageSources(int timeIndex, EField E){
        for (int i=0; i<this.vs.length; i++){
            this.vs[i].updateVoltageSourceEfiled(timeIndex,E);
        }
    }
    
    private void samplingEfield(int timeIndex, EField E){
        for (int i=0; i<eFieldProbCount; i++){
            this.eProbObject[i].CaptureEField(timeIndex,E);
        }
    }
    
    private void samplingVoltageAtPorts(int timeIndex, EField E){
        for (int i=0; i<voltageProbObject.length; i++){
            this.voltageProbObject[i].CaptureVoltage(E,timeIndex);
        }
    }
    
    private void timeMarchingLoop(){
        for (int timeIndex = 0; timeIndex < timeSteps; timeIndex++){
            H.updateHField(E);
            
            pml.applyCPML2Hfield(E,H);
            
            //this.updateHfieldFromCurrentSources(timeIndex,H);
            
            //this.samplingHfield(timeIndex, H);
            
            this.samplingCurrentAtPorts(timeIndex, H);
            
            E.updateEField(H);
            
            pbc.updatePBCinXY(E,H);
            
            pml.applyCPML2Efield(E,H);
            
            this.updateEfieldFromVoltageSources(timeIndex, E);
            
            this.samplingEfield(timeIndex, E);
            
            this.samplingVoltageAtPorts(timeIndex, E);
            
            //ff.CalculateJandM(timeIndex,E,H);
            //if(timeIndex%5==0)
            //    buildTheGeometryToTest();
        }
    }
    
    public void run(){
        H.setBoundary(b);
        E.setBoundary(b);
        ff.setFarFieldDisable(true);
        timeMarchingLoop();
        this.postProcessing();
    }
    
    private void postProcessing(){
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
        ff.calTotalRadiatedPower();
        ff.setupTheta(-180, 179, 360);
        ff.setupPhi(0, 0, 360);
        ff.calFarField();
    }
    
    public double calcFitness(){
        double[] fitness = new double[4];
        double totalFitness = 0.0;
        
        fitness[0] = 10*calcFitnessSparam(2.78e9, 2.82e9,"S11",-13);
        fitness[1] = 1*calcFitnessSparam(2.74e9, 2.86e9,"S11",-0.0001);
        fitness[2] = 1*calcFitnessSparam(2.70e9, 2.90e9,"S21",-35);
        fitness[3] = 10*calcFitnessSparam(2.76e9, 2.84e9,"S21",-50);
        
        for(int index = 0; index < fitness.length;index++){
            totalFitness = totalFitness + fitness[index];
        }
        
        return totalFitness/22;
    }
    
    /*
    public double calcFitness(){
        double fitness1 = 0.0;
        double[] designningBandwidth1 = freqDomain.getSubFreqArray(2.78e9, 2.82e9);
        double[][] s1=this.ports[0].getSparamAbs(designningBandwidth1, this.portCount);
        for(int freqIndex=0; freqIndex<designningBandwidth1.length;freqIndex++){
            if(20*Math.log10(Math.abs(s1[0][freqIndex]))<=(-13)){
                fitness1 = fitness1 + 1.0;
            }
        }
        fitness1 = fitness1/(2*designningBandwidth1.length);
        
        double fitness2 = 0.0;
        double[] designningBandwidth2 = freqDomain.getSubFreqArray(2.7e9, 2.9e9);
        double[][] s2=this.ports[0].getSparamAbs(designningBandwidth2, this.portCount);
        for(int freqIndex=0; freqIndex<designningBandwidth2.length;freqIndex++){
            if(20*Math.log10(Math.abs(s2[1][freqIndex]))<=(-60)){
                fitness2 = fitness2 + 1.0;
            }
        }
        fitness2 = fitness2/(2*designningBandwidth2.length);
        
        double fitness3 = 0.0;
        double[] designningBandwidth3 = freqDomain.getSubFreqArray(2.5e9, 3.1e9);
        double[][] s3=this.ports[0].getSparamAbs(designningBandwidth3, this.portCount);
        for(int freqIndex=0; freqIndex<designningBandwidth3.length;freqIndex++){
            if(20*Math.log10(Math.abs(s3[0][freqIndex]))<=0){
                fitness3 = fitness3 + 1.0;
            }
        }
        fitness3 = fitness3/(2*designningBandwidth3.length);
        
        double fitness = (fitness1 + fitness2 + fitness3)/3;
        return fitness;
    }
    */
    
    public double calcFitnessSparam(double first, double last, String Label,double threshold){
        double fitness = 0.0;
        double[][] s;
        double[] designningBandwidth = freqDomain.getSubFreqArray(first, last);
        double calcValue = 0.0;
        
        if(Label.equals("S11")){
            s=this.ports[0].getSparamAbs(designningBandwidth, this.portCount);
            for(int freqIndex=0; freqIndex<designningBandwidth.length;freqIndex++){
                calcValue = 20*Math.log10(Math.abs(s[0][freqIndex]));
                if(calcValue<=threshold){
                    fitness = fitness + 1.0;
                }
                else if(calcValue<=0){
                    fitness = fitness + calcValue/threshold;
                }
            }
        }
        
        else if(Label.equals("S21")){
            s=this.ports[0].getSparamAbs(designningBandwidth, this.portCount);
            for(int freqIndex=0; freqIndex<designningBandwidth.length;freqIndex++){
                calcValue = 20*Math.log10(Math.abs(s[1][freqIndex]));
                if(calcValue<=threshold){
                    fitness = fitness + 1.0;
                }
                else if(calcValue<=0){
                    fitness = fitness + calcValue/threshold;
                }
            }
        }
        
        else{
            System.out.println("ERROR : Not a valide label.. Exiting the simulation...");
            System.exit(0);
        }
        
        return fitness/(designningBandwidth.length);
    }
    
    public void saveSParam(int i){
        String H=String.valueOf(i).concat("H");
        String V=String.valueOf(i).concat("V");
        ports[0].saveSParam(H);
        ports[1].saveSParam(V);
    }

    public void saveData(){
        //eProbObject[0].saveEfiled();
        
        //voltageProbObject[0].saveTimeDomainValue();
        //voltageProbObject[0].saveFreqDomainValue();
        //voltageProbObject[0].saveTime();
        //voltageProbObject[0].saveFreq();
        
        //currentProbObject[0].saveTimeDomainValue();
        //currentProbObject[0].saveFreqDomainValue();
        //currentProbObject[0].saveTime();
        //currentProbObject[0].saveFreq();
        
        //voltageProbObject[1].saveTimeDomainValue();
        //voltageProbObject[1].saveFreqDomainValue();
        //voltageProbObject[1].saveTime();
        //voltageProbObject[1].saveFreq();
        
        //currentProbObject[1].saveTimeDomainValue();
        //currentProbObject[1].saveFreqDomainValue();
        currentProbObject[1].saveTime();
        currentProbObject[1].saveFreq();
        
        ports[0].saveSParam("H_Very_Best");
        ports[1].saveSParam("V_Very_Best");
        //ff.saveTotalRadiatedPower();
        //ff.saveFarFieldThetaPhiPattern();

    }
}