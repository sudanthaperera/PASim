package pasim;

public class East implements Runnable {
    private int timeSteps;
    private int elementCount;
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
    
    public East(Cell c){
        initializeProblemSpace(c);
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
    
    public void start(){
        t = new Thread(this);
        t.setName("Finite-By-Infinite");
        t.start();
    }
    
    public FarField getFF(){
        return ff;
    }
    
    private void initializeProblemSpace(Cell c){
        timeSteps = 10000;
        courantFactor =0.9;
        operatingFrequency = 2.8e9; // Hz
        farfieldFrequencyMin = 2.5e9; // Hz
        farfieldFrequencyMax = 3.1e9; // Hz
        voltageSourceCount = 2*elementCount;
        currentSourceCount = 0;
        elementCount = 2;
        this.c = c;
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
        b.setAirBuffer(10,0,10);
        b.setCPML(true, false, true);
        b.setPBC(false, true, false);
        b.setCpmlCellNumber(8, 0, 8);
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
    
    private void defineGeometry(){
        brickCount =3 + elementCount;
        bricks = new Brick[brickCount];
        
        bricks[0] = CreateSimulationBox();
        
        bricks[1] = new Brick();
        bricks[1].Material(this.m[4]);
        bricks[1].Xmin(-elementCount*27.0e-3);
        bricks[1].Ymin(-27.0e-3);
        bricks[1].Zmin(-1.5e-3);
        bricks[1].Xmax(elementCount*27.0e-3);
        bricks[1].Ymax(27.0e-3);
        bricks[1].Zmax(0);
        
        bricks[2] = new Brick();
        bricks[2].Material(this.m[5]);
        bricks[2].Xmin(-elementCount*27.0e-3);
        bricks[2].Ymin(-27.0e-3);
        bricks[2].Zmin(-3.575e-3);
        bricks[2].Xmax(elementCount*27.0e-3);
        bricks[2].Ymax(27.0e-3);
        bricks[2].Zmax(-1.5e-3);
        
        bricks[3] = new Brick();
        bricks[3].Material(this.m[5]);
        bricks[3].Xmin(-17.0e-3);
        bricks[3].Ymin(-17.0e-3);
        bricks[3].Zmin(0);
        bricks[3].Xmax(17.0e-3);
        bricks[3].Ymax(17.0e-3);
        bricks[3].Zmax(0);
        
        for(int i=1;i<(elementCount+1)/2;i++){
            bricks[2+i*2] = new Brick(bricks[3]);
            bricks[2+i*2].move(-2*i*27.0e-3, 0, 0);
                
            bricks[3+i*2] = new Brick(bricks[3]);
            bricks[3+i*2].move(+2*i*27.0e-3, 0, 0);
        }
    }
    
    private Brick CreateSimulationBox(){
        double dimX = c.getDeltaX()*((int)Math.round(((elementDimension(0.5)+c.getDeltaX())/2)/c.getDeltaX()));
        double dimZ = c.getDeltaZ()*((int)Math.round(((elementDimension(0.5)+c.getDeltaZ())/2)/c.getDeltaZ()));

        Brick b = new Brick();
        b.Material(this.m[1]);
        b.Xmin(-dimX - (27.0e-3)*elementCount);
        b.Ymin(-27.0e-3);
        b.Zmin(-dimZ - 1.5e-3);
        b.Xmax(dimX + (27.0e-3)*elementCount);
        b.Ymax(27.0e-3);
        b.Zmax(dimZ);
        return b;
    }
    
    private void defineSources(){
        vs = new VoltageSource[2*elementCount];
        
        vs[0] = this.createVoltageSouraceH(0, 0, 1);
        vs[1] = this.createVoltageSouraceV(0, 0, 0);
        
        for(int i=1;i<(elementCount+1)/2;i++){            
            vs[2+4*(i-1)] = this.createVoltageSouraceH(-2*i*27.0e-3, 0, 1);
            vs[3+4*(i-1)] = this.createVoltageSouraceV(-2*i*27.0e-3, 0, 0);
        
            vs[4+4*(i-1)] = this.createVoltageSouraceH(2*i*27.0e-3, 0, 1);
            vs[5+4*(i-1)] = this.createVoltageSouraceV(2*i*27.0e-3, 0, 0);
        }
        /*
        vs[2] = this.createVoltageSouraceH(-2*27.0e-3, 0, 1);
        vs[3] = this.createVoltageSouraceV(-2*27.0e-3, 0, 0);
        
        vs[4] = this.createVoltageSouraceH(2*27.0e-3, 0, 1);
        vs[5] = this.createVoltageSouraceV(2*27.0e-3, 0, 0);
        
        vs[6] = this.createVoltageSouraceH(-4*27.0e-3, 0, 1);
        vs[7] = this.createVoltageSouraceV(-4*27.0e-3, 0, 0);
        
        vs[8] = this.createVoltageSouraceH(4*27.0e-3, 0, 1);
        vs[9] = this.createVoltageSouraceV(4*27.0e-3, 0, 0);
        
        vs[10] = this.createVoltageSouraceH(-6*27.0e-3, 0, 1);
        vs[11] = this.createVoltageSouraceV(-6*27.0e-3, 0, 0);
        
        vs[12] = this.createVoltageSouraceH(6*27.0e-3, 0, 1);
        vs[13] = this.createVoltageSouraceV(6*27.0e-3, 0, 0);
        
        vs[14] = this.createVoltageSouraceH(-8*27.0e-3, 0, 1);
        vs[15] = this.createVoltageSouraceV(-8*27.0e-3, 0, 0);
        
        vs[16] = this.createVoltageSouraceH(8*27.0e-3, 0, 1);
        vs[17] = this.createVoltageSouraceV(8*27.0e-3, 0, 0);
        
        vs[18] = this.createVoltageSouraceH(-10*27.0e-3, 0, 1);
        vs[19] = this.createVoltageSouraceV(-10*27.0e-3, 0, 0);
        
        vs[20] = this.createVoltageSouraceH(10*27.0e-3, 0, 1);
        vs[21] = this.createVoltageSouraceV(10*27.0e-3, 0, 0);
        
        vs[22] = this.createVoltageSouraceH(-12*27.0e-3, 0, 1);
        vs[23] = this.createVoltageSouraceV(-12*27.0e-3, 0, 0);
        
        vs[24] = this.createVoltageSouraceH(12*27.0e-3, 0, 1);
        vs[25] = this.createVoltageSouraceV(12*27.0e-3, 0, 0);
        
        vs[26] = this.createVoltageSouraceH(-14*27.0e-3, 0, 1);
        vs[27] = this.createVoltageSouraceV(-14*27.0e-3, 0, 0);
        
        vs[28] = this.createVoltageSouraceH(14*27.0e-3, 0, 1);
        vs[29] = this.createVoltageSouraceV(14*27.0e-3, 0, 0);
        
        vs[30] = this.createVoltageSouraceH(-16*27.0e-3, 0, 1);
        vs[31] = this.createVoltageSouraceV(-16*27.0e-3, 0, 0);
        
        vs[32] = this.createVoltageSouraceH(16*27.0e-3, 0, 1);
        vs[33] = this.createVoltageSouraceV(16*27.0e-3, 0, 0);
        
        vs[34] = this.createVoltageSouraceH(-18*27.0e-3, 0, 1);
        vs[35] = this.createVoltageSouraceV(-18*27.0e-3, 0, 0);
        
        vs[36] = this.createVoltageSouraceH(18*27.0e-3, 0, 1);
        vs[37] = this.createVoltageSouraceV(18*27.0e-3, 0, 0);
                */
    }
    
    private VoltageSource createVoltageSouraceH(double shiftX, double shiftY, double mag){
        VoltageSource vsource = new VoltageSource(0, 7.0e-3, -1.5e-3, 0, 7.0e-3, 0, c);
        vsource.move(shiftX, shiftY, 0);
        vsource.setDirection(Constants.ZP);
        vsource.setResistance(50.0);
        vsource.setMagnitude(mag);
        
        return vsource;
    }
    
    private VoltageSource createVoltageSouraceV(double shiftX, double shiftY, double mag){
        VoltageSource vsource = new VoltageSource(7.0e-3, 0, -1.5e-3, 7.0e-3, 0, 0, c);
        vsource.move(shiftX, shiftY, 0);
        vsource.setDirection(Constants.ZP);
        vsource.setResistance(50.0);
        vsource.setMagnitude(mag);
        
        return vsource;
    }
    
    private void defineWaveform(){
        gaussianWaveformCount = 2;
        sinusoidalWaveformCount = 0;
        cosineModulatedGaussianWaveformCount = 0;
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
        portCount = 2*elementCount;
        
        currentProbPos = new Position[2];
        voltageProbPos = new Position[2];
        eFieldProbPos = new Position[1];
        hFieldProbPos = new Position[1];
        
        currentProbPos[0] = new Position(0, 7.0e-3, -1e-3);
        voltageProbPos[0] = new Position(0, 7.0e-3, 0);
        currentProbPos[1] = new Position(7.0e-3, 0, -1e-3);
        voltageProbPos[1] = new Position(7.0e-3, 0, 0);
        eFieldProbPos[0] = new Position(20e-3, 20e-3, 0e-3);
        hFieldProbPos[0] = new Position(20e-3, 20e-3, 0e-3);
        
        voltageProbObject = new VoltageProb[voltageProbCount];
        currentProbObject = new CurrentProb[currentProbCount];
        ports = new Port[portCount];
        freqDomain = new FrequencyDomain(2e9,3.5e9,1501);
        ff = new FarField(2.7e9,2.9e9,3,c);
        ff.setCellCountOuterBoundary(13);
        
        if(voltageProbCount > 0){
            for (int i=0;i<voltageProbCount;i++){
                voltageProbObject[i] = new VoltageProb(voltageProbPos[i],Constants.ZP,-1.5e-3);
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
        ps.initMaterialGrid(false);
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
        for(int i=0;i<vs.length;i++){
            vs[i].updatingCoefficients(ps,E,H);
        }
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
        System.out.println("Starting the time marching loop...");
        long startTimeMs, taskTimeMs;
        for (int timeIndex = 0; timeIndex < timeSteps; timeIndex++){
            startTimeMs = System.currentTimeMillis( );
            
            H.updateHField(E);
            
            pml.applyCPML2Hfield(E,H);
            
            this.samplingHfield(timeIndex, H);
            
            this.samplingCurrentAtPorts(timeIndex, H);
            
            E.updateEField(H);
            
            pbc.updatePBCinX(E,H);
            
            pml.applyCPML2Efield(E,H);
            
            this.updateEfieldFromVoltageSources(timeIndex, E);
            
            this.samplingEfield(timeIndex, E);
            
            this.samplingVoltageAtPorts(timeIndex, E);
            
            ff.CalculateJandM(timeIndex,E,H,b);
            
            if(timeIndex%(timeSteps/100) == 0){
                taskTimeMs = System.currentTimeMillis( ) - startTimeMs;
                System.out.printf("****Finite-By-Infinite - %d %% Done! ****  %d Seconds Remaining\n" , (100*timeIndex/timeSteps) , taskTimeMs*(timeSteps-timeIndex)/1000 );
            }
        }
    }
    
    public void run(){
        H.setBoundary(b);
        E.setBoundary(b);
        ff.setFarFieldDisable(false);
        timeMarchingLoop();
        E=null;
        H=null;        
        ps = null;
        pml = null;
        c = null;
        b = null;
        m = null;
        bricks = null;
        vs = null;
        w = null;
        sw = null;
        cgw = null;
        ndgw = null;
        currentProbPos = null;
        voltageProbPos = null;
        eFieldProbPos = null;
        hFieldProbPos = null;
        voltageProbObject = null;
        currentProbObject = null;
        hProbObject = null;
        eProbObject = null;
        ports = null;
        freqDomain = null;
        pbc = null;
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
    
    public void saveSParam(int i){
        String H=String.valueOf(i).concat("H");
        String V=String.valueOf(i).concat("V");
        ports[0].saveSParam(H);
        ports[1].saveSParam(V);
    }

    public void saveData(){
        eProbObject[0].saveEfiled();
        
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
        ff.saveTotalRadiatedPower();
        ff.saveFarFieldThetaPhiPattern();

    }
}