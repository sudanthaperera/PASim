package pasim;


/**
 *
 * @author Sudantha
 */
public class Pasim {
    private int timeSteps;
    private double courantFactor;
    private int brickCount, eFieldProbCount, hFieldProbCount, voltageProbCount, currentProbCount, portCount, gaussianWaveformCount;
    private double speed;
    private double eta0;
    private double[] time;
    
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
    private Position[] currentProbPos, voltageProbPos, eFieldProbPos, hFieldProbPos;
    private VoltageProb[] voltageProbObject;
    private CurrentProb[] currentProbObject;
    private HFieldProb[] hProbObject;
    private EFieldProb[] eProbObject;
    private Port ports[];
    private FarField ff;
    private FrequencyDomain freqDomain;
    
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
        System.out.println("Defining the problem space parameters...");
        timeSteps = 10000;
        courantFactor =0.9;
        double operatingFrequency = 2.8e9;
        double velocityOfLightInVac = 3e8;
        double radious = 11.7*(velocityOfLightInVac/operatingFrequency);
        VoltageSource.voltageSourceCount = 2;
        CurrentSource.currentSourceCount = 0;
        c = new Cell(0.525e-3,0.525e-3/radious,0.525e-3,radious);
        b = new Boundary();
        b.setAirBuffer(10,10,10,10,10,10);
        b.setBoundaryType(true, true, true, true, true, true);
        b.setCpmlCellNumber(8, 8, 8, 8, 8, 8);
        b.setCpmlParam(3, 1.3, 7.0, 0.0, 0.05);
        
        m = new Material[5];
        m[0] = new Material(0, 1.0, 1.0, 0, 0); // Vac
        m[1] = new Material(1, 1.0, 1.0, 0, 0); // air
        m[2] = new Material(2, 1.0, 1.0, 1.0e10, 0); // PEC
        m[3] = new Material(3, 1.0, 1.0, 0, 1.0e10); // PMC
        m[4] = new Material(4, 2.2, 1.0, 0, 0); // Dielectric
    }
    
    private void defineGeometry(){
        System.out.println("Defining the problem geometry...");
        
        brickCount =4;
        bricks = new Brick[brickCount];
        
        bricks[0] = new Brick(c.getR()-1.575e-3, -39.375e-3/c.getR(), -39.375e-3, c.getR() + 0, 39.375e-3/c.getR(), 39.375e-3, this.m[1]); //Air
        bricks[1] = new Brick(c.getR()-1.575e-3, -39.375e-3/c.getR(), -39.375e-3, c.getR() + 0, 39.375e-3/c.getR(), 39.375e-3, this.m[4]); //Dielectric
        bricks[2] = new Brick(c.getR()-1.575e-3, -39.375e-3/c.getR(), -39.375e-3, c.getR() + -1.575e-3, 39.375e-3/c.getR(), 39.375e-3, this.m[2]); //Ground
        bricks[3] = new Brick(c.getR()+0, -17.775e-3/c.getR(), -17.775e-3, c.getR() + 0, 17.775e-3/c.getR(), 17.775e-3, this.m[2]); //Patch
    }
    
    private void defineSources(){
        System.out.println("Defining sources and lumped element components...");
	if (VoltageSource.voltageSourceCount > 0){
            vs = new VoltageSource[VoltageSource.voltageSourceCount];
            for (int i = 0; i<VoltageSource.voltageSourceCount;i++){
                vs[i] = new VoltageSource(c.getR() -1.575e-3, 5.25e-3/c.getR(),0 , c.getR() + 0, 5.25e-3/c.getR(), 0, c);
                vs[i].setDirection(Constants.RP);
                vs[i].setResistance(50.0);
                vs[i].setMagnitude(1.0);
            }
            vs[1].setMagnitude(0.0);
	}

	if (CurrentSource.currentSourceCount > 0){
            cs = new CurrentSource[CurrentSource.currentSourceCount];
	}
    }
    
    private void defineWaveform(){
        System.out.println("Defining source waveform types and parameters...");
        gaussianWaveformCount = 2;
        w = new Gaussian[gaussianWaveformCount];
        
        w[0] = new Gaussian(20,1.0,0.0,timeSteps);
        w[1] = new Gaussian(20,1.0,0.0,timeSteps);
    }
    
    private void defineOutput(){
        System.out.println("Defining output parameters...");

        eFieldProbCount  = 1;
        hFieldProbCount  = 0;
        voltageProbCount  = 2;
        currentProbCount  = 2;
        portCount = 2;
        
        currentProbPos = new Position[2];
        voltageProbPos = new Position[2];
        eFieldProbPos = new Position[1];
        hFieldProbPos = new Position[1];
        
        currentProbPos[0] = new Position(c.getR()-1e-3, 5.25e-3/c.getR(), 0);
        voltageProbPos[0] = new Position(c.getR(), 5.25e-3/c.getR(), 0);
        currentProbPos[1] = new Position(c.getR()-1e-3, 0, 5.25e-3);
        voltageProbPos[1] = new Position(c.getR(), 0, 5.25e-3);
        eFieldProbPos[0] = new Position(c.getR() + 1e-3, 38.375e-3/c.getR(), 38.375e-3);
        hFieldProbPos[0] = new Position(c.getR() + 1e-3, 38.375e-3/c.getR(), 38.375e-3);
        
        voltageProbObject = new VoltageProb[voltageProbCount];
        currentProbObject = new CurrentProb[currentProbCount];
        ports = new Port[portCount];
        freqDomain = new FrequencyDomain(1e9,5e9,40001);
        ff = new FarField(2.5e9,3.1e9,7);
        ff.setCellCountOuterBoundary(13);
        
        if(voltageProbCount > 0){
            for (int i=0;i<voltageProbCount;i++){
                voltageProbObject[i] = new VoltageProb(voltageProbPos[i],Constants.RP,-1.575e-3);
            }
        }
        
        if(currentProbCount > 0){
            for (int i=0;i<currentProbCount;i++){
                currentProbObject[i] = new CurrentProb(currentProbPos[i],Constants.RP);
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
	System.out.printf("Number of cells in R direction = %d\n" , ps.getNR());
	System.out.printf("Number of cells in Alpha direction = %d\n" , ps.getNA());
	System.out.printf("Number of cells in Z direction = %d\n" , ps.getNZ());
	System.out.printf("Number of cells = %d \n" , ps.getNR()*ps.getNA()*ps.getNZ());
        
        System.out.println("Initializing FDTD material grid...");
        ps.initMaterialGrid();
    }
    
    private void initFDTDparams(){
        speed = 1/Math.sqrt(Constants.MU0*Constants.EPS0);
        eta0 = Math.sqrt(Constants.MU0/Constants.EPS0);
        c.setDeltaT(courantFactor,speed,(vs[0].getAmax() - vs[0].getAmin())/2);
        time = new double[this.timeSteps];
        
	for(int index = 0; index < this.timeSteps; index++){
            time[index] = (0.5 + index)*c.getDeltaT();
	}        
    }

    private void initFDTDarrays(){
        H = new HField(this.ps.getNR(),this.ps.getNA(),this.ps.getNZ(),c);
        E = new EField(this.ps.getNR(),this.ps.getNA(),this.ps.getNZ(),c);
        
        H.setAllR(0.0);
        H.setAllA(0.0);
        H.setAllZ(0.0);
        
        E.setAllR(0.0);
        E.setAllA(0.0);
        E.setAllZ(0.0);
    }
    
    private void initSource(){
        double insValue;
        double[] R4Prob = new double[2];
        R4Prob[0] = (vs[0].getAmax() - vs[0].getAmin())/2;
        R4Prob[1] = (vs[1].getAmax() - vs[1].getAmin())/2;
        
	for (int ind = 0; ind < this.gaussianWaveformCount ; ind++){
            if (w[ind].getCellsPerWavelength() == 0){
		w[ind].setCellsPerWavelength(Constants.CELLS_PER_WAVELENGTH);
            }

            w[ind].setMaxFreq(speed/(w[ind].getCellsPerWavelength()*c.getMax(R4Prob[ind])));
            w[ind].setTau((w[ind].getCellsPerWavelength()*c.getMax(R4Prob[ind])/(2*speed)));
            w[ind].setT0(4.5*w[ind].getTau());

            for (int time_ind = 0; time_ind < this.timeSteps; time_ind++){
                insValue = Math.exp(-((time[time_ind] - w[ind].getT0())/w[ind].getTau())* ((time[time_ind] - w[ind].getT0())/w[ind].getTau()));
		w[ind].set(insValue, time_ind);
            }
	}
    }
    
    private void initWaveForm(){
	int is,js,ks,ie,je,ke;
	for(int i = 0; i < VoltageSource.voltageSourceCount; i++){
            is = (int)Math.round((this.vs[i].getRmin() - ps.getRmin())/c.getDeltaR());
            js = (int)Math.round((this.vs[i].getAmin() - ps.getAmin())/c.getDeltaA());
            ks = (int)Math.round((this.vs[i].getZmin() - ps.getZmin())/c.getDeltaZ());
            ie = (int)Math.round((this.vs[i].getRmax() - ps.getRmin())/c.getDeltaR());
            je = (int)Math.round((this.vs[i].getAmax() - ps.getAmin())/c.getDeltaA());
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
	
            if(this.vs[i].getDirection() == Constants.RN){
		n_fields = ie - is;
		r_magnitude_factor = (double)((1 + je - js) * (1 + ke - ks))/(ie - is);
		v_magnitude_factor = (-1*vs[i].getMagnitude())/n_fields;
            }
            else if(vs[i].getDirection() == Constants.RP){
		n_fields = ie - is;
		r_magnitude_factor = (double)((1 + je - js) * (1 + ke - ks))/(ie - is);
		v_magnitude_factor = (1*vs[i].getMagnitude())/n_fields;
            }
            else if(vs[i].getDirection() == Constants.AN){
		n_fields = je - js;
		r_magnitude_factor = (double)((1 + ie - is) * (1 + ke - ks))/(je - js);
		v_magnitude_factor = (-1*vs[i].getMagnitude())/n_fields;
            }
            else if(vs[i].getDirection() == Constants.AP){
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

            for (int j=0;j<this.timeSteps;j++){
		vs[i].setVoltagePerEfield(j, (v_magnitude_factor)*(w[0].get(j)));
		vs[i].setWaveform(j, (v_magnitude_factor)*(w[0].get(j)*(n_fields)));
            }
	}
    }
    
    private void genCoefficients(){
        E.updatingCoefficients(ps);
        H.updatingCoefficients(ps);
        vs[0].updatingCoefficients(ps);
    }
    
    private void initParam(){
        System.out.println("initializing the CPML parameters...");
        pml = new PMLmedia(ps,b,c);
        pml.initAllCPMLboundary();
        
	System.out.println("initializing the output parameters...");
	
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
            
            if(timeIndex%(timeSteps/10) == 0){
                taskTimeMs = System.currentTimeMillis( ) - startTimeMs;
                System.out.printf("**** %d %% Done! ****  %d Seconds Remaining\n" , (100*timeIndex/timeSteps) , taskTimeMs*(timeSteps-timeIndex)/1000 );
            }
        }   
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

    private void saveData(){
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
    
    public static void main(String[] args) {
        long start = System.currentTimeMillis();
        Pasim pas = new Pasim();
        pas.run();
        pas.postProcessing();
        pas.saveData();
        System.out.println("Done!!");
        long totalTime = System.currentTimeMillis() - start;
        System.out.printf("The time taken : %d ms\n",totalTime);
    }
}