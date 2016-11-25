package pasim;

public class VoltageSource extends Source {
    private double resistance;
    private double magnitude;
    private int waveformType;
    private int waveformIndex;
    private double resistancePerComponent = 1e-20;
    private double[] voltagePerEfield;
    private double[] waveform;
    private double[][][] Cezs;
    private double[][][] Ceys;
    private double[][][] Cexs;
    private Complex[] freqDomainValue;
    private double[] frequencies;
    private double[] time;
    
    public VoltageSource(double Xmin, double Ymin, double Zmin, double Xmax, double Ymax, double Zmax, Cell c){
        super(Xmax,Ymax,Zmax,Xmin,Ymin,Zmin,c);
        this.resistance = 50.0;
        this.magnitude = 1.0;
    }
    
    public VoltageSource(VoltageSource vs){
        super(0,0,0,0,0,0,null);
        this.Xmin = vs.getXmin();
        this.Ymin = vs.getYmin();
        this.Zmin = vs.getZmin();
        this.Xmax = vs.getXmax();
        this.Ymax = vs.getYmax();
        this.Zmax = vs.getZmax();
        this.c = vs.getCell();
        this.resistance = vs.getResistance();
        this.magnitude = vs.getMagnitude();
    }
    
    public Cell getCell(){
        return this.c;
    }
    
    public void move(double x, double y, double z){
        this.Xmin = this.Xmin + x;
        this.Xmax = this.Xmax + x;
        this.Ymin = this.Ymin + y;
        this.Ymax = this.Ymax + y;
        this.Zmin = this.Zmin + z;
        this.Zmax = this.Zmax + z;
    }
    
    public void calFreqDomainValues(){
        this.freqDomainValue =  Common.timeDomain2frequencyDomain(waveform, time, this.frequencies, 0);
    }
    
    public void saveWaveform(String fileName){
        Common.save1DArray(waveform, fileName.concat("TimeDomain"));
        Common.save1DComplexArray(freqDomainValue, fileName.concat("FreqDomain"));
        Common.save1DArray(time, "Time".concat(fileName));
        Common.save1DArray(this.frequencies, "Freq".concat(fileName));
    }
    
    public void setTime(double[] time){
        this.time = time;
    }
    
    public void setFreq(double[] freq){
        this.frequencies = freq;
    }
    
    public void updateVoltageSourceEfiled(int timeIndex, EField E){	
	if (direction == Constants.XN || direction == Constants.XP){
	    for (int i=is; i < ie; i++){
		for (int j=js; j < je; j++){
		    for (int k=ks; k < ke; k++){
			E.setEX(i,j,k, E.getEX(i,j,k) + Cexs[i][j][k]*voltagePerEfield[timeIndex]);
		    }
		}
	    }
	}
	if (direction == Constants.YN || direction == Constants.YP){
	    for (int i=is; i < ie; i++){
		for (int j=js; j < je; j++){
		    for (int k=ks; k < ke; k++){
			E.setEY(i,j,k, E.getEY(i,j,k) + Ceys[i][j][k]*voltagePerEfield[timeIndex]);
		    }
		}
	    }
	}
	if (direction == Constants.ZN || direction == Constants.ZP){
	    for (int i=is; i <= ie; i++){
		for (int j=js; j <= je; j++){
		    for (int k=ks; k < ke; k++){
			E.setEZ(i,j,k, E.getEZ(i,j,k) + Cezs[i-is][j-js][k-ks]*voltagePerEfield[timeIndex]);
		    }
		}
	    }
	}
    }
    
    
    public void initVoltagePerEfield(int timeSteps){
        this.voltagePerEfield = new double[timeSteps];
        for(int i=0;i<timeSteps;i++){
            this.voltagePerEfield[i] = 0.0;
        }
    }
    
    public void initWaveform(int timeSteps){
        this.waveform = new double[timeSteps];
        for(int i=0;i<timeSteps;i++){
            this.waveform[i] = 0.0;
        }
    }
    
    public void setVoltagePerEfield(int index, double value){
        this.voltagePerEfield[index] = value;
    }
    
    public void setWaveform(int index, double value){
        this.waveform[index] = value;
    }
    
    public void setIS(int is){
        this.is = is;
    }
    
    public void setJS(int js){
        this.js = js;
    }
    
    public void setKS(int ks){
        this.ks = ks;
    }
    
    public void setIE(int ie){
        this.ie = ie;
    }
    
    public void setJE(int je){
        this.je = je;
    }
    
    public void setKE(int ke){
        this.ke = ke;
    }
    
    public void setDirection(int direction){
        this.direction = direction;
    }
    
    public int getDirection(){
        return this.direction;
    }
    
    public double getMagnitude(){
        return this.magnitude;
    }
    
    public double getXmax(){
        return this.Xmax;
    }
    
    public double getXmin(){
        return this.Xmin;
    }
    
    public double getYmax(){
        return this.Ymax;
    }
    
    public double getYmin(){
        return this.Ymin;
    }
    
    public double getZmax(){
        return this.Zmax;
    }
    
    public double getZmin(){
        return this.Zmin;
    }
    
    public double getResistancePerComponent(){
        return this.resistancePerComponent;
    }
    
    public void setResistancePerComponent(double resistance_per_component){
        this.resistancePerComponent =resistance_per_component;
    }
    
    public void setResistance(double resistance){
        this.resistance =resistance;
    }
    
    public void setMagnitude(double magnitude){
        this.magnitude = magnitude;
    }
    
    public double getResistance(){
        return this.resistance;
    }
    
    public void updatingCoefficients(ProblemSpace ps,EField E, HField H){
        MaterialGrid mg = ps.getMaterialGrid();
	double temp; 
	double rpc;
	rpc = resistancePerComponent;
	
	if (direction==Constants.XP || direction==Constants.XN){
            temp = (c.getDeltaT()*c.getDeltaX())/(rpc*c.getDeltaY()*c.getDeltaZ());
            Cexs = Common.genDouble3DArray(ie-is, je+1-js, ke+1-ks, 0.0);
            for(int i = is;i<=ie-1;i++){
		for(int j = js;j<=je;j++){
                    for(int k = ks;k<=ke;k++){
			E.setCexe(i,j,k, (2*Constants.EPS0*mg.getEpsRX(i, j, k) - c.getDeltaT()*mg.getSigmaEX(i, j, k) - temp)/(2*Constants.EPS0*mg.getEpsRX(i, j, k) + c.getDeltaT()*mg.getSigmaEX(i, j, k) + temp));
			E.setCexhz(i,j,k, (2*c.getDeltaT()/c.getDeltaY())/(2*mg.getEpsRX(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEX(i, j, k) + temp));
			E.setCexhy(i,j,k, -(2*c.getDeltaT()/c.getDeltaZ())/(2*mg.getEpsRX(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEX(i, j, k) + temp));
			Cexs[i-is][j-js][k-ks] = -(2*c.getDeltaT()/(rpc*c.getDeltaZ()*c.getDeltaZ()))/(2*mg.getEpsRX(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEX(i, j, k) + temp);
                    }
                }
            }
        }
	
	if (direction==Constants.YP || direction==Constants.YN){
            temp = (c.getDeltaT()*c.getDeltaY())/(rpc*c.getDeltaZ()*c.getDeltaX());
            Ceys = Common.genDouble3DArray(ie+1-is, je-js, ke+1-ks, 0.0);
            for(int i = is;i<=ie;i++){
		for(int j = js;j<=je-1;j++){
                    for(int k = ks;k<=ke;k++){
			E.setCeye(i,j,k, (2*Constants.EPS0*mg.getEpsRY(i, j, k) - c.getDeltaT()*mg.getSigmaEY(i, j, k) - temp)/(2*Constants.EPS0*mg.getEpsRY(i, j, k) + c.getDeltaT()*mg.getSigmaEY(i, j, k) + temp));
			E.setCeyhx(i,j,k, (2*c.getDeltaT()/c.getDeltaZ())/(2*mg.getEpsRY(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEY(i, j, k) + temp));
			E.setCeyhz(i,j,k, -(2*c.getDeltaT()/c.getDeltaX())/(2*mg.getEpsRY(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEY(i, j, k) + temp));
			Ceys[i-is][j-js][k-ks] = -(2*c.getDeltaT()/(rpc*c.getDeltaZ()*c.getDeltaX()))/(2*mg.getEpsRY(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEY(i, j, k) + temp);
                    }
		}
            }
	}
	
	if (direction==Constants.ZP || direction==Constants.ZN){
            temp = (c.getDeltaT()*c.getDeltaZ())/(rpc*c.getDeltaX()*c.getDeltaY());
            Cezs = Common.genDouble3DArray(ie+1-is, je+1-js, ke-ks, 0.0);
            for(int i = is;i<=ie;i++){
		for(int j = js;j<=je;j++){
                    for(int k = ks;k<=ke-1;k++){
			E.setCeze(i,j,k, (2*Constants.EPS0*mg.getEpsRZ(i, j, k) - c.getDeltaT()*mg.getSigmaEZ(i, j, k) - temp)/(2*Constants.EPS0*mg.getEpsRZ(i, j, k) + c.getDeltaT()*mg.getSigmaEZ(i, j, k) + temp));
			E.setCezhy(i,j,k, (2*c.getDeltaT()/c.getDeltaX())/ (2*mg.getEpsRZ(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEZ(i, j, k) + temp));
			E.setCezhx(i,j,k, -(2*c.getDeltaT()/c.getDeltaY())/ (2*mg.getEpsRZ(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEZ(i, j, k) + temp));      
			Cezs[i-is][j-js][k-ks] = -(2*c.getDeltaT()/(rpc*c.getDeltaX()*c.getDeltaY()))/(2*mg.getEpsRZ(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEZ(i, j, k) + temp);
                    }
                }
            }
        }
    }
}
