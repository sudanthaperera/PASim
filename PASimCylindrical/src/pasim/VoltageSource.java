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
public class VoltageSource extends Source {
    private double resistance;
    private double magnitude;
    private int waveformType;
    private int waveformIndex;
    private double resistancePerComponent = 1e-20;
    private double[] voltagePerEfield;
    private double[] waveform;
    private double[][][] Cezs;
    private double[][][] Ceas;
    private double[][][] Cers;
    private Complex[] freqDomainValue;
    private double[] frequencies;
    
    public static int voltageSourceCount = 0;
    
    public VoltageSource(double Rmin, double Amin, double Zmin, double Rmax, double Amax, double Zmax, Cell c){
        super(Rmax,Amax,Zmax,Rmin,Amin,Zmin,c);
        this.resistance = 50.0;
        this.magnitude = 1.0;
    }
    
    public void updateVoltageSourceEfiled(int timeIndex){	
	if (direction == Constants.RN || direction == Constants.RP){
	    for (int i=is; i < ie; i++){
		for (int j=js; j < je; j++){
		    for (int k=ks; k < ke; k++){
			Er[i][j][k] = Er[i][j][k] + Cers[i][j][k]*voltagePerEfield[timeIndex];
		    }
		}
	    }
	}
	if (direction == Constants.AN || direction == Constants.AP){
	    for (int i=is; i < ie; i++){
		for (int j=js; j < je; j++){
		    for (int k=ks; k < ke; k++){
			Ea[i][j][k] = Ea[i][j][k] + Ceas[i][j][k]*voltagePerEfield[timeIndex];
		    }
		}
	    }
	}
	if (direction == Constants.ZN || direction == Constants.ZP){
	    for (int i=is; i <= ie; i++){
		for (int j=js; j <= je; j++){
		    for (int k=ks; k < ke; k++){
			Ez[i][j][k] = Ez[i][j][k] + Cezs[i-is][j-js][k-ks]*voltagePerEfield[timeIndex];
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
    
    public double getRmax(){
        return this.Rmax;
    }
    
    public double getRmin(){
        return this.Rmin;
    }
    
    public double getAmax(){
        return this.Amax;
    }
    
    public double getAmin(){
        return this.Amin;
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
    
    public void updatingCoefficients(ProblemSpace ps){
        MaterialGrid mg = ps.getMaterialGrid();
	double temp; 
	double rpc;
	rpc = resistancePerComponent;
	
	if (direction==Constants.RP || direction==Constants.RN){
            temp = (c.getDeltaT()*c.getDeltaR())/(rpc*c.getDeltaA()*c.getDeltaZ());
            Cers = Common.genDouble3DArray(ie-is, je+1-js, ke+1-ks, 0.0);
            for(int i = is;i<=ie-1;i++){
		for(int j = js;j<=je;j++){
                    for(int k = ks;k<=ke;k++){
			Cere[i][j][k] = (2*Constants.EPS0*mg.getEpsRR(i, j, k) - c.getDeltaT()*mg.getSigmaER(i, j, k) - temp)/(2*Constants.EPS0*mg.getEpsRR(i, j, k) + c.getDeltaT()*mg.getSigmaER(i, j, k) + temp);
			Cerhz[i][j][k] = (2*c.getDeltaT()/c.getDeltaA())/(2*mg.getEpsRR(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaER(i, j, k) + temp);
			Cerha[i][j][k] = -(2*c.getDeltaT()/c.getDeltaZ())/(2*mg.getEpsRR(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaER(i, j, k) + temp);
			Cers[i-is][j-js][k-ks] = -(2*c.getDeltaT()/(rpc*c.getDeltaZ()*c.getDeltaZ()))/(2*mg.getEpsRR(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaER(i, j, k) + temp);
                    }
                }
            }
        }
	
	if (direction==Constants.AP || direction==Constants.AN){
            temp = (c.getDeltaT()*c.getDeltaA())/(rpc*c.getDeltaZ()*c.getDeltaR());
            Ceas = Common.genDouble3DArray(ie+1-is, je-js, ke+1-ks, 0.0);
            for(int i = is;i<=ie;i++){
		for(int j = js;j<=je-1;j++){
                    for(int k = ks;k<=ke;k++){
			Ceae[i][j][k] = (2*Constants.EPS0*mg.getEpsRA(i, j, k) - c.getDeltaT()*mg.getSigmaEA(i, j, k) - temp)/(2*Constants.EPS0*mg.getEpsRA(i, j, k) + c.getDeltaT()*mg.getSigmaEA(i, j, k) + temp);
			Ceahr[i][j][k] = (2*c.getDeltaT()/c.getDeltaZ())/(2*mg.getEpsRA(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEA(i, j, k) + temp);
			Ceahz[i][j][k] = -(2*c.getDeltaT()/c.getDeltaR())/(2*mg.getEpsRA(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEA(i, j, k) + temp);
			Ceas[i-is][j-js][k-ks] = -(2*c.getDeltaT()/(rpc*c.getDeltaZ()*c.getDeltaR()))/(2*mg.getEpsRA(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEA(i, j, k) + temp);
                    }
		}
            }
	}
	
	if (direction==Constants.ZP || direction==Constants.ZN){
            temp = (c.getDeltaT()*c.getDeltaZ())/(rpc*c.getDeltaR()*c.getDeltaA());
            Cezs = Common.genDouble3DArray(ie+1-is, je+1-js, ke-ks, 0.0);
            for(int i = is;i<=ie;i++){
		for(int j = js;j<=je;j++){
                    for(int k = ks;k<=ke-1;k++){
			Ceze[i][j][k] = (2*Constants.EPS0*mg.getEpsRZ(i, j, k) - c.getDeltaT()*mg.getSigmaEZ(i, j, k) - temp)/(2*Constants.EPS0*mg.getEpsRZ(i, j, k) + c.getDeltaT()*mg.getSigmaEZ(i, j, k) + temp);
			Cezha[i][j][k] = (2*c.getDeltaT()/c.getDeltaR())/ (2*mg.getEpsRZ(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEZ(i, j, k) + temp);
			Cezhr[i][j][k] = -(2*c.getDeltaT()/c.getDeltaA())/ (2*mg.getEpsRZ(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEZ(i, j, k) + temp);      
			Cezs[i-is][j-js][k-ks] = -(2*c.getDeltaT()/(rpc*c.getDeltaR()*c.getDeltaA()))/(2*mg.getEpsRZ(i, j, k)*Constants.EPS0 + c.getDeltaT()*mg.getSigmaEZ(i, j, k) + temp);
                    }
                }
            }
        }
    }
}
