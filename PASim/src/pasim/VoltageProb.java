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
public class VoltageProb extends Prob {
    private int direction;
    private int is;
    private int js;
    private int ks;
    private int ie;
    private int je;
    private int ke;
    private double[] timeDomainValue;
    private double Csvf;
    private double[] time;
    private Complex[] frequencyDomainValue;
    private double[] frequencies;

    
    public VoltageProb(Position probPos,int direction, double d){
        super(probPos);
        this.direction = direction;
        switch(direction){
            case Constants.XP: this.Xmin = d; break;
            case Constants.XN: this.Xmax = d; break;
            case Constants.YP: this.Ymin = d; break;
            case Constants.YN: this.Ymax = d; break;
            case Constants.ZP: this.Zmin = d; break;
            case Constants.ZN: this.Zmax = d; break;
        }
    }
    
    public void initVoltageProb(ProblemSpace ps, Cell c, int timeSteps){
        is = (int)Math.round((Xmin - ps.getXmin())/c.getDeltaX());
	js = (int)Math.round((Ymin - ps.getYmin())/c.getDeltaY());
	ks = (int)Math.round((Zmin - ps.getZmin())/c.getDeltaZ());
	ie = (int)Math.round((Xmax - ps.getXmin())/c.getDeltaX());
	je = (int)Math.round((Ymax - ps.getYmin())/c.getDeltaY());
	ke = (int)Math.round((Zmax - ps.getZmin())/c.getDeltaZ());
	timeDomainValue = Common.genDouble1DArray(timeSteps,0.0);
	time = Common.genDouble1DArray(timeSteps,0.0);
	
	if (direction == Constants.XN || direction == Constants.XP){
            Csvf = -c.getDeltaX()/((je - js + 1)*(ke - ks + 1));
	}
	
        if (direction == Constants.YN || direction == Constants.YP){
            Csvf = -c.getDeltaY()/((ke - ks + 1)*(ie - is + 1));
	}
	
        if (direction == Constants.ZN || direction == Constants.ZP){
            Csvf = -c.getDeltaZ()/((ie - is + 1)*(je - js + 1));
	}
	
        if (direction == Constants.XN || direction == Constants.YN || direction == Constants.ZN ){
            Csvf =  -1 * Csvf;
	}
	
	for (int i = 0; i < timeSteps; i++){
            time[i] = (i + 1)*c.getDeltaT();
	}
    }
    
    public void CaptureVoltage(EField E, int timeIndex){
	double sample = 0.0;
	
	if (direction == Constants.XN || direction == Constants.XP){
	    for (int i=is;i<=ie;i++){
		for (int j=js;j<=je;j++){
		    for (int k=ks;k<ke;k++){
			sample = sample + Csvf*E.getEx(i,j,k);
		    }
		}
	    }
	}
	if (direction == Constants.YN || direction == Constants.YP){
	    for (int i=is;i<=ie;i++){
		for (int j=js;j<=je;j++){
		    for (int k=ks;k<ke;k++){
			sample = sample + Csvf*E.getEy(i,j,k);
		    }
		}
	    }
	}
	if (direction == Constants.ZN || direction == Constants.ZP){
	    for (int i=is;i<=ie;i++){
		for (int j=js;j<=je;j++){
		    for (int k=ks;k<ke;k++){
			sample = sample + Csvf*E.getEz(i,j,k);
		    }
		}
	    }
	}
	
	timeDomainValue[timeIndex] = sample;
    }
    
    public void saveTimeDomainValue(){
        Common.save1DArray(timeDomainValue, "VProbTimeDomainValue");
    }
    
    public void saveFreqDomainValue(){
        Common.save1DComplexArray(frequencyDomainValue, "VProbFrequencyDomainValue");
    }
    
    public void time2freq(double[] freqArray){
        frequencies = freqArray;
        frequencyDomainValue = Common.timeDomain2frequencyDomain(timeDomainValue, time, frequencies, 0.0);
    }
    
    public void saveTime(){
        Common.save1DArray(time, "VProbTime");
    }
    
    public void saveFreq(){
        Common.save1DArray(frequencies, "VProbFrequencies");
    }
    
    public Complex getFreqDomainValue(int index){
        return this.frequencyDomainValue[index];
    }
}
