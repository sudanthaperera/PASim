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
            case Constants.RP: this.Rmin = d; break;
            case Constants.RN: this.Rmax = d; break;
            case Constants.AP: this.Amin = d; break;
            case Constants.AN: this.Amax = d; break;
            case Constants.ZP: this.Zmin = d; break;
            case Constants.ZN: this.Zmax = d; break;
        }
    }
    
    public void initVoltageProb(ProblemSpace ps, Cell c, int timeSteps){
        is = (int)Math.round((Rmin - ps.getRmin())/c.getDeltaR());
	js = (int)Math.round((Amin - ps.getAmin())/c.getDeltaA());
	ks = (int)Math.round((Zmin - ps.getZmin())/c.getDeltaZ());
	ie = (int)Math.round((Rmax - ps.getRmin())/c.getDeltaR());
	je = (int)Math.round((Amax - ps.getAmin())/c.getDeltaA());
	ke = (int)Math.round((Zmax - ps.getZmin())/c.getDeltaZ());
	timeDomainValue = Common.genDouble1DArray(timeSteps,0.0);
	time = Common.genDouble1DArray(timeSteps,0.0);
	
	if (direction == Constants.RN || direction == Constants.RP){
            Csvf = -c.getDeltaR()/((je - js + 1)*(ke - ks + 1));
	}
	
        if (direction == Constants.AN || direction == Constants.AP){
            Csvf = -c.getDeltaA()/((ke - ks + 1)*(ie - is + 1));
	}
	
        if (direction == Constants.ZN || direction == Constants.ZP){
            Csvf = -c.getDeltaZ()/((ie - is + 1)*(je - js + 1));
	}
	
        if (direction == Constants.RN || direction == Constants.AN || direction == Constants.ZN ){
            Csvf =  -1 * Csvf;
	}
	
	for (int i = 0; i < timeSteps; i++){
            time[i] = (i + 1)*c.getDeltaT();
	}
    }
    
    public void CaptureVoltage(EField E, int timeIndex){
	double sample = 0.0;
	
	if (direction == Constants.RN || direction == Constants.RP){
	    for (int i=is;i<=ie;i++){
		for (int j=js;j<=je;j++){
		    for (int k=ks;k<ke;k++){
			sample = sample + Csvf*E.getEr(i,j,k);
		    }
		}
	    }
	}
	if (direction == Constants.AN || direction == Constants.AP){
	    for (int i=is;i<=ie;i++){
		for (int j=js;j<=je;j++){
		    for (int k=ks;k<ke;k++){
			sample = sample + Csvf*E.getEa(i,j,k);
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
