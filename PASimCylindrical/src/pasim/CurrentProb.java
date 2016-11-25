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
public class CurrentProb extends Prob{
    private int direction;
    private int is;
    private int js;
    private int ks;
    private int ie;
    private int je;
    private int ke;
    private double[] timeDomainValue;
    private double[] time;
    private Complex[] frequencyDomainValue;
    private double[] frequencies;
    
    public CurrentProb(Position probPos,int direction){
        super(probPos);
        this.direction = direction;
    }
    
    public void initCurrentProb(ProblemSpace ps, Cell c, int timeSteps){
        this.c = c;
        is = (int)Math.round((Rmin - ps.getRmin())/c.getDeltaR());
	js = (int)Math.round((Amin - ps.getAmin())/c.getDeltaA());
	ks = (int)Math.round((Zmin - ps.getZmin())/c.getDeltaZ());
	ie = (int)Math.round((Rmax - ps.getRmin())/c.getDeltaR());
	je = (int)Math.round((Amax - ps.getAmin())/c.getDeltaA());
	ke = (int)Math.round((Zmax - ps.getZmin())/c.getDeltaZ());
	timeDomainValue = Common.genDouble1DArray(timeSteps,0.0);
	time = Common.genDouble1DArray(timeSteps,0.0);
	
        for (int i = 0; i < timeSteps; i++){
            time[i] = (i + 1)*c.getDeltaT();
	}
    }
    
    public void CaptureCurrent(HField H, int timeIndex){
    	double sample = 0.0;
	double sample1 = 0.0;
	double sample2 = 0.0;
	double sample3 = 0.0;
	double sample4 = 0.0;
	
	if (direction == Constants.RN || direction == Constants.RP){
	    for (int j=js;j<=je;j++){
		sample1 = sample1 + c.getDeltaA()*H.getHa(ie-1,j,ks-1);
		sample3 = sample3 + c.getDeltaA()*H.getHa(ie-1,j,ke);
	    }
	    for (int k=ks; k<=ke; k++){
		sample2 = sample2 + c.getDeltaZ()*H.getHz(ie-1,je,k);
		sample4 = sample4 + c.getDeltaZ()*H.getHz(ie-1,js-1,k);
	    }
	}
	
	else if (direction == Constants.AN || direction == Constants.AP){
	    for (int i=is;i<=ie;i++){
		sample1 = sample1 + c.getDeltaR()*H.getHr(i,je-1,ke);
		sample3 = sample3 + c.getDeltaR()*H.getHr(i,je-1,ks-1);
	    }
	    for (int k=ks; k<=ke; k++){
		sample2 = sample2 + c.getDeltaZ()*H.getHz(is-1,je-1,k);
		sample4 = sample4 + c.getDeltaZ()*H.getHz(ie,je-1,k);
	    }
	}
	
	else if (direction == Constants.ZN || direction == Constants.ZP){
	    for (int i=is;i<=ie;i++){
		sample1 = sample1 + c.getDeltaR()*H.getHr(i,js-1,ke-1);
		sample3 = sample3 + c.getDeltaR()*H.getHr(i,je,ke-1);
	    }
	    for (int j=js; j<=je; j++){
		sample2 = sample2 + c.getDeltaA()*H.getHa(ie,j,ke-1);
		sample4 = sample4 + c.getDeltaA()*H.getHa(is-1,j,ke-1);
	    }
	}
	
	sample = sample1 + sample2 - sample3 - sample4;
	
	if (direction == Constants.RN || direction == Constants.AN || direction == Constants.ZN){
	    sample = -1 * sample;
	}
	this.timeDomainValue[timeIndex] = sample;
    }
    
    public void saveTimeDomainValue(){
        Common.save1DArray(timeDomainValue, "IProbTimeDomainValue");
    }
    
    public void saveFreqDomainValue(){
        Common.save1DComplexArray(frequencyDomainValue, "IProbFrequencyDomainValue");
    }
    
    public void time2freq(double[] freqArray){
        frequencies = freqArray;
        frequencyDomainValue = Common.timeDomain2frequencyDomain(timeDomainValue, time, frequencies, -this.c.getDeltaT()/2);
    }
    
    public void saveTime(){
        Common.save1DArray(time, "IProbTime");
    }
    
    public void saveFreq(){
        Common.save1DArray(frequencies, "IProbFrequencies");
    }
    
    public Complex getFreqDomainValue(int index){
        return this.frequencyDomainValue[index];
    }    
}
