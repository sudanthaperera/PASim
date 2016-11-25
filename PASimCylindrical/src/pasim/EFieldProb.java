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
public class EFieldProb extends Prob {
    private int is;
    private int js;
    private int ks;
    private double[] timeDomainSamplesR, timeDomainSamplesA,timeDomainSamplesZ;
    private int component;
    private double[] time;
    private Complex[] freqDomainSamples;
    private double[] frequency;
    
    public EFieldProb(Position probPos,ProblemSpace ps, Cell c, int timeSteps){
        super(probPos);
        is = (int)Math.round((r - ps.getRmin())/c.getDeltaR());
	js = (int)Math.round((a - ps.getAmin())/c.getDeltaA());
	ks = (int)Math.round((z - ps.getZmin())/c.getDeltaZ());
	timeDomainSamplesR = Common.genDouble1DArray(timeSteps,0.0);
        timeDomainSamplesA = Common.genDouble1DArray(timeSteps,0.0);
        timeDomainSamplesZ = Common.genDouble1DArray(timeSteps,0.0);
	time = Common.genDouble1DArray(timeSteps,0.0);
	
	for (int sample_time_ind = 0; sample_time_ind < timeSteps; sample_time_ind++){
            time[sample_time_ind] = (sample_time_ind + 1)*c.getDeltaT();
	}
    }
    
    public void CaptureEField(int timeIndex){
        timeDomainSamplesR[timeIndex]=Er[is][js][ks];
        timeDomainSamplesA[timeIndex]=Ea[is][js][ks];
        timeDomainSamplesZ[timeIndex]=Ez[is][js][ks];
    }
    
    public void saveEfiled(){
        Common.save1DArray(timeDomainSamplesR,"ER");
        Common.save1DArray(timeDomainSamplesA,"EY");
        Common.save1DArray(timeDomainSamplesZ,"EZ");
    }
}
