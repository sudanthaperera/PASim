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
    private double[] timeDomainSamplesX, timeDomainSamplesY,timeDomainSamplesZ;
    private int component;
    private double[] time;
    private Complex[] freqDomainSamples;
    private double[] frequency;
    
    public EFieldProb(Position probPos,ProblemSpace ps, Cell c, int timeSteps){
        super(probPos);
        is = (int)Math.round((x - ps.getXmin())/c.getDeltaX());
	js = (int)Math.round((y - ps.getYmin())/c.getDeltaY());
	ks = (int)Math.round((z - ps.getZmin())/c.getDeltaZ());
	timeDomainSamplesX = Common.genDouble1DArray(timeSteps,0.0);
        timeDomainSamplesY = Common.genDouble1DArray(timeSteps,0.0);
        timeDomainSamplesZ = Common.genDouble1DArray(timeSteps,0.0);
	time = Common.genDouble1DArray(timeSteps,0.0);
	
	for (int sample_time_ind = 0; sample_time_ind < timeSteps; sample_time_ind++){
            time[sample_time_ind] = (sample_time_ind + 1)*c.getDeltaT();
	}
    }
    
    public void CaptureEField(int timeIndex){
        timeDomainSamplesX[timeIndex]=Ex[is][js][ks];
        timeDomainSamplesY[timeIndex]=Ey[is][js][ks];
        timeDomainSamplesZ[timeIndex]=Ez[is][js][ks];
    }
    
    public void saveEfiled(){
        Common.save1DArray(timeDomainSamplesX,"EX");
        Common.save1DArray(timeDomainSamplesY,"EY");
        Common.save1DArray(timeDomainSamplesZ,"EZ");
    }
}
