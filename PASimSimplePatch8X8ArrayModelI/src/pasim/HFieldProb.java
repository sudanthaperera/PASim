package pasim;

public class HFieldProb extends Prob{
    private int is;
    private int js;
    private int ks;
    private double[] timeDomainSamples;
    private int component;
    private double[] time;
    private Complex[] freqDomainSamples;
    private double[] frequency;
    
    public HFieldProb(Position probPos, ProblemSpace ps, Cell c, int timeSteps){
        super(probPos);
	is = (int)Math.round((x - ps.getXmin())/c.getDeltaX());
	js = (int)Math.round((y - ps.getYmin())/c.getDeltaY());
	ks = (int)Math.round((z - ps.getZmin())/c.getDeltaZ());
	timeDomainSamples = Common.genDouble1DArray(timeSteps,0.0);
	time = Common.genDouble1DArray(timeSteps,0.0);
	
	for (int sample_time_ind = 0; sample_time_ind < timeSteps; sample_time_ind++){
            time[sample_time_ind] = (sample_time_ind + 0.5)*c.getDeltaT();
	}        
    }
    
    public void CaptureHField(int timeIndex, HField H){
        //to be use
    }        
}
