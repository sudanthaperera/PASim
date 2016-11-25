package pasim;

public class Port {
    private VoltageProb probV;
    private CurrentProb probI;
    private final Complex impedance;
    private final boolean isSourcePort;
    private Complex[][] S;
    private Complex[] a,b;
    private double[] frequencies;
    
    public Complex[][] getSparam(){
        return S;
    }
    
    public Complex[][] getSparam(double[] frequencies, int portCount){
        Complex[][] S = new Complex[2][frequencies.length];
        
        for(int portIndex = 0; portIndex<portCount;portIndex++){
            int freqIndexCopy = 0;
            int freqIndexOriginal = 0;
            while(freqIndexOriginal < this.frequencies.length || freqIndexCopy < frequencies.length){
                if (frequencies[freqIndexCopy] == this.frequencies[freqIndexOriginal]){
                    S[portIndex][freqIndexCopy] = this.S[portIndex][freqIndexOriginal];
                    freqIndexCopy++;
                    freqIndexOriginal++;
                }
                else{
                    freqIndexOriginal++;
                }
            }
        }
        return S;
    }
    
    public double[][] getSparamAbs(double[] frequencies, int portCount){
        double S[][] = new double[2][frequencies.length];
        
        for(int portIndex = 0; portIndex<portCount;portIndex++){
            int freqIndexCopy = 0;
            int freqIndexOriginal = 0;
            while(freqIndexOriginal < this.frequencies.length){
                if(freqIndexCopy < frequencies.length){
                    if (frequencies[freqIndexCopy] == this.frequencies[freqIndexOriginal]){
                        S[portIndex][freqIndexCopy] = this.S[portIndex][freqIndexOriginal].abs();
                        freqIndexCopy++;
                        freqIndexOriginal++;
                    }
                    else{
                        freqIndexOriginal++;
                    }
                }
                else{
                    freqIndexOriginal++;
                }
            }
        }
        return S;
    }
    
    public double[][] getSparamAngle(double[] frequencies, int portCount){
        double S[][] = new double[2][frequencies.length];
        
        for(int portIndex = 0; portIndex<portCount;portIndex++){
            int freqIndexCopy = 0;
            int freqIndexOriginal = 0;
            while(freqIndexOriginal < this.frequencies.length || freqIndexCopy < frequencies.length){
                if (frequencies[freqIndexCopy] == this.frequencies[freqIndexOriginal]){
                    S[portIndex][freqIndexCopy] = this.S[portIndex][freqIndexOriginal].angle();
                    freqIndexCopy++;
                    freqIndexOriginal++;
                }
                else{
                    freqIndexOriginal++;
                }
            }
        }
        return S;
    }
    
    public void setVProb(VoltageProb probV){
        this.probV = probV;
    }
    
    public boolean isSourcePort(){
        return isSourcePort;
    }
    
    public void setIProb(CurrentProb probI){
        this.probI = probI;
    }
        
    public Port(double impedance,boolean isSourcePort){
        this.impedance = new Complex(impedance,0);
        this.isSourcePort = isSourcePort;
    }
    
    public Port(Complex impedance,boolean isSourcePort){
        this.impedance = impedance;
        this.isSourcePort = isSourcePort;
    }
    
    public void saveSParam(String portName){
        Common.save2DComplexArray(S, "SParam".concat(portName));
    }
    
    public void calAB(double[] frequencies){
        this.frequencies = frequencies;
        a = Common.genComplex1DArray(frequencies.length, new Complex());
        b = Common.genComplex1DArray(frequencies.length, new Complex());
        Complex tempV = new Complex();
        Complex tempI = new Complex();
    
        for (int ind = 0; ind < frequencies.length; ind++){
            tempV.set(probV.getFreqDomainValue(ind));
            tempI.set(probI.getFreqDomainValue(ind));
            a[ind].set(((tempV.plus((impedance.times(tempI)))).divide(Math.sqrt(impedance.real()))).times(0.5));
            b[ind].set(((tempV.minus((impedance.conj()).times(tempI))).divide(Math.sqrt(impedance.real()))).times(0.5));   
        }
    }
    /*
    public double[] getFrequency(){
        return this.frequencies;
    }
    
    public Complex getA(int index){
        return a[index];
    }
    
    public Complex getB(int index){
        return b[index];
    }
    */
    public void calSparam(Port[] ports){
        if (this.isSourcePort()) {
            S = Common.genComplex2DArray(ports.length,frequencies.length, new Complex());
            for (int portIndex=0; portIndex < ports.length; portIndex++){
                for (int ind = 0; ind < frequencies.length; ind++){
                    if (a[ind].isEqual(new Complex())){
                        S[portIndex][ind].set(new Complex(1,0));
                        System.out.println("WARNING!!   S parramater is not Calculated Correctly");
                    }
                    else{
                        S[portIndex][ind].set(ports[portIndex].b[ind].divide(this.a[ind]));
                    }
		}
            }
	}
    }
}
