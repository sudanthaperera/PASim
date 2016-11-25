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
public class FrequencyDomain {
    private double[] frequencies;
    private double start;
    private double step;
    private double end;
    private int freqCount;
    
    public FrequencyDomain(double start,double end,int freqCount){
        this.start = start;
        this.end = end;
        this.freqCount = freqCount;
        this.step = (end-start)/(freqCount-1);
        this.frequencies = Common.genDouble1DArray(freqCount, 0.0);
        for(int i=0;i<freqCount;i++){
            this.frequencies[i] = this.start + i*this.step; 
        }
    }
    
    public void setMonotone(double frequency){
        this.start = frequency;
        this.end = frequency;
        this.freqCount = 1;
        this.step = 0;
        this.frequencies = new double[1];
    }
    public void resetFrequency(double start,double end,int freqCount){
        this.start = start;
        this.end = end;
        this.freqCount = freqCount;
        this.step = (end-start)/(freqCount-1);
        this.frequencies = new double[freqCount];
    }
    
    public double[] getFreqArray(){
        return this.frequencies;
    }
    
    public int getFreqCount(){
        return this.freqCount;
    }
}

