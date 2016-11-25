/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pasim;

import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLDouble;
import java.util.ArrayList;

/**
 *
 * @author brah3093
 */
public class Common {
    
    public static Complex complexExp(Complex number){
	double r = number.real();
	double i = number.imag();        
	return (new Complex(Math.cos(i), Math.sin(i))).times(Math.exp(r));
    }
    
    public static <T extends Number> T[][][] gen3DArray(int nx,int ny,int nz,T initValue){
        nx = Math.abs(nx);
        ny = Math.abs(ny);
        nz = Math.abs(nz);
        
        T[][][] array = (T[][][])new Number[nx][ny][nz];
        
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz; k++){ 
                    array[i][j][k] = initValue;
                }
            }
        }
        return array;
    }
    
    public static <T extends Number> T[][] gen2DArray(int nx,int ny,T initValue){
        nx = Math.abs(nx);
        ny = Math.abs(ny);
        T[][] array = (T[][])new Number[nx][ny];
        
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){ 
                array[i][j] = initValue;
            }
        }
        return array;
    }
    
    public static <T extends Number> T[] gen1DArray(int nx,T initValue){
        nx = Math.abs(nx);
        T[] array = (T[])new Number[nx];
        
        for(int i=0; i<nx; i++){ 
            array[i] = initValue;
        }
        return array;
    }
    
    public static Complex[] genComplex1DArray(int nx,Complex initValue){
        nx = Math.abs(nx);
        Complex[] array = new Complex[nx];
        
        for(int i=0; i<nx; i++){
            array[i] = new Complex(initValue);
        }
        return array;
    }
    
    public static Complex[][] genComplex2DArray(int nx,int ny,Complex initValue){
        nx = Math.abs(nx);
        ny = Math.abs(ny);

        Complex[][] array = new Complex[nx][ny];
        
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                array[i][j] = new Complex(initValue);
            }
        }
        return array;
    }
    
    public static Complex[][][][] genComplex4DArray(int nx,int ny,int nz,int nw,Complex initValue){
        nx = Math.abs(nx);
        ny = Math.abs(ny);
        nz = Math.abs(nz);
        nw = Math.abs(nw);
        Complex[][][][] array = new Complex[nx][ny][nz][nw];
        
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz; k++){ 
                    for(int l=0; l<nw; l++){
                        array[i][j][k][l] = new Complex(initValue);
                    }
                }
            }
        }
        return array;
    }
    
    public static Complex[][][] genComplex3DArray(int nx,int ny,int nz,Complex initValue){
        nx = Math.abs(nx);
        ny = Math.abs(ny);
        nz = Math.abs(nz);
        Complex[][][] array = new Complex[nx][ny][nz];
        
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz; k++){ 
                    array[i][j][k] = new Complex(initValue);
                }
            }
        }
        return array;
    }
    
    
    public static double[][][][] genDouble4DArray(int nx,int ny,int nz,int nw,double initValue){
        nx = Math.abs(nx);
        ny = Math.abs(ny);
        nz = Math.abs(nz);
        nw = Math.abs(nw);
        double[][][][] array = new double[nx][ny][nz][nw];
        
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz; k++){ 
                    for(int l=0; l<nw; l++){
                        array[i][j][k][l] = initValue;
                    }
                }
            }
        }
        return array;
    }
    
    public static double[][][] genDouble3DArray(int nx,int ny,int nz,double initValue){
        nx = Math.abs(nx);
        ny = Math.abs(ny);
        nz = Math.abs(nz);
        double[][][] array = new double[nx][ny][nz];
        
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz; k++){ 
                    array[i][j][k] = initValue;
                }
            }
        }
        return array;
    }
    
    public static double[][] genDouble2DArray(int nx,int ny,double initValue){
        nx = Math.abs(nx);
        ny = Math.abs(ny);
        double[][] array = new double[nx][ny];
        
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){ 
                array[i][j] = initValue;
            }
        }
        return array;
    }
    
    public static double[] genDouble1DArray(int nx,double initValue){
        nx = Math.abs(nx);
        double[] array = new double[nx];
        
        for(int i=0; i<nx; i++){ 
            array[i] = initValue;
        }
        return array;
    }
    
    public static int[][][] genInt3DArray(int nx,int ny,int nz,int initValue){
        nx = Math.abs(nx);
        ny = Math.abs(ny);
        nz = Math.abs(nz);
        int[][][] array = new int[nx][ny][nz];
        
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz; k++){ 
                    array[i][j][k] = initValue;
                }
            }
        }
        return array;
    }
    
    public static void save3DArray(double[][][] a,String fileName){
        int iCount = a.length;
        int jCount = a[0].length;
        int kCount = a[0][0].length;
        
        ArrayList list = new ArrayList();
        MLDouble[] Array = new MLDouble[kCount];
        double[][] tempArray = new double[iCount][jCount];
        
        for (int k=0;k<kCount;k++){
            for (int i=0;i<iCount;i++){
                for (int j=0;j<jCount;j++){
                    tempArray[i][j] = a[i][j][k];
                }
            }
            
            Array[k] = new MLDouble( "a".concat(String.valueOf(k)), tempArray);
            list.add(Array[k]);
        }
        
        try{
            new MatFileWriter( fileName.concat(".mat"), list);
        }
        catch (java.io.IOException e){
            System.out.println("****The data file was not saved****");
            System.out.println(e.getLocalizedMessage());
        }
    }

    public static void save3DArray(int[][][] a,String fileName){
        int iCount = a.length;
        int jCount = a[0].length;
        int kCount = a[0][0].length;
        
        ArrayList list = new ArrayList();
        MLDouble[] Array = new MLDouble[kCount];
        double[][] tempArray = new double[iCount][jCount];
        
        for (int k=0;k<kCount;k++){
            for (int i=0;i<iCount;i++){
                for (int j=0;j<jCount;j++){
                    tempArray[i][j] = (double)a[i][j][k];
                }
            }
            
            Array[k] = new MLDouble( "a".concat(String.valueOf(k)), tempArray);
            list.add(Array[k]);
        }
        
        try{
            new MatFileWriter( fileName.concat(".mat"), list);
        }
        catch (java.io.IOException e){
            System.out.println("****The data file was not saved****");
            System.out.println(e.getLocalizedMessage());
        }
    }
    
    public static void save2DComplexArray(Complex[][] a,String fileName){
        ArrayList list = new ArrayList();
        MLDouble realArray;
        MLDouble imagArray;
        double[][] real = Common.isolateRealFrom2DArray(a);
        double[][] imag = Common.isolateImagFrom2DArray(a);
        
            
        realArray = new MLDouble("Real", real);
        imagArray = new MLDouble("Imag", imag);
        list.add(realArray);
        list.add(imagArray);
        
        try{
            new MatFileWriter( fileName.concat(".mat"), list);
        }
        catch (java.io.IOException e){
            System.out.println("****The data file was not saved****");
            System.out.println(e.getLocalizedMessage());
        }
    }    
    
    public static void save2DArray(double[][] a,String fileName){
        ArrayList list = new ArrayList();
        MLDouble Array;
            
        Array = new MLDouble( fileName, a);
        list.add(Array);
        
        try{
            new MatFileWriter( fileName.concat(".mat"), list);
        }
        catch (java.io.IOException e){
            System.out.println("****The data file was not saved****");
            System.out.println(e.getLocalizedMessage());
        }
    }
           
    public static void save1DArray(double[] a,String fileName){
        ArrayList list = new ArrayList();
        MLDouble Array;
            
        Array = new MLDouble( fileName, a, 1);
        list.add(Array);
        
        try{
            new MatFileWriter( fileName.concat(".mat"), list);
        }
        catch (java.io.IOException e){
            System.out.println("****The data file was not saved****");
            System.out.println(e.getLocalizedMessage());
        }
    }
    
    public static void save1DComplexArray(Complex[] a,String fileName){
        ArrayList list = new ArrayList();
        MLDouble RealArray;
        MLDouble ImagArray;
            
        RealArray = new MLDouble( "Real", isolateRealFrom1DArray(a), 1);
        ImagArray = new MLDouble( "Imag", isolateImagFrom1DArray(a), 1);
        list.add(RealArray);
        list.add(ImagArray);
        
        try{
            new MatFileWriter( fileName.concat(".mat"), list);
        }
        catch (java.io.IOException e){
            System.out.println("****The data file was not saved****");
            System.out.println(e.getLocalizedMessage());
        }
    }
    
    public static double[][] isolateRealFrom2DArray(Complex[][] a){
        double[][] array = new double[a.length][a[0].length];
        for(int i=0;i<a.length;i++){
            for(int j=0;j<a[0].length;j++){
                array[i][j] = a[i][j].real();
            }
        }
        return array;
    }
    
    public static double[][] isolateImagFrom2DArray(Complex[][] a){
        double[][] array = new double[a.length][a[0].length];
        for(int i=0;i<a.length;i++){
            for(int j=0;j<a[0].length;j++){
                array[i][j] = a[i][j].imag();
            }
        }
        return array;
    }    
    
    public static double[] isolateRealFrom1DArray(Complex[] a){
        double[] array = new double[a.length];
        for(int i=0;i<a.length;i++){
            array[i] = a[i].real();
        }
        return array;
    }
    
    public static double[] isolateImagFrom1DArray(Complex[] a){
        double[] array = new double[a.length];
        for(int i=0;i<a.length;i++){
            array[i] = a[i].imag();
        }
        return array;
    }
    
    public static Complex[] timeDomain2frequencyDomain(double[] x, double[] time, double[] frequencies, double timeShift){
	Complex[] array = new Complex[frequencies.length];
        double angle;
        double dt = time[1]-time[0];

	for (int freqIndex = 0; freqIndex < frequencies.length; freqIndex++){
            array[freqIndex] = new Complex();
            for (int timeIndex = 0; timeIndex < time.length; timeIndex++){
                angle = -(2*Math.PI*frequencies[freqIndex])*(time[timeIndex] + timeShift);
		array[freqIndex].set(((new Complex(angle)).times(x[timeIndex])).plus(array[freqIndex]));
            }
            array[freqIndex].set(array[freqIndex].times(dt));
	}
	return array;
    }
}
