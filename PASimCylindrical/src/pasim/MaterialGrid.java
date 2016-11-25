/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pasim;

/**
 *
 * @author brah3093
 */
public class MaterialGrid extends EMobject {
    private Material[] m;
    private int[][][] index;
    private double[][][] epsrR,epsrA,epsrZ,murR,murA,murZ,sigmaeR,sigmaeA,sigmaeZ,sigmamR,sigmamA,sigmamZ;
    
    public MaterialGrid(int nr, int na, int nz, Cell c , Material[] m, int[][][] index){
        super(nr,na,nz,c);
        this.index = index;
        this.m = m;
        epsrR = new double[nr][na+1][nz+1];
        epsrA = new double[nr+1][na][nz+1];
        epsrZ = new double[nr+1][na+1][nz];
        murR = new double[nr+1][na][nz];
        murA = new double[nr][na+1][nz];
        murZ = new double[nr][na][nz+1];
        sigmaeR = new double[nr][na+1][nz+1];
        sigmaeA = new double[nr+1][na][nz+1];
        sigmaeZ = new double[nr+1][na+1][nz];
        sigmamR = new double[nr+1][na][nz];
        sigmamA = new double[nr][na+1][nz];
        sigmamZ = new double[nr][na][nz+1];
    }
    
    public void saveAllArrays(){
        Common.save3DArray(epsrR,"epsrR");
        Common.save3DArray(epsrA,"epsrA");
        Common.save3DArray(epsrZ,"epsrZ");
        
        Common.save3DArray(murR,"murR");
        Common.save3DArray(murA,"murA");
        Common.save3DArray(murZ,"murZ");

        Common.save3DArray(sigmaeR,"sigmaeR");
        Common.save3DArray(sigmaeA,"sigmaeA");
        Common.save3DArray(sigmaeZ,"sigmaeZ");
        
        Common.save3DArray(sigmamR,"sigmamR");
        Common.save3DArray(sigmamA,"sigmamA");
        Common.save3DArray(sigmamZ,"sigmamZ");
        
        Common.save3DArray(index,"index");
    }
    
    public void averageAll(){
        this.averageEpsR();
        this.averageEpsA();
        this.averageEpsZ();
        
        this.averageMuR();
        this.averageMuA();
        this.averageMuZ();
        
        this.averageSigmaER();
        this.averageSigmaEA();
        this.averageSigmaEZ();
        
        this.averageSigmaMR();
        this.averageSigmaMA();
        this.averageSigmaMZ();        
    }
    
    public void setAll(){
        this.setAllEpsRR(1.0);
        this.setAllEpsRA(1.0);
        this.setAllEpsRZ(1.0);
        
        this.setAllMuRR(1.0);
        this.setAllMuRA(1.0);
        this.setAllMuRZ(1.0);
        
        this.setAllSigmaER(0.0);
        this.setAllSigmaEA(0.0);
        this.setAllSigmaEZ(0.0);
        
        this.setAllSigmaMR(0.0);
        this.setAllSigmaMA(0.0);
        this.setAllSigmaMZ(0.0);
    }
    
    public void setAllEpsRR(double val){
        for(int i=0; i<nr; i++){
            for(int j=0; j<na+1; j++){
                for(int k=0; k<nz+1; k++){         
                    this.epsrR[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllEpsRA(double val){
        for(int i=0; i<nr+1; i++){
            for(int j=0; j<na; j++){
                for(int k=0; k<nz+1; k++){  
                    this.epsrA[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllEpsRZ(double val){
        for(int i=0; i<nr+1; i++){
            for(int j=0; j<na+1; j++){
                for(int k=0; k<nz; k++){  
                    this.epsrZ[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllMuRR(double val){
        for(int i=0; i<nr+1; i++){
            for(int j=0; j<na; j++){
                for(int k=0; k<nz; k++){  
                    this.murR[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllMuRA(double val){
        for(int i=0; i<nr; i++){
            for(int j=0; j<na+1; j++){
                for(int k=0; k<nz; k++){
                    this.murA[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllMuRZ(double val){
        for(int i=0; i<nr; i++){
            for(int j=0; j<na; j++){
                for(int k=0; k<nz+1; k++){  
                    this.murZ[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllSigmaER(double val){
        for(int i=0; i<nr; i++){
            for(int j=0; j<na+1; j++){
                for(int k=0; k<nz+1; k++){  
                    this.sigmaeR[i][j][k] = val;
                }
            }
        }
    }

    public void setAllSigmaEA(double val){
        for(int i=0; i<nr+1; i++){
            for(int j=0; j<na; j++){
                for(int k=0; k<nz+1; k++){  
                    this.sigmaeA[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllSigmaEZ(double val){
        for(int i=0; i<nr+1; i++){
            for(int j=0; j<na+1; j++){
                for(int k=0; k<nz; k++){  
                    this.sigmaeZ[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllSigmaMR(double val){
        for(int i=0; i<nr+1; i++){
            for(int j=0; j<na; j++){
                for(int k=0; k<nz; k++){  
                    this.sigmamR[i][j][k] = val;
                }
            }
        }
    }

    public void setAllSigmaMA(double val){
        for(int i=0; i<nr; i++){
            for(int j=0; j<na+1; j++){
                for(int k=0; k<nz; k++){  
                    this.sigmamA[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllSigmaMZ(double val){
        for(int i=0; i<nr; i++){
            for(int j=0; j<na; j++){
                for(int k=0; k<nz+1; k++){  
                    this.sigmamZ[i][j][k] = val;
                }
            }
        }
    }
    
    public double getEpsRR(int r, int a, int z){
        return this.epsrR[r][a][z];
    }
    
    public double getEpsRA(int r, int a, int z){
        return this.epsrA[r][a][z];
    }
    
    public double getEpsRZ(int r, int a, int z){
        return this.epsrZ[r][a][z];
    }
    
    public double getMuRR(int r, int a, int z){
        return this.murR[r][a][z];
    }
    
    public double getMuRA(int r, int a, int z){
        return this.murA[r][a][z];
    }
    
    public double getMuRZ(int r, int a, int z){
        return this.murZ[r][a][z];
    }
    
    public double getSigmaER(int r, int a, int z){
        return this.sigmaeR[r][a][z];
    }

    public double getSigmaEA(int r, int a, int z){
        return this.sigmaeA[r][a][z];
    }
    
    public double getSigmaEZ(int r, int a, int z){
        return this.sigmaeZ[r][a][z];
    }
    
    public double getSigmaMR(int r, int a, int z){
        return this.sigmamR[r][a][z];
    }

    public double getSigmaMA(int r, int a, int z){
        return this.sigmamA[r][a][z];
    }
    
    public double getSigmaMZ(int r, int a, int z){
        return this.sigmamZ[r][a][z];
    }
    
    public void setEpsRR(int r, int a, int z, double val){
        this.epsrR[r][a][z] = val;
    }
    
    public void setEpsRA(int r, int a, int z, double val){
        this.epsrA[r][a][z] = val;
    }
    
    public void setEpsRZ(int r, int a, int z, double val){
        this.epsrZ[r][a][z] = val;
    }
    
    public void setMuRR(int r, int a, int z, double val){
        this.murR[r][a][z] = val;
    }
    
    public void setMuRA(int r, int a, int z, double val){
        this.murA[r][a][z] = val;
    }
    
    public void setMuRZ(int r, int a, int z, double val){
        this.murZ[r][a][z] = val;
    }
    
    public void setSigmaER(int r, int a, int z, double val){
        this.sigmaeR[r][a][z] = val;
    }

    public void setSigmaEA(int r, int a, int z, double val){
        this.sigmaeA[r][a][z] = val;
    }
    
    public void setSigmaEZ(int r, int a, int z, double val){
        this.sigmaeZ[r][a][z] = val;
    }
    
    public void setSigmaMR(int r, int a, int z, double val){
        this.sigmamR[r][a][z] = val;
    }

    public void setSigmaMA(int r, int a, int z, double val){
        this.sigmamA[r][a][z] = val;
    }
    
    public void setSigmaMZ(int r, int a, int z, double val){
        this.sigmamZ[r][a][z] = val;
    }

    public void averageMuR(){
        int rStart = 1, rEnd = this.nr;
        int aStart = 0, aEnd = this.na;
        int zStart = 0, zEnd = this.nz;
	
	for(int i=rStart;i<rEnd;i++){
            for(int j=aStart;j<aEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.murR[i][j][k] = (2.0*(this.m[this.index[i][j][k]].muR())*(this.m[this.index[i-1][j][k]].muR()))/(this.m[this.index[i][j][k]].muR() + this.m[this.index[i-1][j][k]].muR());
		}
            }
	}
    }

    public void averageMuA(){
        int rStart = 0, rEnd = this.nr;
        int aStart = 1, aEnd = this.na;
        int zStart = 0, zEnd = this.nz;
        
	for(int i=rStart;i<rEnd;i++){
            for(int j=aStart;j<aEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.murA[i][j][k] = (2.0*(this.m[this.index[i][j][k]].muR())*(this.m[this.index[i][j-1][k]].muR()))/(this.m[this.index[i][j][k]].muR() + this.m[this.index[i][j-1][k]].muR());
		}
            }
	}
    }

    public void averageMuZ(){
        int rStart = 0, rEnd = this.nr;
        int aStart = 0, aEnd = this.na;
        int zStart = 1, zEnd = this.nz;
	
	for(int i=rStart;i<rEnd;i++){
            for(int j=aStart;j<aEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.murZ[i][j][k] = (2.0*(this.m[this.index[i][j][k]].muR())*(this.m[this.index[i][j][k-1]].muR()))/(this.m[this.index[i][j][k]].muR() + this.m[this.index[i][j][k-1]].muR());
		}
            }
	}
    }

    public void averageSigmaMR(){
        int rStart = 1, rEnd = this.nr;
        int aStart = 0, aEnd = this.na;
        int zStart = 0, zEnd = this.nz;
	
	for(int i=rStart;i<rEnd;i++){
            for(int j=aStart;j<aEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.sigmamR[i][j][k] = (2.0*(this.m[this.index[i][j][k]].sigmaM())*(this.m[this.index[i-1][j][k]].sigmaM()))/(this.m[this.index[i][j][k]].sigmaM() + this.m[this.index[i-1][j][k]].sigmaM());
		}
            }
	}
    }

    public void averageSigmaMA(){
        int rStart = 0, rEnd = this.nr;
        int aStart = 1, aEnd = this.na;
        int zStart = 0, zEnd = this.nz;
	
	for(int i=rStart;i<rEnd;i++){
            for(int j=aStart;j<aEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.sigmamA[i][j][k] = (2.0*(this.m[this.index[i][j][k]].sigmaM())*(this.m[this.index[i][j-1][k]].sigmaM()))/(this.m[this.index[i][j][k]].sigmaM() + this.m[this.index[i][j-1][k]].sigmaM());
		}
            }
	}
    }

    public void averageSigmaMZ(){
        int rStart = 0, rEnd = this.nr;
        int aStart = 0, aEnd = this.na;
        int zStart = 1, zEnd = this.nz;
	
	for(int i=rStart;i<rEnd;i++){
            for(int j=aStart;j<aEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.sigmamZ[i][j][k] = (2.0*(this.m[this.index[i][j][k]].sigmaM())*(this.m[this.index[i][j][k-1]].sigmaM()))/(this.m[this.index[i][j][k]].sigmaM() + this.m[this.index[i][j][k-1]].sigmaM());
		}
            }
	}
    }

    public void averageEpsR(){
        int rStart = 0, rEnd = this.nr;
        int aStart = 1, aEnd = this.na;
        int zStart = 1, zEnd = this.nz;
	
	for(int i=rStart;i<rEnd;i++){
            for(int j=aStart;j<aEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.epsrR[i][j][k] = 0.25*(this.m[this.index[i][j][k]].epsR() + this.m[this.index[i][j-1][k]].epsR() + this.m[this.index[i][j][k-1]].epsR() + this.m[this.index[i][j-1][k-1]].epsR());
		}
            }
	}
    }

    public void averageEpsA(){
        int rStart = 1, rEnd = this.nr;
        int aStart = 0, aEnd = this.na;
        int zStart = 1, zEnd = this.nz;
	
	for(int i=rStart;i<rEnd;i++){
            for(int j=aStart;j<aEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.epsrA[i][j][k] = 0.25*(this.m[this.index[i][j][k]].epsR() + this.m[this.index[i-1][j][k]].epsR()+ this.m[this.index[i][j][k-1]].epsR() + this.m[this.index[i-1][j][k-1]].epsR());
		}
            }
	}
    }
 
    public void averageEpsZ(){
        int rStart = 1, rEnd = this.nr;
        int aStart = 1, aEnd = this.na;
        int zStart = 0, zEnd = this.nz;
	
	for(int i=rStart;i<rEnd;i++){
            for(int j=aStart;j<aEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.epsrZ[i][j][k] = 0.25*(this.m[this.index[i][j][k]].epsR() + this.m[this.index[i-1][j][k]].epsR() + this.m[this.index[i][j-1][k]].epsR() + this.m[this.index[i-1][j-1][k]].epsR());
		}
            }
	}
    }
 
    public void averageSigmaER(){
        int rStart = 0, rEnd = this.nr;
        int aStart = 1, aEnd = this.na;
        int zStart = 1, zEnd = this.nz;
	
	for(int i=rStart;i<rEnd;i++){
            for(int j=aStart;j<aEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.sigmaeR[i][j][k] = 0.25*(this.m[this.index[i][j][k]].sigmaE() + this.m[this.index[i][j-1][k]].sigmaE() + this.m[this.index[i][j][k-1]].sigmaE() + this.m[this.index[i][j-1][k-1]].sigmaE());
		}
            }
	}
    }
 
    void averageSigmaEA(){
        int rStart = 1, rEnd = this.nr;
        int aStart = 0, aEnd = this.na;
        int zStart = 1, zEnd = this.nz;
	
	for(int i=rStart;i<rEnd;i++){
            for(int j=aStart;j<aEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.sigmaeA[i][j][k] = 0.25*(this.m[this.index[i][j][k]].sigmaE() + this.m[this.index[i-1][j][k]].sigmaE() + this.m[this.index[i][j][k-1]].sigmaE() + this.m[this.index[i-1][j][k-1]].sigmaE());
		}
            }
	}
    }
 
    public void averageSigmaEZ(){
        int rStart = 1, rEnd = this.nr;
        int aStart = 1, aEnd = this.na;
        int zStart = 0, zEnd = this.nz;
	
	for(int i=rStart;i<rEnd;i++){
            for(int j=aStart;j<aEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.sigmaeZ[i][j][k] = 0.25*(this.m[this.index[i][j][k]].sigmaE() + this.m[this.index[i-1][j][k]].sigmaE() + this.m[this.index[i][j-1][k]].sigmaE() + this.m[this.index[i-1][j-1][k]].sigmaE());
		}
            }
	}
    }
}
