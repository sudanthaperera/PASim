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
    private double[][][] epsrX,epsrY,epsrZ,murX,murY,murZ,sigmaeX,sigmaeY,sigmaeZ,sigmamX,sigmamY,sigmamZ;
    
    public MaterialGrid(int nx, int ny, int nz, Cell c , Material[] m, int[][][] index){
        super(nx,ny,nz,c);
        this.index = index;
        this.m = m;
        epsrX = new double[nx][ny+1][nz+1];
        epsrY = new double[nx+1][ny][nz+1];
        epsrZ = new double[nx+1][ny+1][nz];
        murX = new double[nx+1][ny][nz];
        murY = new double[nx][ny+1][nz];
        murZ = new double[nx][ny][nz+1];
        sigmaeX = new double[nx][ny+1][nz+1];
        sigmaeY = new double[nx+1][ny][nz+1];
        sigmaeZ = new double[nx+1][ny+1][nz];
        sigmamX = new double[nx+1][ny][nz];
        sigmamY = new double[nx][ny+1][nz];
        sigmamZ = new double[nx][ny][nz+1];
    }
    
    public void saveAllArrays(){
        Common.save3DArray(epsrX,"epsrX");
        Common.save3DArray(epsrY,"epsrY");
        Common.save3DArray(epsrZ,"epsrZ");
        
        Common.save3DArray(murX,"murX");
        Common.save3DArray(murY,"murY");
        Common.save3DArray(murZ,"murZ");

        Common.save3DArray(sigmaeX,"sigmaeX");
        Common.save3DArray(sigmaeY,"sigmaeY");
        Common.save3DArray(sigmaeZ,"sigmaeZ");
        
        Common.save3DArray(sigmamX,"sigmamX");
        Common.save3DArray(sigmamY,"sigmamY");
        Common.save3DArray(sigmamZ,"sigmamZ");
        
        Common.save3DArray(index,"index");
    }
    
    public void averageAll(){
        this.averageEpsX();
        this.averageEpsY();
        this.averageEpsZ();
        
        this.averageMuX();
        this.averageMuY();
        this.averageMuZ();
        
        this.averageSigmaEX();
        this.averageSigmaEY();
        this.averageSigmaEZ();
        
        this.averageSigmaMX();
        this.averageSigmaMY();
        this.averageSigmaMZ();        
    }
    
    public void setAll(){
        this.setAllEpsRX(1.0);
        this.setAllEpsRY(1.0);
        this.setAllEpsRZ(1.0);
        
        this.setAllMuRX(1.0);
        this.setAllMuRY(1.0);
        this.setAllMuRZ(1.0);
        
        this.setAllSigmaEX(0.0);
        this.setAllSigmaEY(0.0);
        this.setAllSigmaEZ(0.0);
        
        this.setAllSigmaMX(0.0);
        this.setAllSigmaMY(0.0);
        this.setAllSigmaMZ(0.0);
    }
    
    public void setAllEpsRX(double val){
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny+1; j++){
                for(int k=0; k<nz+1; k++){         
                    this.epsrX[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllEpsRY(double val){
        for(int i=0; i<nx+1; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz+1; k++){  
                    this.epsrY[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllEpsRZ(double val){
        for(int i=0; i<nx+1; i++){
            for(int j=0; j<ny+1; j++){
                for(int k=0; k<nz; k++){  
                    this.epsrZ[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllMuRX(double val){
        for(int i=0; i<nx+1; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz; k++){  
                    this.murX[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllMuRY(double val){
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny+1; j++){
                for(int k=0; k<nz; k++){
                    this.murY[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllMuRZ(double val){
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz+1; k++){  
                    this.murZ[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllSigmaEX(double val){
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny+1; j++){
                for(int k=0; k<nz+1; k++){  
                    this.sigmaeX[i][j][k] = val;
                }
            }
        }
    }

    public void setAllSigmaEY(double val){
        for(int i=0; i<nx+1; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz+1; k++){  
                    this.sigmaeY[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllSigmaEZ(double val){
        for(int i=0; i<nx+1; i++){
            for(int j=0; j<ny+1; j++){
                for(int k=0; k<nz; k++){  
                    this.sigmaeZ[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllSigmaMX(double val){
        for(int i=0; i<nx+1; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz; k++){  
                    this.sigmamX[i][j][k] = val;
                }
            }
        }
    }

    public void setAllSigmaMY(double val){
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny+1; j++){
                for(int k=0; k<nz; k++){  
                    this.sigmamY[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllSigmaMZ(double val){
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz+1; k++){  
                    this.sigmamZ[i][j][k] = val;
                }
            }
        }
    }
    
    public double getEpsRX(int x, int y, int z){
        return this.epsrX[x][y][z];
    }
    
    public double getEpsRY(int x, int y, int z){
        return this.epsrY[x][y][z];
    }
    
    public double getEpsRZ(int x, int y, int z){
        return this.epsrZ[x][y][z];
    }
    
    public double getMuRX(int x, int y, int z){
        return this.murX[x][y][z];
    }
    
    public double getMuRY(int x, int y, int z){
        return this.murY[x][y][z];
    }
    
    public double getMuRZ(int x, int y, int z){
        return this.murZ[x][y][z];
    }
    
    public double getSigmaEX(int x, int y, int z){
        return this.sigmaeX[x][y][z];
    }

    public double getSigmaEY(int x, int y, int z){
        return this.sigmaeY[x][y][z];
    }
    
    public double getSigmaEZ(int x, int y, int z){
        return this.sigmaeZ[x][y][z];
    }
    
    public double getSigmaMX(int x, int y, int z){
        return this.sigmamX[x][y][z];
    }

    public double getSigmaMY(int x, int y, int z){
        return this.sigmamY[x][y][z];
    }
    
    public double getSigmaMZ(int x, int y, int z){
        return this.sigmamZ[x][y][z];
    }
    
    public void setEpsRX(int x, int y, int z, double val){
        this.epsrX[x][y][z] = val;
    }
    
    public void setEpsRY(int x, int y, int z, double val){
        this.epsrY[x][y][z] = val;
    }
    
    public void setEpsRZ(int x, int y, int z, double val){
        this.epsrZ[x][y][z] = val;
    }
    
    public void setMuRX(int x, int y, int z, double val){
        this.murX[x][y][z] = val;
    }
    
    public void setMuRY(int x, int y, int z, double val){
        this.murY[x][y][z] = val;
    }
    
    public void setMuRZ(int x, int y, int z, double val){
        this.murZ[x][y][z] = val;
    }
    
    public void setSigmaEX(int x, int y, int z, double val){
        this.sigmaeX[x][y][z] = val;
    }

    public void setSigmaEY(int x, int y, int z, double val){
        this.sigmaeY[x][y][z] = val;
    }
    
    public void setSigmaEZ(int x, int y, int z, double val){
        this.sigmaeZ[x][y][z] = val;
    }
    
    public void setSigmaMX(int x, int y, int z, double val){
        this.sigmamX[x][y][z] = val;
    }

    public void setSigmaMY(int x, int y, int z, double val){
        this.sigmamY[x][y][z] = val;
    }
    
    public void setSigmaMZ(int x, int y, int z, double val){
        this.sigmamZ[x][y][z] = val;
    }

    public void averageMuX(){
        int xStart = 1, xEnd = this.nx;
        int yStart = 0, yEnd = this.ny;
        int zStart = 0, zEnd = this.nz;
	
	for(int i=xStart;i<xEnd;i++){
            for(int j=yStart;j<yEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.murX[i][j][k] = (2.0*(this.m[this.index[i][j][k]].muR())*(this.m[this.index[i-1][j][k]].muR()))/(this.m[this.index[i][j][k]].muR() + this.m[this.index[i-1][j][k]].muR());
		}
            }
	}
    }

    public void averageMuY(){
        int xStart = 0, xEnd = this.nx;
        int yStart = 1, yEnd = this.ny;
        int zStart = 0, zEnd = this.nz;
        
	for(int i=xStart;i<xEnd;i++){
            for(int j=yStart;j<yEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.murY[i][j][k] = (2.0*(this.m[this.index[i][j][k]].muR())*(this.m[this.index[i][j-1][k]].muR()))/(this.m[this.index[i][j][k]].muR() + this.m[this.index[i][j-1][k]].muR());
		}
            }
	}
    }

    public void averageMuZ(){
        int xStart = 0, xEnd = this.nx;
        int yStart = 0, yEnd = this.ny;
        int zStart = 1, zEnd = this.nz;
	
	for(int i=xStart;i<xEnd;i++){
            for(int j=yStart;j<yEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.murZ[i][j][k] = (2.0*(this.m[this.index[i][j][k]].muR())*(this.m[this.index[i][j][k-1]].muR()))/(this.m[this.index[i][j][k]].muR() + this.m[this.index[i][j][k-1]].muR());
		}
            }
	}
    }

    public void averageSigmaMX(){
        int xStart = 1, xEnd = this.nx;
        int yStart = 0, yEnd = this.ny;
        int zStart = 0, zEnd = this.nz;
	
	for(int i=xStart;i<xEnd;i++){
            for(int j=yStart;j<yEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.sigmamX[i][j][k] = (2.0*(this.m[this.index[i][j][k]].sigmaM())*(this.m[this.index[i-1][j][k]].sigmaM()))/(this.m[this.index[i][j][k]].sigmaM() + this.m[this.index[i-1][j][k]].sigmaM());
		}
            }
	}
    }

    public void averageSigmaMY(){
        int xStart = 0, xEnd = this.nx;
        int yStart = 1, yEnd = this.ny;
        int zStart = 0, zEnd = this.nz;
	
	for(int i=xStart;i<xEnd;i++){
            for(int j=yStart;j<yEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.sigmamY[i][j][k] = (2.0*(this.m[this.index[i][j][k]].sigmaM())*(this.m[this.index[i][j-1][k]].sigmaM()))/(this.m[this.index[i][j][k]].sigmaM() + this.m[this.index[i][j-1][k]].sigmaM());
		}
            }
	}
    }

    public void averageSigmaMZ(){
        int xStart = 0, xEnd = this.nx;
        int yStart = 0, yEnd = this.ny;
        int zStart = 1, zEnd = this.nz;
	
	for(int i=xStart;i<xEnd;i++){
            for(int j=yStart;j<yEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.sigmamZ[i][j][k] = (2.0*(this.m[this.index[i][j][k]].sigmaM())*(this.m[this.index[i][j][k-1]].sigmaM()))/(this.m[this.index[i][j][k]].sigmaM() + this.m[this.index[i][j][k-1]].sigmaM());
		}
            }
	}
    }

    public void averageEpsX(){
        int xStart = 0, xEnd = this.nx;
        int yStart = 1, yEnd = this.ny;
        int zStart = 1, zEnd = this.nz;
	
	for(int i=xStart;i<xEnd;i++){
            for(int j=yStart;j<yEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.epsrX[i][j][k] = 0.25*(this.m[this.index[i][j][k]].epsR() + this.m[this.index[i][j-1][k]].epsR() + this.m[this.index[i][j][k-1]].epsR() + this.m[this.index[i][j-1][k-1]].epsR());
		}
            }
	}
    }

    public void averageEpsY(){
        int xStart = 1, xEnd = this.nx;
        int yStart = 0, yEnd = this.ny;
        int zStart = 1, zEnd = this.nz;
	
	for(int i=xStart;i<xEnd;i++){
            for(int j=yStart;j<yEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.epsrY[i][j][k] = 0.25*(this.m[this.index[i][j][k]].epsR() + this.m[this.index[i-1][j][k]].epsR()+ this.m[this.index[i][j][k-1]].epsR() + this.m[this.index[i-1][j][k-1]].epsR());
		}
            }
	}
    }
 
    public void averageEpsZ(){
        int xStart = 1, xEnd = this.nx;
        int yStart = 1, yEnd = this.ny;
        int zStart = 0, zEnd = this.nz;
	
	for(int i=xStart;i<xEnd;i++){
            for(int j=yStart;j<yEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.epsrZ[i][j][k] = 0.25*(this.m[this.index[i][j][k]].epsR() + this.m[this.index[i-1][j][k]].epsR() + this.m[this.index[i][j-1][k]].epsR() + this.m[this.index[i-1][j-1][k]].epsR());
		}
            }
	}
    }
 
    public void averageSigmaEX(){
        int xStart = 0, xEnd = this.nx;
        int yStart = 1, yEnd = this.ny;
        int zStart = 1, zEnd = this.nz;
	
	for(int i=xStart;i<xEnd;i++){
            for(int j=yStart;j<yEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.sigmaeX[i][j][k] = 0.25*(this.m[this.index[i][j][k]].sigmaE() + this.m[this.index[i][j-1][k]].sigmaE() + this.m[this.index[i][j][k-1]].sigmaE() + this.m[this.index[i][j-1][k-1]].sigmaE());
		}
            }
	}
    }
 
    void averageSigmaEY(){
        int xStart = 1, xEnd = this.nx;
        int yStart = 0, yEnd = this.ny;
        int zStart = 1, zEnd = this.nz;
	
	for(int i=xStart;i<xEnd;i++){
            for(int j=yStart;j<yEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.sigmaeY[i][j][k] = 0.25*(this.m[this.index[i][j][k]].sigmaE() + this.m[this.index[i-1][j][k]].sigmaE() + this.m[this.index[i][j][k-1]].sigmaE() + this.m[this.index[i-1][j][k-1]].sigmaE());
		}
            }
	}
    }
 
    public void averageSigmaEZ(){
        int xStart = 1, xEnd = this.nx;
        int yStart = 1, yEnd = this.ny;
        int zStart = 0, zEnd = this.nz;
	
	for(int i=xStart;i<xEnd;i++){
            for(int j=yStart;j<yEnd;j++){
                for(int k=zStart;k<zEnd;k++){
                    this.sigmaeZ[i][j][k] = 0.25*(this.m[this.index[i][j][k]].sigmaE() + this.m[this.index[i-1][j][k]].sigmaE() + this.m[this.index[i][j-1][k]].sigmaE() + this.m[this.index[i-1][j-1][k]].sigmaE());
		}
            }
	}
    }
}
