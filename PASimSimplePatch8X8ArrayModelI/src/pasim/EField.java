package pasim;

public class EField{
    private double[][][] Ex,Ey,Ez;
    private double[][][] Cexe,Cexhz,Cexhy,Ceye,Ceyhx,Ceyhz,Ceze,Cezhy,Cezhx;
    private int nx, ny, nz;
    private Cell c;
    private double dt;
    private Boundary b;
    
    public EField(int nx,int ny,int nz,Cell c){
        this.nx = nx;
        this.ny = ny;
        this.nz = nz;
        this.c = c;
        this.dt = c.getDeltaT();
        
        Ex = new double[nx][ny+1][nz+1];
        Ey = new double[nx+1][ny][nz+1];
        Ez = new double[nx+1][ny+1][nz];
        
        Cexe =  Common.genDouble3DArray(nx,ny+1,nz+1,0.0);
        Cexhz =  Common.genDouble3DArray(nx,ny+1,nz+1,0.0);
        Cexhy =  Common.genDouble3DArray(nx,ny+1,nz+1,0.0);
        Ceye =  Common.genDouble3DArray(nx+1,ny,nz+1,0.0);
        Ceyhx =  Common.genDouble3DArray(nx+1,ny,nz+1,0.0);
        Ceyhz =  Common.genDouble3DArray(nx+1,ny,nz+1,0.0);
        Ceze =  Common.genDouble3DArray(nx+1,ny+1,nz,0.0);
        Cezhy =  Common.genDouble3DArray(nx+1,ny+1,nz,0.0);
        Cezhx =  Common.genDouble3DArray(nx+1,ny+1,nz,0.0);
    }
    
    public void setBoundary(Boundary b){
        this.b=b;
    }
    
    public double[][] getEFieldArray(int d){
        return Ez[d];
    }
    
    public double getCexe(int x,int y, int z){
        return this.Cexe[x][y][z];
    }
    
    public double getCexhz(int x,int y, int z){
        return this.Cexhz[x][y][z];
    }
    
    public double getCexhy(int x,int y, int z){
        return this.Cexhy[x][y][z];
    }
    
    public double getCeye(int x,int y, int z){
        return this.Ceye[x][y][z];
    }
    
    public double getCeyhx(int x,int y, int z){
        return this.Ceyhx[x][y][z];
    }
    
    public double getCeyhz(int x,int y, int z){
        return this.Ceyhz[x][y][z];
    }
    
    public double getCeze(int x,int y, int z){
        return this.Ceze[x][y][z];
    }
    
    public double getCezhy(int x,int y, int z){
        return this.Cezhy[x][y][z];
    }
    
    public double getCezhx(int x,int y, int z){
        return this.Cezhx[x][y][z];
    }
    
    public void setCexe(int x,int y, int z, double value){
        this.Cexe[x][y][z] = value;
    }
    
    public void setCexhz(int x,int y, int z, double value){
        this.Cexhz[x][y][z] = value;
    }
    
    public void setCexhy(int x,int y, int z, double value){
        this.Cexhy[x][y][z] = value;
    }
    
    public void setCeye(int x,int y, int z, double value){
        this.Ceye[x][y][z] = value;
    }
    
    public void setCeyhx(int x,int y, int z, double value){
        this.Ceyhx[x][y][z] = value;
    }
    
    public void setCeyhz(int x,int y, int z, double value){
        this.Ceyhz[x][y][z] = value;
    }
    
    public void setCeze(int x,int y, int z, double value){
        this.Ceze[x][y][z] = value;
    }
    
    public void setCezhy(int x,int y, int z, double value){
        this.Cezhy[x][y][z] = value;
    }
    
    public void setCezhx(int x,int y, int z, double value){
        this.Cezhx[x][y][z] = value;
    }
        
    public void setEX(int xIndex, int yIndex, int zIndex, double val){
        Ex[xIndex][yIndex][zIndex] = val;
    }

    public void setEY(int xIndex, int yIndex, int zIndex, double val){
        Ey[xIndex][yIndex][zIndex] = val;
    }
    
    public void setEZ(int xIndex, int yIndex, int zIndex, double val){
        Ez[xIndex][yIndex][zIndex] = val;
    }
    
    public double getEX(int xIndex, int yIndex, int zIndex){
        return Ex[xIndex][yIndex][zIndex];
    }

    public double getEY(int xIndex, int yIndex, int zIndex){
        return Ey[xIndex][yIndex][zIndex];
    }
    
    public double getEZ(int xIndex, int yIndex, int zIndex){
        return Ez[xIndex][yIndex][zIndex];
    }
    
    public void updatingCoefficients(ProblemSpace ps){
        
        MaterialGrid mg = ps.getMaterialGrid();
	for(int i=0;i<nx;i++){
            for(int j=0;j<ny+1;j++){
		for(int k=0;k<nz+1;k++){
                    Cexe[i][j][k]  =  (2*mg.getEpsRX(i, j, k)*Constants.EPS0 - dt*mg.getSigmaEX(i, j, k))/(2*mg.getEpsRX(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEX(i, j, k));
                    Cexhz[i][j][k]=  (2*dt/c.getDeltaY())/(2*mg.getEpsRX(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEX(i, j, k));
                    Cexhy[i][j][k]= -(2*dt/c.getDeltaZ())/(2*mg.getEpsRX(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEX(i, j, k));
		}
            }
	}
	
	for(int i=0;i<nx+1;i++){
            for(int j=0;j<ny;j++){
                for(int k=0;k<nz+1;k++){
                    Ceye[i][j][k]  =  (2*mg.getEpsRY(i, j, k)*Constants.EPS0 - dt*mg.getSigmaEY(i, j, k))/(2*mg.getEpsRY(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEY(i, j, k));
                    Ceyhx[i][j][k] =  (2*dt/c.getDeltaZ())/(2*mg.getEpsRY(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEY(i, j, k));
                    Ceyhz[i][j][k] = -(2*dt/c.getDeltaX())/(2*mg.getEpsRY(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEY(i, j, k));
		}
            }
	}
	
	for(int i=0;i<nx+1;i++){
            for(int j=0;j<ny+1;j++){
		for(int k=0;k<nz;k++){
                    Ceze[i][j][k]  =  (2*mg.getEpsRZ(i, j, k)*Constants.EPS0 - dt*mg.getSigmaEZ(i, j, k))/(2*mg.getEpsRZ(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEZ(i, j, k));
                    Cezhy[i][j][k] =  (2*dt/c.getDeltaX())/(2*mg.getEpsRZ(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEZ(i, j, k));
                    Cezhx[i][j][k] = -(2*dt/c.getDeltaY())/(2*mg.getEpsRZ(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEZ(i, j, k));
		}
            }
	}        
    }
    
    public void setAllX(double val){
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny+1; j++){
                for(int k=0; k<nz+1; k++){ 
                    Ex[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllY(double val){
        for(int i=0; i<nx+1; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz+1; k++){ 
                    Ey[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllZ(double val){
        for(int i=0; i<nx+1; i++){
            for(int j=0; j<ny+1; j++){
                for(int k=0; k<nz; k++){
                    Ez[i][j][k] = val;
                }
            }
        }
    }
    
    public void updateEField(HField H){
	for (int i=0;i<nx;i++){
	    for (int j=1;j<ny;j++){
		for (int k=1;k<nz;k++){
		    Ex[i][j][k] = Cexe[i][j][k]*Ex[i][j][k] + Cexhz[i][j][k]*(H.getHZ(i,j,k) - H.getHZ(i,j-1,k)) + Cexhy[i][j][k]*(H.getHY(i,j,k) - H.getHY(i,j,k-1));
		}
	    }
	}
	for (int i=1;i<nx;i++){
	    for (int j=0;j<ny;j++){
		for (int k=1;k<nz;k++){
		    Ey[i][j][k] = Ceye[i][j][k]*Ey[i][j][k] + Ceyhx[i][j][k]*(H.getHX(i,j,k) - H.getHX(i,j,k-1)) + Ceyhz[i][j][k]*(H.getHZ(i,j,k) - H.getHZ(i-1,j,k)); 
		}
	    }
	}
	for (int i=1;i<nx;i++){
	    for (int j=1;j<ny;j++){
		for (int k=0;k<nz;k++){
		    Ez[i][j][k] = Ceze[i][j][k]*Ez[i][j][k] + Cezhy[i][j][k]*(H.getHY(i,j,k) - H.getHY(i-1,j,k)) + Cezhx[i][j][k]*(H.getHX(i,j,k) - H.getHX(i,j-1,k));
		}
	    }
	}	
    }
}
