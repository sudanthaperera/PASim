package pasim;

public class HField {
    private double[][][] Hx,Hy,Hz;
    private double[][][] Chxh,Chxez,Chxey,Chyh,Chyex,Chyez,Chzh,Chzey,Chzex;
    private int nx, ny, nz;
    private Cell c;
    private double dt;
    private Boundary b;
    
    public HField(int nx,int ny,int nz,Cell c){
        this.nx = nx;
        this.ny = ny;
        this.nz = nz;
        this.c = c;
        this.dt = c.getDeltaT();
        
        Hx = new double[nx+1][ny][nz];
        Hy = new double[nx][ny+1][nz];
        Hz = new double[nx][ny][nz+1];
	
        Chxh =  Common.genDouble3DArray(nx+1, ny, nz, 0.0);
        Chxez =  Common.genDouble3DArray(nx+1, ny, nz, 0.0);
        Chxey =  Common.genDouble3DArray(nx+1, ny, nz, 0.0);
        Chyh =  Common.genDouble3DArray(nx, ny+1, nz, 0.0);
        Chyex =  Common.genDouble3DArray(nx, ny+1, nz, 0.0);
        Chyez =  Common.genDouble3DArray(nx, ny+1, nz, 0.0);
        Chzh =  Common.genDouble3DArray(nx, ny, nz+1, 0.0);
        Chzey =  Common.genDouble3DArray(nx, ny, nz+1, 0.0);
        Chzex =  Common.genDouble3DArray(nx, ny, nz+1, 0.0);
    }
    
    public void setBoundary(Boundary b){
        this.b=b;
    }
    
    public double getChxh(int x,int y, int z){
        return this.Chxh[x][y][z];
    }
    
    public double getChxez(int x,int y, int z){
        return this.Chxez[x][y][z];
    }
    
    public double getChxey(int x,int y, int z){
        return this.Chxey[x][y][z];
    }
    
    public double getChyh(int x,int y, int z){
        return this.Chyh[x][y][z];
    }
    
    public double getChyex(int x,int y, int z){
        return this.Chyex[x][y][z];
    }
    
    public double getChyez(int x,int y, int z){
        return this.Chyez[x][y][z];
    }
    
    public double getChzh(int x,int y, int z){
        return this.Chzh[x][y][z];
    }
    
    public double getChzey(int x,int y, int z){
        return this.Chzey[x][y][z];
    }
    
    public double getChzex(int x,int y, int z){
        return this.Chzex[x][y][z];
    }
    
    public void setChxh(int x,int y, int z, double value){
        this.Chxh[x][y][z] = value;
    }
    
    public void setChxez(int x,int y, int z, double value){
        this.Chxez[x][y][z] = value;
    }
    
    public void setChxey(int x,int y, int z, double value){
        this.Chxey[x][y][z] = value;
    }
    
    public void setChyh(int x,int y, int z, double value){
        this.Chyh[x][y][z] = value;
    }
    
    public void setChyex(int x,int y, int z, double value){
        this.Chyex[x][y][z] = value;
    }
    
    public void setChyez(int x,int y, int z, double value){
        this.Chyez[x][y][z] = value;
    }
    
    public void setChzh(int x,int y, int z, double value){
        this.Chzh[x][y][z] = value;
    }
    
    public void setChzey(int x,int y, int z, double value){
        this.Chzey[x][y][z] = value;
    }
    
    public void setChzex(int x,int y, int z, double value){
        this.Chzex[x][y][z] = value;
    }
    
    public void setHX(int xIndex, int yIndex, int zIndex, double val){
        Hx[xIndex][yIndex][zIndex] = val;
    }

    public void setHY(int xIndex, int yIndex, int zIndex, double val){
        Hy[xIndex][yIndex][zIndex] = val;
    }
    
    public void setHZ(int xIndex, int yIndex, int zIndex, double val){
        Hz[xIndex][yIndex][zIndex] = val;
    }
    
    public double getHX(int xIndex, int yIndex, int zIndex){
        return Hx[xIndex][yIndex][zIndex];
    }

    public double getHY(int xIndex, int yIndex, int zIndex){
        return Hy[xIndex][yIndex][zIndex];
    }
    
    public double getHZ(int xIndex, int yIndex, int zIndex){
        return Hz[xIndex][yIndex][zIndex];
    }
    
    public void updatingCoefficients(ProblemSpace ps){
        
        MaterialGrid mg = ps.getMaterialGrid();
 	for(int i=0;i<nx+1;i++){
	    for(int j=0;j<ny;j++){
		for(int k=0;k<nz;k++){
		    Chxh[i][j][k]  =  (2*mg.getMuRX(i, j, k)*Constants.MU0 - dt*mg.getSigmaMX(i, j, k))/(2*mg.getMuRX(i, j, k)*Constants.MU0 + dt*mg.getSigmaMX(i, j, k));
		    Chxez[i][j][k] = -(2*dt/c.getDeltaY())/(2*mg.getMuRX(i, j, k)*Constants.MU0 + dt*mg.getSigmaMX(i, j, k));
		    Chxey[i][j][k] =  (2*dt/c.getDeltaZ())/(2*mg.getMuRX(i, j, k)*Constants.MU0 + dt*mg.getSigmaMX(i, j, k));
		}
	    }
	}
	
	for(int i=0;i<nx;i++){
	    for(int j=0;j<ny+1;j++){
		for(int k=0;k<nz;k++){
		    Chyh[i][j][k]  =  (2*mg.getMuRY(i, j, k)*Constants.MU0 - dt*mg.getSigmaMY(i, j, k))/(2*mg.getMuRY(i, j, k)*Constants.MU0 + dt*mg.getSigmaMY(i, j, k));
		    Chyex[i][j][k] = -(2*dt/c.getDeltaZ())/(2*mg.getMuRY(i, j, k)*Constants.MU0 + dt*mg.getSigmaMY(i, j, k));
		    Chyez[i][j][k] =  (2*dt/c.getDeltaX())/(2*mg.getMuRY(i, j, k)*Constants.MU0 + dt*mg.getSigmaMY(i, j, k));
		}
	    }
	}
	
	for(int i=0;i<nx;i++){
	    for(int j=0;j<ny;j++){
		for(int k=0;k<nz+1;k++){
		    Chzh[i][j][k]  =  (2*mg.getMuRZ(i, j, k)*Constants.MU0 - dt*mg.getSigmaMZ(i, j, k))/(2*mg.getMuRZ(i, j, k)*Constants.MU0 + dt*mg.getSigmaMZ(i, j, k));
		    Chzey[i][j][k] = -(2*dt/c.getDeltaX())/(2*mg.getMuRZ(i, j, k)*Constants.MU0 + dt*mg.getSigmaMZ(i, j, k));
		    Chzex[i][j][k] =  (2*dt/c.getDeltaY())/(2*mg.getMuRZ(i, j, k)*Constants.MU0 + dt*mg.getSigmaMZ(i, j, k));
		}
	    }
	}       
    }
    
    public void setAllX(double val){
        for(int i=0; i<nx+1; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz; k++){ 
                    Hx[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllY(double val){
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny+1; j++){
                for(int k=0; k<nz; k++){ 
                    Hy[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllZ(double val){
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz+1; k++){
                    Hz[i][j][k] = val;
                }
            }
        }
    }
    
    public void updateHField(EField E){
	for (int i=0; i<nx+1; i++){
	    for (int j=0; j<ny; j++){
		for (int k=0; k<nz; k++){	   
		    Hx[i][j][k] = Chxh[i][j][k]*Hx[i][j][k] + Chxey[i][j][k]*(E.getEY(i,j,k+1) - E.getEY(i,j,k)) + Chxez[i][j][k]*(E.getEZ(i,j+1,k) - E.getEZ(i,j,k));
		}
	    }
	}
	for (int i=0; i<nx; i++){
	    for (int j=0; j<ny+1; j++){
		for (int k=0; k<nz; k++){
		    Hy[i][j][k] = Chyh[i][j][k]*Hy[i][j][k] + Chyez[i][j][k]*(E.getEZ(i+1,j,k) - E.getEZ(i,j,k)) + Chyex[i][j][k]*(E.getEX(i,j,k+1) - E.getEX(i,j,k));
		}
	    }
	}
	for (int i=0; i<nx; i++){
	    for (int j=0; j<ny; j++){
		for (int k=0; k<nz+1; k++){
		    Hz[i][j][k] = Chzh[i][j][k]*Hz[i][j][k] + Chzex[i][j][k]*(E.getEX(i,j+1,k) - E.getEX(i,j,k)) + Chzey[i][j][k]*(E.getEY(i+1,j,k) - E.getEY(i,j,k));
		}
	    }
	}
    }
}
