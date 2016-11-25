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
public class HField extends Field{
    
    public HField(int nx,int ny,int nz,Cell c){
        super(nx,ny,nz,c);
        Hx = new double[nx+1][ny][nz];
        Hy = new double[nx][ny+1][nz];
        Hz = new double[nx][ny][nz+1];
	
        Chxh = Common.genDouble3DArray(nx+1, ny, nz, 0.0);
        Chxez = Common.genDouble3DArray(nx+1, ny, nz, 0.0);
        Chxey = Common.genDouble3DArray(nx+1, ny, nz, 0.0);
        Chyh = Common.genDouble3DArray(nx, ny+1, nz, 0.0);
        Chyex = Common.genDouble3DArray(nx, ny+1, nz, 0.0);
        Chyez = Common.genDouble3DArray(nx, ny+1, nz, 0.0);
        Chzh = Common.genDouble3DArray(nx, ny, nz+1, 0.0);
        Chzey = Common.genDouble3DArray(nx, ny, nz+1, 0.0);
        Chzex = Common.genDouble3DArray(nx, ny, nz+1, 0.0);
    }
    
    public void updatingCoefficients(ProblemSpace ps){
        System.out.println("General magnetic field updating coefficients....");
        MaterialGrid mg = ps.getMaterialGrid();
 	for(int i=0;i<nx+1;i++){
	    for(int j=0;j<ny;j++){
		for(int k=0;k<nz;k++){
		    //Coeffiecients updating Hx
		    Chxh[i][j][k]  =  (2*mg.getMuRX(i, j, k)*Constants.MU0 - dt*mg.getSigmaMX(i, j, k))/(2*mg.getMuRX(i, j, k)*Constants.MU0 + dt*mg.getSigmaMX(i, j, k));
		    Chxez[i][j][k] = -(2*dt/c.getDeltaY())/(2*mg.getMuRX(i, j, k)*Constants.MU0 + dt*mg.getSigmaMX(i, j, k));
		    Chxey[i][j][k] =  (2*dt/c.getDeltaZ())/(2*mg.getMuRX(i, j, k)*Constants.MU0 + dt*mg.getSigmaMX(i, j, k));
		}
	    }
	}
	
	for(int i=0;i<nx;i++){
	    for(int j=0;j<ny+1;j++){
		for(int k=0;k<nz;k++){
		    //Coeffiecients updating Hy
		    Chyh[i][j][k]  =  (2*mg.getMuRY(i, j, k)*Constants.MU0 - dt*mg.getSigmaMY(i, j, k))/(2*mg.getMuRY(i, j, k)*Constants.MU0 + dt*mg.getSigmaMY(i, j, k));
		    Chyex[i][j][k] = -(2*dt/c.getDeltaZ())/(2*mg.getMuRY(i, j, k)*Constants.MU0 + dt*mg.getSigmaMY(i, j, k));
		    Chyez[i][j][k] =  (2*dt/c.getDeltaX())/(2*mg.getMuRY(i, j, k)*Constants.MU0 + dt*mg.getSigmaMY(i, j, k));
		}
	    }
	}
	
	for(int i=0;i<nx;i++){
	    for(int j=0;j<ny;j++){
		for(int k=0;k<nz+1;k++){
		    //Coeffiecients updating Hz
		    Chzh[i][j][k]  =  (2*mg.getMuRZ(i, j, k)*Constants.MU0 - dt*mg.getSigmaMZ(i, j, k))/(2*mg.getMuRZ(i, j, k)*Constants.MU0 + dt*mg.getSigmaMZ(i, j, k));
		    Chzey[i][j][k] = -(2*dt/c.getDeltaX())/(2*mg.getMuRZ(i, j, k)*Constants.MU0 + dt*mg.getSigmaMZ(i, j, k));
		    Chzex[i][j][k] =  (2*dt/c.getDeltaY())/(2*mg.getMuRZ(i, j, k)*Constants.MU0 + dt*mg.getSigmaMZ(i, j, k));
		}
	    }
	}       
    }
    
    public double getHx(int i, int j, int k){
        return Hx[i][j][k];
    }
    
    public double getHy(int i, int j, int k){
        return Hy[i][j][k];
    }
        
    public double getHz(int i, int j, int k){
        return Hz[i][j][k];
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
    
    public void updateHField(){
	for (int i=0; i<nx+1; i++){
	    for (int j=0; j<ny; j++){
		for (int k=0; k<nz; k++){	   
		    Hx[i][j][k] = Chxh[i][j][k]*Hx[i][j][k] + Chxey[i][j][k]*(Ey[i][j][k+1] - Ey[i][j][k]) + Chxez[i][j][k]*(Ez[i][j+1][k] - Ez[i][j][k]);
		}
	    }
	}
	for (int i=0; i<nx; i++){
	    for (int j=0; j<ny+1; j++){
		for (int k=0; k<nz; k++){
		    Hy[i][j][k] = Chyh[i][j][k]*Hy[i][j][k] + Chyez[i][j][k]*(Ez[i+1][j][k] - Ez[i][j][k]) + Chyex[i][j][k]*(Ex[i][j][k+1] - Ex[i][j][k]);
		}
	    }
	}
	for (int i=0; i<nx; i++){
	    for (int j=0; j<ny; j++){
		for (int k=0; k<nz+1; k++){
		    Hz[i][j][k] = Chzh[i][j][k]*Hz[i][j][k] + Chzex[i][j][k]*(Ex[i][j+1][k] - Ex[i][j][k]) + Chzey[i][j][k]*(Ey[i+1][j][k] - Ey[i][j][k]);
		}
	    }
	}
    }
}
