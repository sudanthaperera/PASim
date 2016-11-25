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
public class EField extends Field {
    
    public EField(int nx,int ny,int nz,Cell c){
        super(nx,ny,nz,c);
        Ex = new double[nx][ny+1][nz+1];
        Ey = new double[nx+1][ny][nz+1];
        Ez = new double[nx+1][ny+1][nz];
        
        Cexe = Common.genDouble3DArray(nx,ny+1,nz+1,0.0);
        Cexhz = Common.genDouble3DArray(nx,ny+1,nz+1,0.0);
        Cexhy = Common.genDouble3DArray(nx,ny+1,nz+1,0.0);
        Ceye = Common.genDouble3DArray(nx+1,ny,nz+1,0.0);
        Ceyhx = Common.genDouble3DArray(nx+1,ny,nz+1,0.0);
        Ceyhz = Common.genDouble3DArray(nx+1,ny,nz+1,0.0);
        Ceze = Common.genDouble3DArray(nx+1,ny+1,nz,0.0);
        Cezhy = Common.genDouble3DArray(nx+1,ny+1,nz,0.0);
        Cezhx = Common.genDouble3DArray(nx+1,ny+1,nz,0.0);
    }
    
    public double[][] getEFieldArray(int d){
        return Ez[d];
    }
    
    public void updatingCoefficients(ProblemSpace ps){
        MaterialGrid mg = ps.getMaterialGrid();
	for(int i=0;i<nx;i++){
            for(int j=0;j<ny+1;j++){
		for(int k=0;k<nz+1;k++){
                    // Coeffiecients updating Ex
                    Cexe[i][j][k] =  (2*mg.getEpsRX(i, j, k)*Constants.EPS0 - dt*mg.getSigmaEX(i, j, k))/(2*mg.getEpsRX(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEX(i, j, k));
                    Cexhz[i][j][k]=  (2*dt/c.getDeltaY())/(2*mg.getEpsRX(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEX(i, j, k));
                    Cexhy[i][j][k]= -(2*dt/c.getDeltaZ())/(2*mg.getEpsRX(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEX(i, j, k));
		}
            }
	}
	
	for(int i=0;i<nx+1;i++){
            for(int j=0;j<ny;j++){
                for(int k=0;k<nz+1;k++){
                    //Coeffiecients updating Ey
                    Ceye[i][j][k]  =  (2*mg.getEpsRY(i, j, k)*Constants.EPS0 - dt*mg.getSigmaEY(i, j, k))/(2*mg.getEpsRY(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEY(i, j, k));
                    Ceyhx[i][j][k] =  (2*dt/c.getDeltaZ())/(2*mg.getEpsRY(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEY(i, j, k));
                    Ceyhz[i][j][k] = -(2*dt/c.getDeltaX())/(2*mg.getEpsRY(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEY(i, j, k));
		}
            }
	}
	
	for(int i=0;i<nx+1;i++){
            for(int j=0;j<ny+1;j++){
		for(int k=0;k<nz;k++){
                    //Coeffiecients updating Ez
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
    
    public void updateEField(){
	for (int i=0;i<nx;i++){
	    for (int j=1;j<ny;j++){
		for (int k=1;k<nz;k++){
		    Ex[i][j][k] = Cexe[i][j][k]*Ex[i][j][k] + Cexhz[i][j][k]*(Hz[i][j][k] - Hz[i][j-1][k]) + Cexhy[i][j][k]*(Hy[i][j][k] - Hy[i][j][k-1]);   
		}
	    }
	}
	for (int i=1;i<nx;i++){
	    for (int j=0;j<ny;j++){
		for (int k=1;k<nz;k++){
		    Ey[i][j][k] = Ceye[i][j][k]*Ey[i][j][k] + Ceyhx[i][j][k]*(Hx[i][j][k] - Hx[i][j][k-1]) + Ceyhz[i][j][k]*(Hz[i][j][k] - Hz[i-1][j][k]); 
		}
	    }
	}
	for (int i=1;i<nx;i++){
	    for (int j=1;j<ny;j++){
		for (int k=0;k<nz;k++){
		    Ez[i][j][k] = Ceze[i][j][k]*Ez[i][j][k] + Cezhy[i][j][k]*(Hy[i][j][k] - Hy[i-1][j][k]) + Cezhx[i][j][k]*(Hx[i][j][k] - Hx[i][j-1][k]);
		}
	    }
	}	
    }

    public void updateEFieldWithPBC_V2(){
        int loc=0;
	for (int i=0;i<nx;i++){
	    for (int j=1;j<ny;j++){
		for (int k=1;k<nz;k++){
                    if (b.getPBCY() && j == ny-1 && k != nz-1) loc = 1;
                    else if (b.getPBCZ() && j != ny-1 && k == nz-1) loc = 2;
                    else if (b.getPBCY() && b.getPBCZ() && j == ny-1 && k == nz-1) loc = 3;
                    else loc = 0;
                    
                    switch (loc){
                        case 0:
                            Ex[i][j][k] = Cexe[i][j][k]*Ex[i][j][k] + Cexhz[i][j][k]*(Hz[i][j][k] - Hz[i][j-1][k]) + Cexhy[i][j][k]*(Hy[i][j][k] - Hy[i][j][k-1]);
                            break;
                        case 1:
                            Ex[i][j][k] = Cexe[i][j][k]*Ex[i][j][k] + Cexhz[i][j][k]*(Hz[i][0][k] - Hz[i][j-1][k]) + Cexhy[i][j][k]*(Hy[i][j][k] - Hy[i][j][k-1]);
                            Ex[i][0][k] = Ex[i][j][k];
                            break;
                        case 2:
                            Ex[i][j][k] = Cexe[i][j][k]*Ex[i][j][k] + Cexhz[i][j][k]*(Hz[i][j][k] - Hz[i][j-1][k]) + Cexhy[i][j][k]*(Hy[i][j][0] - Hy[i][j][k-1]);
                            Ex[i][j][0] = Ex[i][j][k];
                            break;
                        case 3:
                            Ex[i][j][k] = Cexe[i][j][k]*Ex[i][j][k] + Cexhz[i][j][k]*(Hz[i][0][k] - Hz[i][j-1][k]) + Cexhy[i][j][k]*(Hy[i][j][0] - Hy[i][j][k-1]);
                            Ex[i][0][0] = Ex[i][j][k];
                            break;
                    }
		}
	    }
	}
	for (int i=1;i<nx;i++){
	    for (int j=0;j<ny;j++){
		for (int k=1;k<nz;k++){
                    if (b.getPBCX() && i == nx-1 && k != nz-1) loc = 1;
                    else if (b.getPBCZ() && i != nx-1 && k == nz-1) loc = 2;
                    else if (b.getPBCX() && b.getPBCZ() && i == nx-1 && k == nz-1) loc = 3;
                    else loc = 0;
                    
                    switch (loc){
                        case 0:
                            Ey[i][j][k] = Ceye[i][j][k]*Ey[i][j][k] + Ceyhx[i][j][k]*(Hx[i][j][k] - Hx[i][j][k-1]) + Ceyhz[i][j][k]*(Hz[i][j][k] - Hz[i-1][j][k]);
                            break;
                        case 1:
                            Ey[i][j][k] = Ceye[i][j][k]*Ey[i][j][k] + Ceyhx[i][j][k]*(Hx[i][j][k] - Hx[i][j][k-1]) + Ceyhz[i][j][k]*(Hz[0][j][k] - Hz[i-1][j][k]);
                            Ey[0][j][k] = Ey[i][j][k];
                            break;
                        case 2:
                            Ey[i][j][k] = Ceye[i][j][k]*Ey[i][j][k] + Ceyhx[i][j][k]*(Hx[i][j][0] - Hx[i][j][k-1]) + Ceyhz[i][j][k]*(Hz[i][j][k] - Hz[i-1][j][k]);
                            Ey[i][j][0] = Ey[i][j][k];
                            break;
                        case 3:
                            Ey[i][j][k] = Ceye[i][j][k]*Ey[i][j][k] + Ceyhx[i][j][k]*(Hx[i][j][0] - Hx[i][j][k-1]) + Ceyhz[i][j][k]*(Hz[0][j][k] - Hz[i-1][j][k]);
                            Ey[0][j][0] = Ey[i][j][k];
                    }
		}
	    }
	}
	for (int i=1;i<nx;i++){
	    for (int j=1;j<ny;j++){
		for (int k=0;k<nz;k++){
                    if (b.getPBCX() && i == nx-1 && j != ny-1) loc = 1;
                    else if (b.getPBCY() && i != nx-1 && j == ny-1) loc = 2;
                    else if (b.getPBCX() && b.getPBCY() && i == nx-1 && j == ny-1) loc = 3;
                    else loc = 0;
                    
                    switch (loc){
                        case 0:
                            Ez[i][j][k] = Ceze[i][j][k]*Ez[i][j][k] + Cezhy[i][j][k]*(Hy[i][j][k] - Hy[i-1][j][k]) + Cezhx[i][j][k]*(Hx[i][j][k] - Hx[i][j-1][k]);
                            break;
                        case 1:
                            Ez[i][j][k] = Ceze[i][j][k]*Ez[i][j][k] + Cezhy[i][j][k]*(Hy[0][j][k] - Hy[i-1][j][k]) + Cezhx[i][j][k]*(Hx[i][j][k] - Hx[i][j-1][k]);
                            Ez[0][j][k] = Ez[i][j][k];
                            break;
                        case 2:
                            Ez[i][j][k] = Ceze[i][j][k]*Ez[i][j][k] + Cezhy[i][j][k]*(Hy[i][j][k] - Hy[i-1][j][k]) + Cezhx[i][j][k]*(Hx[i][0][k] - Hx[i][j-1][k]);
                            Ez[i][0][k] = Ez[i][j][k];
                            break;
                        case 3:
                            Ez[i][j][k] = Ceze[i][j][k]*Ez[i][j][k] + Cezhy[i][j][k]*(Hy[0][j][k] - Hy[i-1][j][k]) + Cezhx[i][j][k]*(Hx[i][0][k] - Hx[i][j-1][k]);
                            Ez[0][0][k] = Ez[i][j][k];
                    }
		}
	    }
	}	
    }
    

    public void updateEFieldWithPBC(){
	for (int i=0;i<nx;i++){
	    for (int j=1;j<ny;j++){
		for (int k=1;k<nz;k++){
                    if (b.getPBCY() && j == ny-1 && k != nz-1) {
                        Ex[i][j][k] = Cexe[i][j][k]*Ex[i][j][k] + Cexhz[i][j][k]*(Hz[i][0][k] - Hz[i][j-1][k]) + Cexhy[i][j][k]*(Hy[i][j][k] - Hy[i][j][k-1]);
                        Ex[i][0][k] = Ex[i][j][k];
                    }
                    else if (b.getPBCZ() && j != ny-1 && k == nz-1){
                        Ex[i][j][k] = Cexe[i][j][k]*Ex[i][j][k] + Cexhz[i][j][k]*(Hz[i][j][k] - Hz[i][j-1][k]) + Cexhy[i][j][k]*(Hy[i][j][0] - Hy[i][j][k-1]);
                        Ex[i][j][0] = Ex[i][j][k];
                    }
                    else if (b.getPBCY() && b.getPBCZ() && j == ny-1 && k == nz-1){
                        Ex[i][j][k] = Cexe[i][j][k]*Ex[i][j][k] + Cexhz[i][j][k]*(Hz[i][0][k] - Hz[i][j-1][k]) + Cexhy[i][j][k]*(Hy[i][j][0] - Hy[i][j][k-1]);
                        Ex[i][0][0] = Ex[i][j][k];
                    }
                    else{
                        Ex[i][j][k] = Cexe[i][j][k]*Ex[i][j][k] + Cexhz[i][j][k]*(Hz[i][j][k] - Hz[i][j-1][k]) + Cexhy[i][j][k]*(Hy[i][j][k] - Hy[i][j][k-1]);
                    }
		}
	    }
	}
	for (int i=1;i<nx;i++){
	    for (int j=0;j<ny;j++){
		for (int k=1;k<nz;k++){
                    if (b.getPBCX() && i == nx-1 && k != nz-1){
                        Ey[i][j][k] = Ceye[i][j][k]*Ey[i][j][k] + Ceyhx[i][j][k]*(Hx[i][j][k] - Hx[i][j][k-1]) + Ceyhz[i][j][k]*(Hz[0][j][k] - Hz[i-1][j][k]);
                        Ey[0][j][k] = Ey[i][j][k];
                    }
                    else if (b.getPBCZ() && i != nx-1 && k == nz-1){
                        Ey[i][j][k] = Ceye[i][j][k]*Ey[i][j][k] + Ceyhx[i][j][k]*(Hx[i][j][0] - Hx[i][j][k-1]) + Ceyhz[i][j][k]*(Hz[i][j][k] - Hz[i-1][j][k]);
                        Ey[i][j][0] = Ey[i][j][k];
                    }
                    else if (b.getPBCX() && b.getPBCZ() && i == nx-1 && k == nz-1){
                        Ey[i][j][k] = Ceye[i][j][k]*Ey[i][j][k] + Ceyhx[i][j][k]*(Hx[i][j][0] - Hx[i][j][k-1]) + Ceyhz[i][j][k]*(Hz[0][j][k] - Hz[i-1][j][k]);
                        Ey[0][j][0] = Ey[i][j][k];
                    }
                    else{
                        Ey[i][j][k] = Ceye[i][j][k]*Ey[i][j][k] + Ceyhx[i][j][k]*(Hx[i][j][k] - Hx[i][j][k-1]) + Ceyhz[i][j][k]*(Hz[i][j][k] - Hz[i-1][j][k]);
                    }
		}
	    }
	}
	for (int i=1;i<nx;i++){
	    for (int j=1;j<ny;j++){
		for (int k=0;k<nz;k++){
                    if (b.getPBCX() && i == nx-1 && j != ny-1){
                        Ez[i][j][k] = Ceze[i][j][k]*Ez[i][j][k] + Cezhy[i][j][k]*(Hy[0][j][k] - Hy[i-1][j][k]) + Cezhx[i][j][k]*(Hx[i][j][k] - Hx[i][j-1][k]);
                        Ez[0][j][k] = Ez[i][j][k];
                    }
                    else if (b.getPBCY() && i != nx-1 && j == ny-1){
                        Ez[i][j][k] = Ceze[i][j][k]*Ez[i][j][k] + Cezhy[i][j][k]*(Hy[i][j][k] - Hy[i-1][j][k]) + Cezhx[i][j][k]*(Hx[i][0][k] - Hx[i][j-1][k]);
                        Ez[i][0][k] = Ez[i][j][k];
                    }
                    else if (b.getPBCX() && b.getPBCY() && i == nx-1 && j == ny-1){
                        Ez[i][j][k] = Ceze[i][j][k]*Ez[i][j][k] + Cezhy[i][j][k]*(Hy[0][j][k] - Hy[i-1][j][k]) + Cezhx[i][j][k]*(Hx[i][0][k] - Hx[i][j-1][k]);
                        Ez[0][0][k] = Ez[i][j][k];
                    }
                    else{
                        Ez[i][j][k] = Ceze[i][j][k]*Ez[i][j][k] + Cezhy[i][j][k]*(Hy[i][j][k] - Hy[i-1][j][k]) + Cezhx[i][j][k]*(Hx[i][j][k] - Hx[i][j-1][k]);
                    }
		}
	    }
	}	
    }
    
    
    public double getEx(int i, int j, int k){
        return Ex[i][j][k];
    }
    
    public double getEy(int i, int j, int k){
        return Ey[i][j][k];
    }
        
    public double getEz(int i, int j, int k){
        return Ez[i][j][k];
    }
}
