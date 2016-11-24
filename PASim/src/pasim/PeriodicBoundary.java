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
public class PeriodicBoundary extends EMobject{
    public PeriodicBoundary(){
        super();
    }
    
    public void updatePBCinY(){
	//Updating E_x at y=0 and y=Py
        double HzAt0;
        for (int i=0;i<nx;i++){
            for (int k=1;k<nz;k++){
                HzAt0 = Hz[i][ny-1][k];
                Ex[i][0][k] = Cexe[i][0][k]*Ex[i][0][k] + Cexhz[i][0][k]*(Hz[i][0][k] - HzAt0) + Cexhy[i][0][k]*(Hy[i][0][k] - Hy[i][0][k-1]);
                Ex[i][ny][k] = Ex[i][0][k];
            }
        }
        
        //Updating Ez at y=0, y=Py
        double HxAt0;
        for (int i=1;i<nx;i++){
            for (int k=0;k<nz;k++){
                HxAt0 = Hx[i][ny-1][k];
                Ez[i][0][k] = Ceze[i][0][k]*Ez[i][0][k] + Cezhy[i][0][k]*(Hy[i][0][k] - Hy[i-1][0][k]) + Cezhx[i][0][k]*(Hx[i][0][k] - HxAt0);
                Ez[i][ny][k] = Ez[i][0][k];
            }
        }
    }
    
    public void updatePBCinX(){
        //Updating Ey at x=0 and x=Px
        for (int j=0;j<ny;j++){
            for (int k=1;k<nz;k++){
                Ey[0][j][k] = Ceye[0][j][k]*Ey[0][j][k] + Ceyhx[0][j][k]*(Hx[0][j][k] - Hx[0][j][k-1]) + Ceyhz[0][j][k]*(Hz[0][j][k] - Hz[nx-1][j][k]);
                Ey[nx][j][k] = Ey[0][j][k];
	    }
	}
        
	//Updating Ez at y=0, y=Py, x=0, and x=Px
        for (int j=1;j<ny;j++){
            for (int k=0;k<nz;k++){
                Ez[0][j][k] = Ceze[0][j][k]*Ez[0][j][k] + Cezhy[0][j][k]*(Hy[0][j][k] - Hy[nx-1][j][k]) + Cezhx[0][j][k]*(Hx[0][j][k] - Hx[0][j-1][k]);
                Ez[nx][j][k] = Ez[0][j][k];
            }
        }
    }
    
    public void updatePBCinXY(){
        //Updating E_x at y=0 and y=Py
        double HzAt0;
        for (int i=0;i<nx;i++){
            for (int k=1;k<nz;k++){
                HzAt0 = Hz[i][ny-1][k];
                Ex[i][0][k] = Cexe[i][0][k]*Ex[i][0][k] + Cexhz[i][0][k]*(Hz[i][0][k] - HzAt0) + Cexhy[i][0][k]*(Hy[i][0][k] - Hy[i][0][k-1]);
                Ex[i][ny][k] = Ex[i][0][k];
            }
        }
        
        //Updating Ey at x=0 and x=Px
        for (int j=0;j<ny;j++){
            for (int k=1;k<nz;k++){
                Ey[0][j][k] = Ceye[0][j][k]*Ey[0][j][k] + Ceyhx[0][j][k]*(Hx[0][j][k] - Hx[0][j][k-1]) + Ceyhz[0][j][k]*(Hz[0][j][k] - Hz[nx-1][j][k]);
                Ey[nx][j][k] = Ey[0][j][k];
	    }
	}
        
        //Updating Ez at y=0, y=Py, x=0, and x=Px
        double HxAt0;
        for (int i=1;i<nx;i++){
            for (int k=0;k<nz;k++){
                HxAt0 = Hx[i][ny-1][k];
                Ez[i][0][k] = Ceze[i][0][k]*Ez[i][0][k] + Cezhy[i][0][k]*(Hy[i][0][k] - Hy[i-1][0][k]) + Cezhx[i][0][k]*(Hx[i][0][k] - HxAt0);
                Ez[i][ny][k] = Ez[i][0][k];
            }
        }
        
        for (int j=1;j<ny;j++){
            for (int k=0;k<nz;k++){
                Ez[0][j][k] = Ceze[0][j][k]*Ez[0][j][k] + Cezhy[0][j][k]*(Hy[0][j][k] - Hy[nx-1][j][k]) + Cezhx[0][j][k]*(Hx[0][j][k] - Hx[0][j-1][k]);
                Ez[nx][j][k] = Ez[0][j][k];
            }
        }
        
        // Update Ez at the cornaer
        for (int k=0;k<nz;k++){
            Ez[0][0][k] = Ceze[0][0][k]*Ez[0][0][k] + Cezhy[0][0][k]*(Hy[0][0][k] - Hy[nx-1][0][k]) + Cezhx[0][0][k]*(Hx[0][0][k] - Hx[0][ny-1][k]);
            //Ez[nx][0][k] = Ez[0][0][k];
            //Ez[0][ny][k] = Ez[nx][0][k];
            //Ez[nx][ny][k] = Ez[0][ny][k];
            Ez[nx][0][k] = Ceze[nx][0][k]*Ez[nx][0][k] + Cezhy[nx][0][k]*(Hy[0][0][k] - Hy[nx-1][0][k]) + Cezhx[nx][0][k]*(Hx[nx][0][k] - Hx[nx][ny-1][k]);
            Ez[0][ny][k] = Ceze[0][ny][k]*Ez[0][ny][k] + Cezhy[0][ny][k]*(Hy[0][ny][k] - Hy[nx-1][ny][k]) + Cezhx[0][0][k]*(Hx[0][0][k] - Hx[0][ny-1][k]);
            Ez[nx][ny][k] = Ceze[nx][ny][k]*Ez[nx][ny][k] + Cezhy[nx][ny][k]*(Hy[0][0][k] - Hy[nx-1][ny][k]) + Cezhx[nx][ny][k]*(Hx[0][0][k] - Hx[nx][ny-1][k]);
        }
    }
    
    /*
    public void updatePBCinXY(){
        int gap = 1;
        //Updating E_x at y=0 and y=Py
        double HzAt0;
        for (int i=0;i<nx;i++){
            for (int k=gap;k<nz+1-gap;k++){
                HzAt0 = Hz[i][ny-1][k];
                Ex[i][0][k] = Cexe[i][0][k]*Ex[i][0][k] + Cexhz[i][0][k]*(Hz[i][0][k] - HzAt0) + Cexhy[i][0][k]*(Hy[i][0][k] - Hy[i][0][k-1]);
                Ex[i][ny][k] = Ex[i][0][k];
            }
        }
        
        //Updating Ey at x=0 and x=Px
        for (int j=0;j<ny;j++){
            for (int k=gap;k<nz+1-gap;k++){
                Ey[0][j][k] = Ceye[0][j][k]*Ey[0][j][k] + Ceyhx[0][j][k]*(Hx[0][j][k] - Hx[0][j][k-1]) + Ceyhz[0][j][k]*(Hz[0][j][k] - Hz[nx-1][j][k]);
                Ey[nx][j][k] = Ey[0][j][k];
	    }
	}
        
        //Updating Ez at y=0, y=Py, x=0, and x=Px
        double HxAt0;
        for (int i=1;i<nx;i++){
            for (int k=gap;k<nz-gap;k++){
                HxAt0 = Hx[i][ny-1][k];
                Ez[i][0][k] = Ceze[i][0][k]*Ez[i][0][k] + Cezhy[i][0][k]*(Hy[i][0][k] - Hy[i-1][0][k]) + Cezhx[i][0][k]*(Hx[i][0][k] - HxAt0);
                Ez[i][ny][k] = Ez[i][0][k];
            }
        }
        
        for (int j=1;j<ny;j++){
            for (int k=gap;k<nz-gap;k++){
                Ez[0][j][k] = Ceze[0][j][k]*Ez[0][j][k] + Cezhy[0][j][k]*(Hy[0][j][k] - Hy[nx-1][j][k]) + Cezhx[0][j][k]*(Hx[0][j][k] - Hx[0][j-1][k]);
                Ez[nx][j][k] = Ez[0][j][k];
            }
        }
        
        // Update Ez at the cornaer
        for (int k=gap;k<nz-gap;k++){
            Ez[0][0][k] = Ceze[0][0][k]*Ez[0][0][k] + Cezhy[0][0][k]*(Hy[0][0][k] - Hy[nx-1][0][k]) + Cezhx[0][0][k]*(Hx[0][0][k] - Hx[0][ny-1][k]);
            Ez[nx][0][k] = Ez[0][0][k];
            Ez[0][ny][k] = Ez[nx][0][k];
            Ez[nx][ny][k] = Ez[0][ny][k];
        }
    }
    */
}
