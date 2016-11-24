package pasim;

public class PeriodicBoundary{
    public int nx,ny,nz;
    
    public PeriodicBoundary(ProblemSpace ps){
        this.nx = ps.getNX();
        this.ny = ps.getNY();
        this.nz = ps.getNZ();
    }
    
    public void updatePBCinY(EField E,HField H){
	//Updating E_x at y=0 and y=Py
        double HzAt0;
        for (int i=0;i<nx;i++){
            for (int k=1;k<nz;k++){
                HzAt0 = H.getHZ(i,ny-1,k);
                E.setEX(i, 0, k, E.getCexe(i,0,k)*E.getEX(i,0,k) + E.getCexhz(i,0,k)*(H.getHZ(i,0,k) - HzAt0) + E.getCexhy(i,0,k)*(H.getHY(i,0,k) - H.getHY(i,0,k-1)));
                E.setEX(i,ny,k, E.getEX(i,0,k));
            }
        }
        
        //Updating Ez at y=0, y=Py
        double HxAt0;
        for (int i=1;i<nx;i++){
            for (int k=0;k<nz;k++){
                HxAt0 = H.getHX(i,ny-1,k);
                E.setEZ(i,0,k, E.getCeze(i,0,k)*E.getEZ(i,0,k) + E.getCezhy(i,0,k)*(H.getHY(i,0,k) - H.getHY(i-1,0,k)) + E.getCezhx(i,0,k)*(H.getHX(i,0,k) - HxAt0));
                E.setEZ(i, ny, k, E.getEZ(i,0,k));
            }
        }
    }
    
    public void updatePBCinX(EField E,HField H){
        //Updating Ey at x=0 and x=Px
        for (int j=0;j<ny;j++){
            for (int k=1;k<nz;k++){
                E.setEY(0, j, k, E.getCeye(0,j,k)*E.getEY(0,j,k) + E.getCeyhx(0,j,k)*(H.getHX(0,j,k) - H.getHX(0,j,k-1)) + E.getCeyhz(0,j,k)*(H.getHZ(0,j,k) - H.getHZ(nx-1,j,k)));
                E.setEY(nx,j,k, E.getEY(0,j,k));
	    }
	}
        
	//Updating Ez at y=0, y=Py, x=0, and x=Px
        for (int j=1;j<ny;j++){
            for (int k=0;k<nz;k++){
                E.setEZ(0, j, k, E.getCeze(0,j,k)*E.getEZ(0,j,k) + E.getCezhy(0,j,k)*(H.getHY(0,j,k) - H.getHY(nx-1,j,k)) + E.getCezhx(0,j,k)*(H.getHX(0,j,k) - H.getHX(0,j-1,k)));
                E.setEZ(nx, j, k, E.getEZ(0,j,k));
            }
        }
    }
    
    public void updatePBCinXY(EField E,HField H){
        //Updating E_x at y=0 and y=Py
        double HzAt0;
        for (int i=0;i<nx;i++){
            for (int k=1;k<nz;k++){
                HzAt0 = H.getHZ(i,ny-1,k);
                E.setEX(i, 0, k, E.getCexe(i,0,k)*E.getEX(i,0,k) + E.getCexhz(i,0,k)*(H.getHZ(i,0,k) - HzAt0) + E.getCexhy(i,0,k)*(H.getHY(i,0,k) - H.getHY(i,0,k-1)));
                E.setEX(i,ny,k, E.getEX(i,0,k));
            }
        }
        
        //Updating Ey at x=0 and x=Px
        for (int j=0;j<ny;j++){
            for (int k=1;k<nz;k++){
                E.setEY(0,j,k, E.getCeye(0,j,k)*E.getEY(0,j,k) + E.getCeyhx(0,j,k)*(H.getHX(0,j,k) - H.getHX(0,j,k-1)) + E.getCeyhz(0,j,k)*(H.getHZ(0,j,k) - H.getHZ(nx-1,j,k)));
                E.setEY(nx,j,k, E.getEY(0,j,k));
	    }
	}
        
        //Updating Ez at y=0, y=Py, x=0, and x=Px
        double HxAt0;
        for (int i=1;i<nx;i++){
            for (int k=0;k<nz;k++){
                HxAt0 = H.getHX(i,ny-1,k);
                E.setEZ(i,0,k, E.getCeze(i,0,k)*E.getEZ(i,0,k) + E.getCezhy(i,0,k)*(H.getHY(i,0,k) - H.getHY(i-1,0,k)) + E.getCezhx(i,0,k)*(H.getHX(i,0,k) - HxAt0));
                E.setEZ(i,ny,k, E.getEZ(i,0,k));
            }
        }
        
        for (int j=1;j<ny;j++){
            for (int k=0;k<nz;k++){
                E.setEZ(0,j,k, E.getCeze(0,j,k)*E.getEZ(0,j,k) + E.getCezhy(0,j,k)*(H.getHY(0,j,k) - H.getHY(nx-1,j,k)) + E.getCezhx(0,j,k)*(H.getHX(0,j,k) - H.getHX(0,j-1,k)));
                E.setEZ(nx,j,k, E.getEZ(0,j,k));
            }
        }
        
        // Update Ez at the cornaer
        for (int k=0;k<nz;k++){
            E.setEZ(0,0,k, E.getCeze(0,0,k)*E.getEZ(0,0,k) + E.getCezhy(0,0,k)*(H.getHY(0,0,k) - H.getHY(nx-1,0,k)) + E.getCezhx(0,0,k)*(H.getHX(0,0,k) - H.getHX(0,ny-1,k)));
            //Ez[nx][0][k] = Ez[0][0][k];
            //Ez[0][ny][k] = Ez[nx][0][k];
            //Ez[nx][ny][k] = Ez[0][ny][k];
            E.setEZ(nx,0,k, E.getCeze(nx,0,k)*E.getEZ(nx,0,k) + E.getCezhy(nx,0,k)*(H.getHY(0,0,k) - H.getHY(nx-1,0,k)) + E.getCezhx(nx,0,k)*(H.getHX(nx,0,k) - H.getHX(nx,ny-1,k)));
            E.setEZ(0,ny,k, E.getCeze(0,ny,k)*E.getEZ(0,ny,k) + E.getCezhy(0,ny,k)*(H.getHY(0,ny,k) - H.getHY(nx-1,ny,k)) + E.getCezhx(0,0,k)*(H.getHX(0,0,k) - H.getHX(0,ny-1,k)));
            E.setEZ(nx,ny,k, E.getCeze(nx,ny,k)*E.getEZ(nx,ny,k) + E.getCezhy(nx,ny,k)*(H.getHY(0,0,k) - H.getHY(nx-1,ny,k)) + E.getCezhx(nx,ny,k)*(H.getHX(0,0,k) - H.getHX(nx,ny-1,k)));
        }
    }
}
