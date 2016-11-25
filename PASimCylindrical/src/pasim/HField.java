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
    
    public HField(int nr,int na,int nz,Cell c){
        super(nr,na,nz,c);
        Hr = new double[nr+1][na][nz];
        Ha = new double[nr][na+1][nz];
        Hz = new double[nr][na][nz+1];
	
        Chrh = Common.genDouble3DArray(nr+1, na, nz, 0.0);
        Chrez = Common.genDouble3DArray(nr+1, na, nz, 0.0);
        Chrea = Common.genDouble3DArray(nr+1, na, nz, 0.0);
        Chah = Common.genDouble3DArray(nr, na+1, nz, 0.0);
        Chaer = Common.genDouble3DArray(nr, na+1, nz, 0.0);
        Chaez = Common.genDouble3DArray(nr, na+1, nz, 0.0);
        Chzh = Common.genDouble3DArray(nr, na, nz+1, 0.0);
        Chzea = Common.genDouble3DArray(nr, na, nz+1, 0.0);
        Chzer = Common.genDouble3DArray(nr, na, nz+1, 0.0);
    }
    
    public void updatingCoefficients(ProblemSpace ps){
        System.out.println("General magnetic field updating coefficients....");
        MaterialGrid mg = ps.getMaterialGrid();
 	for(int i=0;i<nr+1;i++){
	    for(int j=0;j<na;j++){
		for(int k=0;k<nz;k++){
		    //Coeffiecients updating Hr
		    Chrh[i][j][k]  =  (2*mg.getMuRR(i, j, k)*Constants.MU0 - dt*mg.getSigmaMR(i, j, k))/(2*mg.getMuRR(i, j, k)*Constants.MU0 + dt*mg.getSigmaMR(i, j, k));
		    Chrez[i][j][k] = -(2*dt/(c.getR()*c.getDeltaA()))/(2*mg.getMuRR(i, j, k)*Constants.MU0 + dt*mg.getSigmaMR(i, j, k));
		    Chrea[i][j][k] =  (2*dt/c.getDeltaZ())/(2*mg.getMuRR(i, j, k)*Constants.MU0 + dt*mg.getSigmaMR(i, j, k));
		}
	    }
	}
	
	for(int i=0;i<nr;i++){
	    for(int j=0;j<na+1;j++){
		for(int k=0;k<nz;k++){
		    //Coeffiecients updating Ha
		    Chah[i][j][k]  =  (2*mg.getMuRA(i, j, k)*Constants.MU0 - dt*mg.getSigmaMA(i, j, k))/(2*mg.getMuRA(i, j, k)*Constants.MU0 + dt*mg.getSigmaMA(i, j, k));
		    Chaer[i][j][k] = -(2*dt/c.getDeltaZ())/(2*mg.getMuRA(i, j, k)*Constants.MU0 + dt*mg.getSigmaMA(i, j, k));
		    Chaez[i][j][k] =  (2*dt/c.getDeltaR())/(2*mg.getMuRA(i, j, k)*Constants.MU0 + dt*mg.getSigmaMA(i, j, k));
		}
	    }
	}
	
	for(int i=0;i<nr;i++){
	    for(int j=0;j<na;j++){
		for(int k=0;k<nz+1;k++){
		    //Coeffiecients updating Hz
		    Chzh[i][j][k]  =  (2*mg.getMuRZ(i, j, k)*Constants.MU0 - dt*mg.getSigmaMZ(i, j, k))/(2*mg.getMuRZ(i, j, k)*Constants.MU0 + dt*mg.getSigmaMZ(i, j, k));
		    Chzea[i][j][k] = -(2*dt/c.getDeltaR())/(2*mg.getMuRZ(i, j, k)*Constants.MU0 + dt*mg.getSigmaMZ(i, j, k));
		    Chzer[i][j][k] =  (2*dt/c.getDeltaA())/(2*mg.getMuRZ(i, j, k)*Constants.MU0 + dt*mg.getSigmaMZ(i, j, k));
		}
	    }
	}       
    }
    
    public double getHr(int i, int j, int k){
        return Hr[i][j][k];
    }
    
    public double getHa(int i, int j, int k){
        return Ha[i][j][k];
    }
        
    public double getHz(int i, int j, int k){
        return Hz[i][j][k];
    }
    
    public void setAllR(double val){
        for(int i=0; i<nr+1; i++){
            for(int j=0; j<na; j++){
                for(int k=0; k<nz; k++){ 
                    Hr[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllA(double val){
        for(int i=0; i<nr; i++){
            for(int j=0; j<na+1; j++){
                for(int k=0; k<nz; k++){ 
                    Ha[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllZ(double val){
        for(int i=0; i<nr; i++){
            for(int j=0; j<na; j++){
                for(int k=0; k<nz+1; k++){
                    Hz[i][j][k] = val;
                }
            }
        }
    }
    
    public void updateHField(){
        double R = c.getR();
        double RdeltaR = c.getR()+c.getDeltaR();
	for (int i=0; i<nr+1; i++){
	    for (int j=0; j<na; j++){
		for (int k=0; k<nz; k++){	   
		    Hr[i][j][k] = Chrh[i][j][k]*Hr[i][j][k] + Chrea[i][j][k]*(Ea[i][j][k+1] - Ea[i][j][k]) + Chrez[i][j][k]*(Ez[i][j+1][k] - Ez[i][j][k]);
		}
	    }
	}
	for (int i=0; i<nr; i++){
	    for (int j=0; j<na+1; j++){
		for (int k=0; k<nz; k++){
		    Ha[i][j][k] = Chah[i][j][k]*Ha[i][j][k] + Chaez[i][j][k]*(Ez[i+1][j][k] - Ez[i][j][k]) + Chaer[i][j][k]*(Er[i][j][k+1] - Er[i][j][k]);
		}
	    }
	}
	for (int i=0; i<nr; i++){
	    for (int j=0; j<na; j++){
		for (int k=0; k<nz+1; k++){
		    Hz[i][j][k] = Chzh[i][j][k]*Hz[i][j][k] + Chzer[i][j][k]*(Er[i][j+1][k] - Er[i][j][k]) + Chzea[i][j][k]*(RdeltaR*Ea[i+1][j][k] - R*Ea[i][j][k]);
		}
	    }
	}
    }
}
