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
    
    public EField(int nr,int na,int nz,Cell c){
        super(nr,na,nz,c);
        Er = new double[nr][na+1][nz+1];
        Ea = new double[nr+1][na][nz+1];
        Ez = new double[nr+1][na+1][nz];
        
        Cere = Common.genDouble3DArray(nr,na+1,nz+1,0.0);
        Cerhz = Common.genDouble3DArray(nr,na+1,nz+1,0.0);
        Cerha = Common.genDouble3DArray(nr,na+1,nz+1,0.0);
        Ceae = Common.genDouble3DArray(nr+1,na,nz+1,0.0);
        Ceahr = Common.genDouble3DArray(nr+1,na,nz+1,0.0);
        Ceahz = Common.genDouble3DArray(nr+1,na,nz+1,0.0);
        Ceze = Common.genDouble3DArray(nr+1,na+1,nz,0.0);
        Cezha = Common.genDouble3DArray(nr+1,na+1,nz,0.0);
        Cezhr = Common.genDouble3DArray(nr+1,na+1,nz,0.0);
    }
    
    public double[][] getEFieldArray(int d){
        return Ez[d];
    }
    
    public void updatingCoefficients(ProblemSpace ps){
        System.out.println("General electric field updating coefficients...");
        MaterialGrid mg = ps.getMaterialGrid();
	for(int i=0;i<nr;i++){
            for(int j=0;j<na+1;j++){
		for(int k=0;k<nz+1;k++){
                    // Coeffiecients updating Er
                    Cere[i][j][k]  =  (2*mg.getEpsRR(i, j, k)*Constants.EPS0 - dt*mg.getSigmaER(i, j, k))/(2*mg.getEpsRR(i, j, k)*Constants.EPS0 + dt*mg.getSigmaER(i, j, k));
                    Cerhz[i][j][k]=  (2*dt/(c.getR()*c.getDeltaA()))/(2*mg.getEpsRR(i, j, k)*Constants.EPS0 + dt*mg.getSigmaER(i, j, k));
                    Cerha[i][j][k]= -(2*dt/c.getDeltaZ())/(2*mg.getEpsRR(i, j, k)*Constants.EPS0 + dt*mg.getSigmaER(i, j, k));
		}
            }
	}
	
	for(int i=0;i<nr+1;i++){
            for(int j=0;j<na;j++){
                for(int k=0;k<nz+1;k++){
                    //Coeffiecients updating Ea
                    Ceae[i][j][k]  =  (2*mg.getEpsRA(i, j, k)*Constants.EPS0 - dt*mg.getSigmaEA(i, j, k))/(2*mg.getEpsRA(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEA(i, j, k));
                    Ceahr[i][j][k] =  (2*dt/c.getDeltaZ())/(2*mg.getEpsRA(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEA(i, j, k));
                    Ceahz[i][j][k] = -(2*dt/c.getDeltaR())/(2*mg.getEpsRA(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEA(i, j, k));
		}
            }
	}
	
	for(int i=0;i<nr+1;i++){
            for(int j=0;j<na+1;j++){
		for(int k=0;k<nz;k++){
                    //Coeffiecients updating Ez
                    Ceze[i][j][k]  =  (2*mg.getEpsRZ(i, j, k)*Constants.EPS0 - dt*mg.getSigmaEZ(i, j, k))/(2*mg.getEpsRZ(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEZ(i, j, k));
                    Cezha[i][j][k] =  (2*dt/(c.getR()*c.getDeltaR()))/(2*mg.getEpsRZ(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEZ(i, j, k));
                    Cezhr[i][j][k] = -(2*dt/(c.getR()*c.getDeltaA()))/(2*mg.getEpsRZ(i, j, k)*Constants.EPS0 + dt*mg.getSigmaEZ(i, j, k));
		}
            }
	}
    }
    
    public void setAllR(double val){
        for(int i=0; i<nr; i++){
            for(int j=0; j<na+1; j++){
                for(int k=0; k<nz+1; k++){
                    Er[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllA(double val){
        for(int i=0; i<nr+1; i++){
            for(int j=0; j<na; j++){
                for(int k=0; k<nz+1; k++){
                    Ea[i][j][k] = val;
                }
            }
        }
    }
    
    public void setAllZ(double val){
        for(int i=0; i<nr+1; i++){
            for(int j=0; j<na+1; j++){
                for(int k=0; k<nz; k++){
                    Ez[i][j][k] = val;
                }
            }
        }
    }
    
    public void updateEField(){
        double R = c.getR();
        double RDeltaR = c.getR()-c.getDeltaR();
	for (int i=0;i<nr;i++){
	    for (int j=1;j<na;j++){
		for (int k=1;k<nz;k++){
		    Er[i][j][k] = Cere[i][j][k]*Er[i][j][k] + Cerhz[i][j][k]*(Hz[i][j][k] - Hz[i][j-1][k]) + Cerha[i][j][k]*(Ha[i][j][k] - Ha[i][j][k-1]);   
		}
	    }
	}
	for (int i=1;i<nr;i++){
	    for (int j=0;j<na;j++){
		for (int k=1;k<nz;k++){
		    Ea[i][j][k] = Ceae[i][j][k]*Ea[i][j][k] + Ceahr[i][j][k]*(Hr[i][j][k] - Hr[i][j][k-1]) + Ceahz[i][j][k]*(Hz[i][j][k] - Hz[i-1][j][k]); 
		}
	    }
	}
	for (int i=1;i<nr;i++){
	    for (int j=1;j<na;j++){
		for (int k=0;k<nz;k++){
		    Ez[i][j][k] = Ceze[i][j][k]*Ez[i][j][k] + Cezha[i][j][k]*(R*Ha[i][j][k] - RDeltaR*Ha[i-1][j][k]) + Cezhr[i][j][k]*(Hr[i][j][k] - Hr[i][j-1][k]);
		}
	    }
	}	
    }
    
    public double getEr(int i, int j, int k){
        return Er[i][j][k];
    }
    
    public double getEa(int i, int j, int k){
        return Ea[i][j][k];
    }
        
    public double getEz(int i, int j, int k){
        return Ez[i][j][k];
    }
}
