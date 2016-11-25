/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pasim;

/**
 *
 * @author Sudantha
 */
public class FarField extends EMobject{
    private double[] frequencies, omegas;;
    private int outerBoundaryCellCount, angleStepsCount;
    private double[] radiatedPower;
    private int li,lj,lk,ui,uj,uk;
    private double[] theta,phi;
    private double[][] farfield_dataTheta;
    private double[][] farfield_dataPhi;
    
    private double[][][][] tjrap, tjrzp, tjarp, tjazp, tjzrp, tjzap, tjran, tjrzn, tjarn, tjazn, tjzrn, tjzan, tmrap, tmrzp, tmarp, tmazp, tmzrp, tmzap, tmran, tmrzn, tmarn, tmazn, tmzrn, tmzan;
    private Complex[][][][] cjrap, cjrzp, cjarp, cjazp, cjzrp, cjzap, cjran, cjrzn, cjarn, cjazn, cjzrn, cjzan, cmrap, cmrzp, cmarp, cmazp, cmzrp, cmzap, cmran, cmrzn, cmarn, cmazn, cmzrn, cmzan;
    
    public FarField(double start,double stop,int count){
        double step = (stop - start)/(count - 1);
        frequencies = new double[count];
        
        for (int i = 0; i < count; i++){
            frequencies[i] = start + i*step;
        }
        
        angleStepsCount = 360;
        theta = Common.genDouble1DArray(angleStepsCount, 0.0);
        phi = Common.genDouble1DArray(angleStepsCount, 0.0);
        
        for(int i=0;i<angleStepsCount;i++){
            theta[i]=i*Math.PI/angleStepsCount;
            phi[i]=i*(2*Math.PI/angleStepsCount)-2*Math.PI;
        }
    }
    
    public void setupFarFieldOutPut(double[] theta, double[] phi){
        this.angleStepsCount = theta.length;
        this.theta = theta;
        this.phi = phi;
    }
    
    public void setupTheta(double start,double end,double stepsCount){
        start = (Math.PI/180)*start;
        end = (Math.PI/180)*end;
        double stepSize = (end-start)/(stepsCount-1);
        for(int i=0;i<angleStepsCount;i++){
            theta[i] = start + stepSize*i;
        }
    }
    
    public void setupPhi(double start,double end,double stepsCount){
        start = (Math.PI/180)*start;
        end = (Math.PI/180)*end;
        double stepSize = (Math.PI/180)*(end-start)/stepsCount;
        for(int i=0;i<angleStepsCount;i++){
            phi[i] = start + stepSize*i;
        }
    }
    
    public void saveFarFieldThetaPhiPattern(){
        Common.save2DArray(farfield_dataPhi, "FarFieldPhiPattern");
        Common.save2DArray(farfield_dataTheta, "FarFieldThetaPattern");
        Common.save1DArray(phi, "FarFIeldPhi");
        Common.save1DArray(theta, "FarFieldTheta");
    }
    
    public void calFarField(){ 
        Complex[] exp_jk_rpr = Common.genComplex1DArray(angleStepsCount, new Complex());
        double[] dr_sinth_cosphi = Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] da_sinth_sinphi = Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dz_costh = Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] da_dz_costh_sinphi = Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] da_dz_sinth = Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] da_dz_cosphi = Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dr_dz_costh_cosphi = Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dr_dz_sinth = Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dr_dz_sinphi = Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dr_da_costh_cosphi = Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dr_da_costh_sinphi = Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dr_da_sinphi = Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dr_da_cosphi = Common.genDouble1DArray(angleStepsCount, 0.0);
        
        //double[][] farfield_dirTheta = Common.genDouble2DArray(frequencies.length,angleStepsCount,0.0);
        //double[][] farfield_dir = Common.genDouble2DArray(frequencies.length,angleStepsCount,0.0);
        //double[][] farfield_dirPhi = Common.genDouble2DArray(frequencies.length,angleStepsCount,0.0);
        
        farfield_dataTheta = Common.genDouble2DArray(frequencies.length,angleStepsCount,0.0);
        farfield_dataPhi = Common.genDouble2DArray(frequencies.length,angleStepsCount,0.0);
                
        for (int i = 0;i<angleStepsCount;i++){
            dr_sinth_cosphi[i] = c.getDeltaR()*Math.sin(theta[i])*Math.cos(phi[i]);
            da_sinth_sinphi[i] = c.getDeltaA()*Math.sin(theta[i])*Math.sin(phi[i]);
            dz_costh[i] = c.getDeltaZ()*Math.cos(theta[i]);
            da_dz_costh_sinphi[i] = c.getDeltaA()*c.getDeltaZ()*Math.cos(theta[i])*Math.sin(phi[i]);
            da_dz_sinth[i] = c.getDeltaA()*c.getDeltaZ()*Math.sin(theta[i]);
            da_dz_cosphi[i] = c.getDeltaA()*c.getDeltaZ()*Math.cos(phi[i]);
            dr_dz_costh_cosphi[i] = c.getDeltaR()*c.getDeltaZ()*Math.cos(theta[i])*Math.cos(phi[i]);
            dr_dz_sinth[i] = c.getDeltaR()*c.getDeltaZ()*Math.sin(theta[i]);
            dr_dz_sinphi[i] = c.getDeltaR()*c.getDeltaZ()*Math.sin(phi[i]);
            dr_da_costh_cosphi[i] = c.getDeltaR()*c.getDeltaA()*Math.cos(theta[i])*Math.cos(phi[i]);
            dr_da_costh_sinphi[i] = c.getDeltaR()*c.getDeltaA()*Math.cos(theta[i])*Math.sin(phi[i]);
            dr_da_sinphi[i] = c.getDeltaR()*c.getDeltaA()*Math.sin(phi[i]);
            dr_da_cosphi[i] = c.getDeltaR()*c.getDeltaA()*Math.cos(phi[i]);
        }
        
        double ci = 0.5*(ui+li);
        double cj = 0.5*(uj+lj);
        double ck = 0.5*(uk+lk);
        
        for (int mi=0; mi<this.frequencies.length; mi++){
                    
            Complex[] Ntheta = Common.genComplex1DArray(angleStepsCount, new Complex());
            Complex[] Ltheta = Common.genComplex1DArray(angleStepsCount, new Complex());
            Complex[] Nphi = Common.genComplex1DArray(angleStepsCount, new Complex());
            Complex[] Lphi = Common.genComplex1DArray(angleStepsCount, new Complex());
            double[] rpr = Common.genDouble1DArray(angleStepsCount,1);
        
            double k;
            
            System.out.println("Calculating directivity for ".concat(String.valueOf(this.frequencies[mi]).concat(" Hz")));
            k = 2*Math.PI*frequencies[mi]*Math.sqrt(Constants.MU0*Constants.EPS0);
            for (int nj=lj;nj<uj;nj++){
                for (int nk=lk;nk<uk;nk++){
                    for (int i=0;i<angleStepsCount;i++){
                        //for +ax direction
                        rpr[i] = (ui - ci)*dr_sinth_cosphi[i] + (nj-cj+0.5)*da_sinth_sinphi[i] + (nk-ck+0.5)*dz_costh[i];
                        exp_jk_rpr[i].set(Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i])));
                        Ntheta[i].set(Ntheta[i].plus((cjarp[0][nj-lj][nk-lk][mi].times(da_dz_costh_sinphi[i]).minus(cjzrp[0][nj-lj][nk-lk][mi].times(da_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Ltheta[i].set(Ltheta[i].plus((cmarp[0][nj-lj][nk-lk][mi].times(da_dz_costh_sinphi[i]).minus(cmzrp[0][nj-lj][nk-lk][mi].times(da_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Nphi[i].set(Nphi[i].plus((cjarp[0][nj-lj][nk-lk][mi].times(da_dz_cosphi[i])).times(exp_jk_rpr[i])));
                        Lphi[i].set(Lphi[i].plus((cmarp[0][nj-lj][nk-lk][mi].times(da_dz_cosphi[i])).times(exp_jk_rpr[i])));
                        
                        //for -ax direction
                        rpr[i] = (li - ci)*dr_sinth_cosphi[i] + (nj-cj+0.5)*da_sinth_sinphi[i] + (nk-ck+0.5)*dz_costh[i];
                        exp_jk_rpr[i].set(Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i])));
                        Ntheta[i].set(Ntheta[i].plus((cjarn[0][nj-lj][nk-lk][mi].times(da_dz_costh_sinphi[i]).minus(cjzrn[0][nj-lj][nk-lk][mi].times(da_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Ltheta[i].set(Ltheta[i].plus(cmarn[0][nj-lj][nk-lk][mi].times(da_dz_costh_sinphi[i]).minus(cmzrn[0][nj-lj][nk-lk][mi].times(da_dz_sinth[i])).times(exp_jk_rpr[i])));
                        Nphi[i].set(Nphi[i].plus(cjarn[0][nj-lj][nk-lk][mi].times(da_dz_cosphi[i]).times(exp_jk_rpr[i])));
                        Lphi[i].set(Lphi[i].plus(cmarn[0][nj-lj][nk-lk][mi].times(da_dz_cosphi[i]).times(exp_jk_rpr[i])));
                    }
                }
            }
           
            for (int ni=li;ni<ui;ni++){
                for (int nk=lk;nk<uk;nk++){
                    for (int i=0;i<angleStepsCount;i++){
                        //for +aa direction
                        rpr[i] = (ni - ci + 0.5)*dr_sinth_cosphi[i] + (uj-cj)*da_sinth_sinphi[i] + (nk-ck+0.5)*dz_costh[i];
                        exp_jk_rpr[i].set(Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i])));
                        Ntheta[i].set( Ntheta[i].plus((cjrap[ni-li][0][nk-lk][mi].times(dr_dz_costh_cosphi[i]).minus(cjzap[ni-li][0][nk-lk][mi].times(dr_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Ltheta[i].set(Ltheta[i].plus((cmrap[ni-li][0][nk-lk][mi].times(dr_dz_costh_cosphi[i]).minus(cmzap[ni-li][0][nk-lk][mi].times(dr_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Nphi[i].set(Nphi[i].plus(((cjrap[ni-li][0][nk-lk][mi].times(-1.0)).times(dr_dz_sinphi[i])).times(exp_jk_rpr[i])));
                        Lphi[i].set(Lphi[i].plus(((cmrap[ni-li][0][nk-lk][mi].times(-1.0)).times(dr_dz_sinphi[i])).times(exp_jk_rpr[i])));

                        //for -aa direction
                        rpr[i] = (ni - ci + 0.5)*dr_sinth_cosphi[i] + (lj-cj)*da_sinth_sinphi[i] + (nk-ck+0.5)*dz_costh[i];
                        exp_jk_rpr[i].set(Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i])));
                        Ntheta[i].set(Ntheta[i].plus(((cjran[ni-li][0][nk-lk][mi].times(dr_dz_costh_cosphi[i])).minus(cjzan[ni-li][0][nk-lk][mi].times(dr_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Ltheta[i].set(Ltheta[i].plus(((cmran[ni-li][0][nk-lk][mi].times(dr_dz_costh_cosphi[i])).minus(cmzan[ni-li][0][nk-lk][mi].times(dr_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Nphi[i].set(Nphi[i].plus(((cjran[ni-li][0][nk-lk][mi].times(-1.0)).times(dr_dz_sinphi[i])).times(exp_jk_rpr[i])));
                        Lphi[i].set(Lphi[i].plus(((cmran[ni-li][0][nk-lk][mi].times(-1.0)).times(dr_dz_sinphi[i])).times(exp_jk_rpr[i])));
                    }
                }
            }
            
            for (int ni=li;ni<ui;ni++){
                for (int nj=lj;nj<uj;nj++){
                    for (int i=0;i<angleStepsCount;i++){
                        //for +az direction

                        rpr[i] = (ni-ci+0.5)*dr_sinth_cosphi[i] + (nj - cj + 0.5)*da_sinth_sinphi[i] + (uk-ck)*dz_costh[i];
                        exp_jk_rpr[i].set(Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i])));
                        Ntheta[i].set(Ntheta[i].plus(((cjrzp[ni-li][nj-lj][0][mi].times(dr_da_costh_cosphi[i])).plus(cjazp[ni-li][nj-lj][0][mi].times(dr_da_costh_sinphi[i]))).times(exp_jk_rpr[i])));
                        Ltheta[i].set(Ltheta[i].plus((((cmrzp[ni-li][nj-lj][0][mi].times(dr_da_costh_cosphi[i])).plus((cmazp[ni-li][nj-lj][0][mi].times(dr_da_costh_sinphi[i])))).times(exp_jk_rpr[i]))));
                        Nphi[i].set(Nphi[i].plus((((cjrzp[ni-li][nj-lj][0][mi].times(-1.0)).times(dr_da_sinphi[i])).plus((cjazp[ni-li][nj-lj][0][mi].times(dr_da_cosphi[i])))).times(exp_jk_rpr[i])));
                        Lphi[i].set(Lphi[i].plus((((cmrzp[ni-li][nj-lj][0][mi].times(-1.0)).times(dr_da_sinphi[i])).plus(cmazp[ni-li][nj-lj][0][mi].times(dr_da_cosphi[i]))).times(exp_jk_rpr[i])));

                        //for -az direction
                        rpr[i] = (ni-ci+0.5)*dr_sinth_cosphi[i] + (nj - cj + 0.5)*da_sinth_sinphi[i] + (lk-ck)*dz_costh[i];
                        exp_jk_rpr[i].set(Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i])));
                        Ntheta[i].set(Ntheta[i].plus(((cjrzn[ni-li][nj-lj][0][mi].times(dr_da_costh_cosphi[i])).plus(cjazn[ni-li][nj-lj][0][mi].times(dr_da_costh_sinphi[i]))).times(exp_jk_rpr[i])));
                        Ltheta[i].set(Ltheta[i].plus(((cmrzn[ni-li][nj-lj][0][mi].times(dr_da_costh_cosphi[i])).plus(cmazn[ni-li][nj-lj][0][mi].times(dr_da_costh_sinphi[i]))).times(exp_jk_rpr[i])));
                        Nphi[i].set(Nphi[i].plus((((cjrzn[ni-li][nj-lj][0][mi].times(-1.0)).times(dr_da_sinphi[i])).plus(cjazn[ni-li][nj-lj][0][mi].times(dr_da_cosphi[i]))).times(exp_jk_rpr[i])));
                        Lphi[i].set(Lphi[i].plus((((cmrzn[ni-li][nj-lj][0][mi].times(-1.0)).times(dr_da_sinphi[i])).plus(cmazn[ni-li][nj-lj][0][mi].times(dr_da_cosphi[i]))).times(exp_jk_rpr[i])));                        
                    }
                }
            }
            
            for (int i=0;i<angleStepsCount;i++){
                farfield_dataTheta[mi][i] = (Math.pow(k, 2)/(8*Math.PI*Constants.ETA0*radiatedPower[mi]))*Math.pow((Lphi[i].plus(Ntheta[i].times(Constants.ETA0))).abs(),2);
                farfield_dataPhi[mi][i] = (Math.pow(k, 2)/(8*Math.PI*Constants.ETA0*radiatedPower[mi]))*Math.pow((Ltheta[i].minus(Nphi[i].times(Constants.ETA0))).abs(),2);
            }
        }
    }
    
    
    public void saveTotalRadiatedPower(){
        Common.save1DArray(radiatedPower, "RadiatedPower");
    }
    
    public void calTotalRadiatedPower(){
	radiatedPower = Common.genDouble1DArray(frequencies.length,0.0);

        Complex[] sum_of_crn = Common.genComplex1DArray(frequencies.length,new Complex());
        Complex[] sum_of_crp = Common.genComplex1DArray(frequencies.length,new Complex());
        Complex[] sum_of_can = Common.genComplex1DArray(frequencies.length,new Complex());
        Complex[] sum_of_cap = Common.genComplex1DArray(frequencies.length,new Complex());
        Complex[] sum_of_czn = Common.genComplex1DArray(frequencies.length,new Complex());
        Complex[] sum_of_czp = Common.genComplex1DArray(frequencies.length,new Complex());
        

	for (int ind=0; ind < frequencies.length; ind++){
            for (int i=0;i<ui-li;i++){
		for (int j=0;j<uj-lj;j++){
                    sum_of_czp[ind].set(sum_of_czp[ind].plus((cmazp[i][j][0][ind].times(cjrzp[i][j][0][ind].conj())).minus(cmrzp[i][j][0][ind].times(cjazp[i][j][0][ind].conj()))));
                    sum_of_czn[ind].set(sum_of_czn[ind].plus((cmazn[i][j][0][ind].times(cjrzn[i][j][0][ind].conj())).minus(cmrzn[i][j][0][ind].times(cjazn[i][j][0][ind].conj()))));
		}
            }

            for (int i=0;i<ui-li;i++){
		for (int k=0;k<uk-lk;k++){
                    sum_of_cap[ind].set(sum_of_cap[ind].plus((cmrap[i][0][k][ind].times(cjzap[i][0][k][ind].conj())).minus(cmzap[i][0][k][ind].times(cjrap[i][0][k][ind].conj()))));
                    sum_of_can[ind].set(sum_of_can[ind].plus((cmran[i][0][k][ind].times(cjzan[i][0][k][ind].conj())).minus(cmzan[i][0][k][ind].times(cjran[i][0][k][ind].conj()))));
		}
            }

            for (int j=0;j<uj-lj;j++){
		for (int k=0;k<uk-lk;k++){
                    sum_of_crp[ind].set(sum_of_crp[ind].plus((cmzrp[0][j][k][ind].times(cjarp[0][j][k][ind].conj())).minus(cmarp[0][j][k][ind].times(cjzrp[0][j][k][ind].conj()))));
                    sum_of_crn[ind].set(sum_of_crn[ind].plus((cmzrn[0][j][k][ind].times(cjarn[0][j][k][ind].conj())).minus(cmarn[0][j][k][ind].times(cjzrn[0][j][k][ind].conj()))));
		}
            }
            
            Complex temp1, temp2, temp3, temp4, temp5, temp6; 
            
            temp1 = sum_of_czp[ind].times(c.getDeltaR()*c.getDeltaA());
            temp2 = sum_of_czn[ind].times(c.getDeltaR()*c.getDeltaA());
            temp3 = sum_of_cap[ind].times(c.getDeltaR()*c.getDeltaZ());
            temp4 = sum_of_can[ind].times(c.getDeltaR()*c.getDeltaZ());
            temp5 = sum_of_crp[ind].times(c.getDeltaA()*c.getDeltaZ());
            temp6 = sum_of_crn[ind].times(c.getDeltaA()*c.getDeltaZ());
            
            Complex temp = new Complex();
            
            temp = temp.plus(temp1);
            temp = temp.minus(temp2);
            temp = temp.plus(temp3);
            temp = temp.minus(temp4);
            temp = temp.plus(temp5);
            temp = temp.minus(temp6);
            
            radiatedPower[ind] = (temp.real())/2;
            //radiatedPower[ind] = 0.5*((temp1.minus(temp2.plus(temp3.minus(temp4.plus(temp5.minus(temp6)))))).real());
	}
    }
    
    public void CalculateJandM(int timeIndex){
        Complex zi = new Complex(0.0, 1.0);
	Complex[] exp_h = Common.genComplex1DArray(frequencies.length, new Complex(0.0,0.0));
	Complex[] exp_e = Common.genComplex1DArray(frequencies.length, new Complex(0.0,0.0));
        
	for (int mi=0; mi<frequencies.length; mi++){
	    exp_h[mi] = (Common.complexExp(new Complex(0,-omegas[mi]*(timeIndex + 0.5)*c.getDeltaT()))).times(c.getDeltaT());
	    exp_e[mi] = (Common.complexExp(new Complex(0,-omegas[mi]*(timeIndex + 1)*c.getDeltaT()))).times(c.getDeltaT());
	}
	
	if (frequencies.length > 0){
	    for(int i = li; i <= ui; i++){
		for (int j = lj; j <= uj; j++){
		    for (int k = lk; k <= uk; k++){
			if (i == ui && j != uj && k != uk){
			    tmarp[0][0][j-lj][k-lk] =  0.5*(Ez[ui][j][k] + Ez[ui][j+1][k]);
			    tmzrp[0][0][j-lj][k-lk] = -0.5*(Ea[ui][j][k] + Ea[ui][j][k+1]);
			    
			    tjarp[0][0][j-lj][k-lk] = -0.25*(Hz[ui][j][k] + Hz[ui][j][k+1] + Hz[ui-1][j][k] + Hz[ui-1][j][k+1]);
			    tjzrp[0][0][j-lj][k-lk] =  0.25*(Ha[ui][j][k] + Ha[ui][j+1][k] + Ha[ui-1][j][k] + Ha[ui-1][j+1][k]);
			    
			    tmarn[0][0][j-lj][k-lk] = -0.5*(Ez[li][j][k] + Ez[li][j+1][k]);
			    tmzrn[0][0][j-lj][k-lk] =  0.5*(Ea[li][j][k] + Ea[li][j][k+1]);
			    
			    tjarn[0][0][j-lj][k-lk] =  0.25*(Hz[li][j][k] + Hz[li][j][k+1] + Hz[li-1][j][k] + Hz[li-1][j][k+1]);
			    tjzrn[0][0][j-lj][k-lk] = -0.25*(Ha[li][j][k] + Ha[li][j+1][k] + Ha[li-1][j][k] + Ha[li-1][j+1][k]);
			    
			    //fourier transform
			    for (int mi=0; mi<frequencies.length; mi++){
				cjarp[0][j-lj][k-lk][mi].set(cjarp[0][j-lj][k-lk][mi].plus(exp_h[mi].times(tjarp[0][0][j-lj][k-lk])));
				cjzrp[0][j-lj][k-lk][mi].set(cjzrp[0][j-lj][k-lk][mi].plus(exp_h[mi].times(tjzrp[0][0][j-lj][k-lk])));
				
				cjarn[0][j-lj][k-lk][mi].set(cjarn[0][j-lj][k-lk][mi].plus(exp_h[mi].times(tjarn[0][0][j-lj][k-lk])));
				cjzrn[0][j-lj][k-lk][mi].set(cjzrn[0][j-lj][k-lk][mi].plus(exp_h[mi].times(tjzrn[0][0][j-lj][k-lk])));
				
				cmarp[0][j-lj][k-lk][mi].set(cmarp[0][j-lj][k-lk][mi].plus(exp_e[mi].times(tmarp[0][0][j-lj][k-lk])));
				cmzrp[0][j-lj][k-lk][mi].set(cmzrp[0][j-lj][k-lk][mi].plus(exp_e[mi].times(tmzrp[0][0][j-lj][k-lk]))); 
				
				cmarn[0][j-lj][k-lk][mi].set(cmarn[0][j-lj][k-lk][mi].plus(exp_e[mi].times(tmarn[0][0][j-lj][k-lk])));
				cmzrn[0][j-lj][k-lk][mi].set(cmzrn[0][j-lj][k-lk][mi].plus(exp_e[mi].times(tmzrn[0][0][j-lj][k-lk]))); 
			    }
			}
			
			if (i != ui && j == uj && k != uk){
			    tmrap[0][i-li][0][k-lk] = -0.5*(Ez[i][uj][k] + Ez[i+1][uj][k]);
			    tmzap[0][i-li][0][k-lk] =  0.5*(Er[i][uj][k] + Er[i][uj][k+1]);
			    
			    tjzap[0][i-li][0][k-lk] = -0.25*(Hr[i][uj][k] + Hr[i+1][uj][k] + Hr[i][uj-1][k] + Hr[i+1][uj-1][k]);
			    tjrap[0][i-li][0][k-lk] =  0.25*(Hz[i][uj][k] + Hz[i][uj][k+1] + Hz[i][uj-1][k] + Hz[i][uj-1][k+1]);
			    
			    tmran[0][i-li][0][k-lk] =  0.5*(Ez[i][lj][k] + Ez[i+1][lj][k]);
			    tmzan[0][i-li][0][k-lk] = -0.5*(Er[i][lj][k] + Er[i][lj][k+1]);
			    
			    tjzan[0][i-li][0][k-lk] =  0.25*(Hr[i][lj][k] + Hr[i+1][lj][k] + Hr[i][lj-1][k] + Hr[i+1][lj-1][k]);
			    tjran[0][i-li][0][k-lk] = -0.25*(Hz[i][lj][k] + Hz[i][lj][k+1] + Hz[i][lj-1][k] + Hz[i][lj-1][k+1]);
			    
			    //fourier transform
			    for (int mi=0; mi<frequencies.length; mi++){
				cjrap[i-li][0][k-lk][mi].set(cjrap[i-li][0][k-lk][mi].plus(exp_h[mi].times(tjrap[0][i-li][0][k-lk])));
				cjzap[i-li][0][k-lk][mi].set(cjzap[i-li][0][k-lk][mi].plus(exp_h[mi].times(tjzap[0][i-li][0][k-lk])));
				
				cjran[i-li][0][k-lk][mi].set(cjran[i-li][0][k-lk][mi].plus(exp_h[mi].times(tjran[0][i-li][0][k-lk])));
				cjzan[i-li][0][k-lk][mi].set(cjzan[i-li][0][k-lk][mi].plus(exp_h[mi].times(tjzan[0][i-li][0][k-lk])));
				
				cmrap[i-li][0][k-lk][mi].set(cmrap[i-li][0][k-lk][mi].plus(exp_e[mi].times(tmrap[0][i-li][0][k-lk]))); 
				cmzap[i-li][0][k-lk][mi].set(cmzap[i-li][0][k-lk][mi].plus(exp_e[mi].times(tmzap[0][i-li][0][k-lk]))); 
				
				cmran[i-li][0][k-lk][mi].set(cmran[i-li][0][k-lk][mi].plus(exp_e[mi].times(tmran[0][i-li][0][k-lk]))); 
				cmzan[i-li][0][k-lk][mi].set(cmzan[i-li][0][k-lk][mi].plus(exp_e[mi].times(tmzan[0][i-li][0][k-lk]))); 
			    }
			}
			
			if (i != ui && j != uj && k == uk){
			    tmrzp[0][i-li][j-lj][0] =  0.5*(Ea[i][j][uk] + Ea[i+1][j][uk]);
			    tmazp[0][i-li][j-lj][0] = -0.5*(Er[i][j][uk] + Er[i][j+1][uk]);
			    
			    tjazp[0][i-li][j-lj][0] =  0.25*(Hr[i][j][uk] + Hr[i+1][j][uk] + Hr[i][j][uk-1] + Hr[i+1][j][uk-1]);
			    tjrzp[0][i-li][j-lj][0] = -0.25*(Ha[i][j][uk] + Ha[i][j+1][uk] + Ha[i][j][uk-1] + Ha[i][j+1][uk-1]);
			    
			    tmrzn[0][i-li][j-lj][0] = -0.5*(Ea[i][j][lk] + Ea[i+1][j][lk]);
			    tmazn[0][i-li][j-lj][0] =  0.5*(Er[i][j][lk] + Er[i][j+1][lk]);
			    
			    tjazn[0][i-li][j-lj][0] = -0.25*(Hr[i][j][lk] + Hr[i+1][j][lk] + Hr[i][j][lk-1] + Hr[i+1][j][lk-1]);
			    tjrzn[0][i-li][j-lj][0] =  0.25*(Ha[i][j][lk] + Ha[i][j+1][lk] + Ha[i][j][lk-1] + Ha[i][j+1][lk-1]);
			    
			    //fourier transform
			    for (int mi=0; mi<frequencies.length; mi++){
				cjrzp[i-li][j-lj][0][mi].set(cjrzp[i-li][j-lj][0][mi].plus(exp_h[mi].times(tjrzp[0][i-li][j-lj][0])));
				cjazp[i-li][j-lj][0][mi].set(cjazp[i-li][j-lj][0][mi].plus(exp_h[mi].times(tjazp[0][i-li][j-lj][0]))); 
				
				cjrzn[i-li][j-lj][0][mi].set(cjrzn[i-li][j-lj][0][mi].plus(exp_h[mi].times(tjrzn[0][i-li][j-lj][0])));
				cjazn[i-li][j-lj][0][mi].set(cjazn[i-li][j-lj][0][mi].plus(exp_h[mi].times(tjazn[0][i-li][j-lj][0])));
				
				cmrzp[i-li][j-lj][0][mi].set(cmrzp[i-li][j-lj][0][mi].plus(exp_e[mi].times(tmrzp[0][i-li][j-lj][0])));
				cmazp[i-li][j-lj][0][mi].set(cmazp[i-li][j-lj][0][mi].plus(exp_e[mi].times(tmazp[0][i-li][j-lj][0])));
				
				cmrzn[i-li][j-lj][0][mi].set(cmrzn[i-li][j-lj][0][mi].plus(exp_e[mi].times(tmrzn[0][i-li][j-lj][0]))); 
				cmazn[i-li][j-lj][0][mi].set(cmazn[i-li][j-lj][0][mi].plus(exp_e[mi].times(tmazn[0][i-li][j-lj][0]))); 
			    }
			}
		    }
		}
	    }
	}        
    }
    
    public double[] getFreqArray(){
        return this.frequencies;
    }    
    
    public double getFreq(int index){
        return frequencies[index];
    }
    
    public int getFreqCount(){
        return this.frequencies.length;
    }
    
    public void setCellCountOuterBoundary(int count){
        this.outerBoundaryCellCount = count;
    }
    
    public void initFarFieldArrays(ProblemSpace ps){
        li = this.outerBoundaryCellCount;
	lj = this.outerBoundaryCellCount;
	lk = this.outerBoundaryCellCount;
	ui = ps.getNR() - this.outerBoundaryCellCount;
	uj = ps.getNA() - this.outerBoundaryCellCount;
	uk = ps.getNZ() - this.outerBoundaryCellCount;
        
        this.omegas = Common.genDouble1DArray(this.frequencies.length, 0.0);
        
	for (int i = 0; i < frequencies.length ; i++){
            omegas[i] = 2*Math.PI*frequencies[i];
	}
        
	tjrap = Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tjrzp = Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tjarp = Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tjazp = Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tjzrp = Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tjzap = Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tjran = Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tjrzn = Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tjarn = Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tjazn = Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tjzrn = Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tjzan = Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tmrap = Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tmrzp = Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tmarp = Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tmazp = Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tmzrp = Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tmzap = Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tmran = Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tmrzn = Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tmarn = Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tmazn = Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tmzrn = Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tmzan = Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
        
	cjrap = Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjrzp = Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cjarp = Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjazp = Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cjzrp = Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjzap = Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjran = Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjrzn = Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cjarn = Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjazn = Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cjzrn = Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjzan = Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmrap = Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmrzp = Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cmarp = Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmazp = Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cmzrp = Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmzap = Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmran = Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmrzn = Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cmarn = Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmazn = Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cmzrn = Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmzan = Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
    }
}
