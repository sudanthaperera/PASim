package pasim;

public class FarField{
    private double[] frequencies, omegas;;
    private int outerBoundaryCellCount, angleStepsCount;
    private double[] radiatedPower;
    private int li,lj,lk,ui,uj,uk;
    private double[] theta,phi;
    private double[][] farfield_dataTheta;
    private double[][] farfield_dataPhi;
    private Cell c;
    private boolean farFieldDisable;
    
    private double[][][][] tjxyp, tjxzp, tjyxp, tjyzp, tjzxp, tjzyp, tjxyn, tjxzn, tjyxn, tjyzn, tjzxn, tjzyn, tmxyp, tmxzp, tmyxp, tmyzp, tmzxp, tmzyp, tmxyn, tmxzn, tmyxn, tmyzn, tmzxn, tmzyn;
    private Complex[][][][] cjxyp, cjxzp, cjyxp, cjyzp, cjzxp, cjzyp, cjxyn, cjxzn, cjyxn, cjyzn, cjzxn, cjzyn, cmxyp, cmxzp, cmyxp, cmyzp, cmzxp, cmzyp, cmxyn, cmxzn, cmyxn, cmyzn, cmzxn, cmzyn;
    
    public FarField(double start,double stop,int count, Cell c){
        this.c = c;
        double step = (stop - start)/(count - 1);
        frequencies = new double[count];
        
        for (int i = 0; i < count; i++){
            frequencies[i] = start + i*step;
        }
        
        angleStepsCount = 360;
        theta =  Common.genDouble1DArray(angleStepsCount, 0.0);
        phi =  Common.genDouble1DArray(angleStepsCount, 0.0);
        
        for(int i=0;i<angleStepsCount;i++){
            theta[i]=i*Math.PI/angleStepsCount;
            phi[i]=i*(2*Math.PI/angleStepsCount)-2*Math.PI;
        }
    }
    
    public void setFarFieldDisable(boolean status){
        this.farFieldDisable = status;
    }
    
    public void setupFarFieldOutPut(double[] theta, double[] phi){
        this.angleStepsCount = theta.length;
        this.theta = theta;
        this.phi = phi;
    }
    
    public void setupTheta(double start,double end,double stepsCount){
        if(this.farFieldDisable){
            return;
        }
        start = (Math.PI/180)*start;
        end = (Math.PI/180)*end;
        double stepSize = (end-start)/(stepsCount-1);
        for(int i=0;i<angleStepsCount;i++){
            theta[i] = start + stepSize*i;
        }
    }
    
    public void setupPhi(double start,double end,double stepsCount){
        if(this.farFieldDisable){
            return;
        }
        start = (Math.PI/180)*start;
        end = (Math.PI/180)*end;
        double stepSize = (Math.PI/180)*(end-start)/stepsCount;
        for(int i=0;i<angleStepsCount;i++){
            phi[i] = start + stepSize*i;
        }
    }
    
    public void saveFarFieldThetaPhiPattern(){
        if(this.farFieldDisable){
            return;
        }
         Common.save2DArray(farfield_dataPhi, "FarFieldPhiPattern");
         Common.save2DArray(farfield_dataTheta, "FarFieldThetaPattern");
         Common.save1DArray(phi, "FarFIeldPhi");
         Common.save1DArray(theta, "FarFieldTheta");
    }
    
    public void calFarField(){
        
        if(this.farFieldDisable){
            return;
        }
        
        Complex[] exp_jk_rpr =  Common.genComplex1DArray(angleStepsCount, new Complex());
        double[] dx_sinth_cosphi =  Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dy_sinth_sinphi =  Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dz_costh =  Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dy_dz_costh_sinphi =  Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dy_dz_sinth =  Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dy_dz_cosphi =  Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dx_dz_costh_cosphi =  Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dx_dz_sinth =  Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dx_dz_sinphi =  Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dx_dy_costh_cosphi =  Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dx_dy_costh_sinphi =  Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dx_dy_sinphi =  Common.genDouble1DArray(angleStepsCount, 0.0);
        double[] dx_dy_cosphi =  Common.genDouble1DArray(angleStepsCount, 0.0);
        
        //double[][] farfield_dirTheta =  Common.genDouble2DArray(frequencies.length,angleStepsCount,0.0);
        //double[][] farfield_dir =  Common.genDouble2DArray(frequencies.length,angleStepsCount,0.0);
        //double[][] farfield_dirPhi =  Common.genDouble2DArray(frequencies.length,angleStepsCount,0.0);
        
        farfield_dataTheta =  Common.genDouble2DArray(frequencies.length,angleStepsCount,0.0);
        farfield_dataPhi =  Common.genDouble2DArray(frequencies.length,angleStepsCount,0.0);
                
        for (int i = 0;i<angleStepsCount;i++){
            dx_sinth_cosphi[i] = c.getDeltaX()*Math.sin(theta[i])*Math.cos(phi[i]);
            dy_sinth_sinphi[i] = c.getDeltaY()*Math.sin(theta[i])*Math.sin(phi[i]);
            dz_costh[i] = c.getDeltaZ()*Math.cos(theta[i]);
            dy_dz_costh_sinphi[i] = c.getDeltaY()*c.getDeltaZ()*Math.cos(theta[i])*Math.sin(phi[i]);
            dy_dz_sinth[i] = c.getDeltaY()*c.getDeltaZ()*Math.sin(theta[i]);
            dy_dz_cosphi[i] = c.getDeltaY()*c.getDeltaZ()*Math.cos(phi[i]);
            dx_dz_costh_cosphi[i] = c.getDeltaX()*c.getDeltaZ()*Math.cos(theta[i])*Math.cos(phi[i]);
            dx_dz_sinth[i] = c.getDeltaX()*c.getDeltaZ()*Math.sin(theta[i]);
            dx_dz_sinphi[i] = c.getDeltaX()*c.getDeltaZ()*Math.sin(phi[i]);
            dx_dy_costh_cosphi[i] = c.getDeltaX()*c.getDeltaY()*Math.cos(theta[i])*Math.cos(phi[i]);
            dx_dy_costh_sinphi[i] = c.getDeltaX()*c.getDeltaY()*Math.cos(theta[i])*Math.sin(phi[i]);
            dx_dy_sinphi[i] = c.getDeltaX()*c.getDeltaY()*Math.sin(phi[i]);
            dx_dy_cosphi[i] = c.getDeltaX()*c.getDeltaY()*Math.cos(phi[i]);
        }
        
        double ci = 0.5*(ui+li);
        double cj = 0.5*(uj+lj);
        double ck = 0.5*(uk+lk);
        
        for (int mi=0; mi<this.frequencies.length; mi++){
                    
            Complex[] Ntheta =  Common.genComplex1DArray(angleStepsCount, new Complex());
            Complex[] Ltheta =  Common.genComplex1DArray(angleStepsCount, new Complex());
            Complex[] Nphi =  Common.genComplex1DArray(angleStepsCount, new Complex());
            Complex[] Lphi =  Common.genComplex1DArray(angleStepsCount, new Complex());
            double[] rpr =  Common.genDouble1DArray(angleStepsCount,1);
        
            double k;
            
            System.out.println("Calculating directivity for ".concat(String.valueOf(this.frequencies[mi]).concat(" Hz")));
            k = 2*Math.PI*frequencies[mi]*Math.sqrt(Constants.MU0*Constants.EPS0);
            for (int nj=lj;nj<uj;nj++){
                for (int nk=lk;nk<uk;nk++){
                    for (int i=0;i<angleStepsCount;i++){
                        //for +ax direction
                        rpr[i] = (ui - ci)*dx_sinth_cosphi[i] + (nj-cj+0.5)*dy_sinth_sinphi[i] + (nk-ck+0.5)*dz_costh[i];
                        exp_jk_rpr[i].set( Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i])));
                        Ntheta[i].set(Ntheta[i].plus((cjyxp[0][nj-lj][nk-lk][mi].times(dy_dz_costh_sinphi[i]).minus(cjzxp[0][nj-lj][nk-lk][mi].times(dy_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Ltheta[i].set(Ltheta[i].plus((cmyxp[0][nj-lj][nk-lk][mi].times(dy_dz_costh_sinphi[i]).minus(cmzxp[0][nj-lj][nk-lk][mi].times(dy_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Nphi[i].set(Nphi[i].plus((cjyxp[0][nj-lj][nk-lk][mi].times(dy_dz_cosphi[i])).times(exp_jk_rpr[i])));
                        Lphi[i].set(Lphi[i].plus((cmyxp[0][nj-lj][nk-lk][mi].times(dy_dz_cosphi[i])).times(exp_jk_rpr[i])));
                        
                        //for -ax direction
                        rpr[i] = (li - ci)*dx_sinth_cosphi[i] + (nj-cj+0.5)*dy_sinth_sinphi[i] + (nk-ck+0.5)*dz_costh[i];
                        exp_jk_rpr[i].set( Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i])));
                        Ntheta[i].set(Ntheta[i].plus((cjyxn[0][nj-lj][nk-lk][mi].times(dy_dz_costh_sinphi[i]).minus(cjzxn[0][nj-lj][nk-lk][mi].times(dy_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Ltheta[i].set(Ltheta[i].plus(cmyxn[0][nj-lj][nk-lk][mi].times(dy_dz_costh_sinphi[i]).minus(cmzxn[0][nj-lj][nk-lk][mi].times(dy_dz_sinth[i])).times(exp_jk_rpr[i])));
                        Nphi[i].set(Nphi[i].plus(cjyxn[0][nj-lj][nk-lk][mi].times(dy_dz_cosphi[i]).times(exp_jk_rpr[i])));
                        Lphi[i].set(Lphi[i].plus(cmyxn[0][nj-lj][nk-lk][mi].times(dy_dz_cosphi[i]).times(exp_jk_rpr[i])));
                    }
                }
            }
           
            for (int ni=li;ni<ui;ni++){
                for (int nk=lk;nk<uk;nk++){
                    for (int i=0;i<angleStepsCount;i++){
                        //for +ay direction
                        rpr[i] = (ni - ci + 0.5)*dx_sinth_cosphi[i] + (uj-cj)*dy_sinth_sinphi[i] + (nk-ck+0.5)*dz_costh[i];
                        exp_jk_rpr[i].set( Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i])));
                        Ntheta[i].set( Ntheta[i].plus((cjxyp[ni-li][0][nk-lk][mi].times(dx_dz_costh_cosphi[i]).minus(cjzyp[ni-li][0][nk-lk][mi].times(dx_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Ltheta[i].set(Ltheta[i].plus((cmxyp[ni-li][0][nk-lk][mi].times(dx_dz_costh_cosphi[i]).minus(cmzyp[ni-li][0][nk-lk][mi].times(dx_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Nphi[i].set(Nphi[i].plus(((cjxyp[ni-li][0][nk-lk][mi].times(-1.0)).times(dx_dz_sinphi[i])).times(exp_jk_rpr[i])));
                        Lphi[i].set(Lphi[i].plus(((cmxyp[ni-li][0][nk-lk][mi].times(-1.0)).times(dx_dz_sinphi[i])).times(exp_jk_rpr[i])));

                        //for -ay direction
                        rpr[i] = (ni - ci + 0.5)*dx_sinth_cosphi[i] + (lj-cj)*dy_sinth_sinphi[i] + (nk-ck+0.5)*dz_costh[i];
                        exp_jk_rpr[i].set( Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i])));
                        Ntheta[i].set(Ntheta[i].plus(((cjxyn[ni-li][0][nk-lk][mi].times(dx_dz_costh_cosphi[i])).minus(cjzyn[ni-li][0][nk-lk][mi].times(dx_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Ltheta[i].set(Ltheta[i].plus(((cmxyn[ni-li][0][nk-lk][mi].times(dx_dz_costh_cosphi[i])).minus(cmzyn[ni-li][0][nk-lk][mi].times(dx_dz_sinth[i]))).times(exp_jk_rpr[i])));
                        Nphi[i].set(Nphi[i].plus(((cjxyn[ni-li][0][nk-lk][mi].times(-1.0)).times(dx_dz_sinphi[i])).times(exp_jk_rpr[i])));
                        Lphi[i].set(Lphi[i].plus(((cmxyn[ni-li][0][nk-lk][mi].times(-1.0)).times(dx_dz_sinphi[i])).times(exp_jk_rpr[i])));
                    }
                }
            }
            
            for (int ni=li;ni<ui;ni++){
                for (int nj=lj;nj<uj;nj++){
                    for (int i=0;i<angleStepsCount;i++){
                        //for +az direction

                        rpr[i] = (ni-ci+0.5)*dx_sinth_cosphi[i] + (nj - cj + 0.5)*dy_sinth_sinphi[i] + (uk-ck)*dz_costh[i];
                        exp_jk_rpr[i].set( Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i])));
                        Ntheta[i].set(Ntheta[i].plus(((cjxzp[ni-li][nj-lj][0][mi].times(dx_dy_costh_cosphi[i])).plus(cjyzp[ni-li][nj-lj][0][mi].times(dx_dy_costh_sinphi[i]))).times(exp_jk_rpr[i])));
                        Ltheta[i].set(Ltheta[i].plus((((cmxzp[ni-li][nj-lj][0][mi].times(dx_dy_costh_cosphi[i])).plus((cmyzp[ni-li][nj-lj][0][mi].times(dx_dy_costh_sinphi[i])))).times(exp_jk_rpr[i]))));
                        Nphi[i].set(Nphi[i].plus((((cjxzp[ni-li][nj-lj][0][mi].times(-1.0)).times(dx_dy_sinphi[i])).plus((cjyzp[ni-li][nj-lj][0][mi].times(dx_dy_cosphi[i])))).times(exp_jk_rpr[i])));
                        Lphi[i].set(Lphi[i].plus((((cmxzp[ni-li][nj-lj][0][mi].times(-1.0)).times(dx_dy_sinphi[i])).plus(cmyzp[ni-li][nj-lj][0][mi].times(dx_dy_cosphi[i]))).times(exp_jk_rpr[i])));

                        //for -az direction
                        rpr[i] = (ni-ci+0.5)*dx_sinth_cosphi[i] + (nj - cj + 0.5)*dy_sinth_sinphi[i] + (lk-ck)*dz_costh[i];
                        exp_jk_rpr[i].set( Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i])));
                        Ntheta[i].set(Ntheta[i].plus(((cjxzn[ni-li][nj-lj][0][mi].times(dx_dy_costh_cosphi[i])).plus(cjyzn[ni-li][nj-lj][0][mi].times(dx_dy_costh_sinphi[i]))).times(exp_jk_rpr[i])));
                        Ltheta[i].set(Ltheta[i].plus(((cmxzn[ni-li][nj-lj][0][mi].times(dx_dy_costh_cosphi[i])).plus(cmyzn[ni-li][nj-lj][0][mi].times(dx_dy_costh_sinphi[i]))).times(exp_jk_rpr[i])));
                        Nphi[i].set(Nphi[i].plus((((cjxzn[ni-li][nj-lj][0][mi].times(-1.0)).times(dx_dy_sinphi[i])).plus(cjyzn[ni-li][nj-lj][0][mi].times(dx_dy_cosphi[i]))).times(exp_jk_rpr[i])));
                        Lphi[i].set(Lphi[i].plus((((cmxzn[ni-li][nj-lj][0][mi].times(-1.0)).times(dx_dy_sinphi[i])).plus(cmyzn[ni-li][nj-lj][0][mi].times(dx_dy_cosphi[i]))).times(exp_jk_rpr[i])));                        
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
        if(this.farFieldDisable){
            return;
        }
         Common.save1DArray(radiatedPower, "RadiatedPower");
    }
    
    public void calTotalRadiatedPower(){
        if(this.farFieldDisable){
            return;
        }
	radiatedPower =  Common.genDouble1DArray(frequencies.length,0.0);

        Complex[] sum_of_cxn =  Common.genComplex1DArray(frequencies.length,new Complex());
        Complex[] sum_of_cxp =  Common.genComplex1DArray(frequencies.length,new Complex());
        Complex[] sum_of_cyn =  Common.genComplex1DArray(frequencies.length,new Complex());
        Complex[] sum_of_cyp =  Common.genComplex1DArray(frequencies.length,new Complex());
        Complex[] sum_of_czn =  Common.genComplex1DArray(frequencies.length,new Complex());
        Complex[] sum_of_czp =  Common.genComplex1DArray(frequencies.length,new Complex());
        

	for (int ind=0; ind < frequencies.length; ind++){
            for (int i=0;i<ui-li;i++){
		for (int j=0;j<uj-lj;j++){
                    sum_of_czp[ind].set(sum_of_czp[ind].plus((cmyzp[i][j][0][ind].times(cjxzp[i][j][0][ind].conj())).minus(cmxzp[i][j][0][ind].times(cjyzp[i][j][0][ind].conj()))));
                    sum_of_czn[ind].set(sum_of_czn[ind].plus((cmyzn[i][j][0][ind].times(cjxzn[i][j][0][ind].conj())).minus(cmxzn[i][j][0][ind].times(cjyzn[i][j][0][ind].conj()))));
		}
            }

            for (int i=0;i<ui-li;i++){
		for (int k=0;k<uk-lk;k++){
                    sum_of_cyp[ind].set(sum_of_cyp[ind].plus((cmxyp[i][0][k][ind].times(cjzyp[i][0][k][ind].conj())).minus(cmzyp[i][0][k][ind].times(cjxyp[i][0][k][ind].conj()))));
                    sum_of_cyn[ind].set(sum_of_cyn[ind].plus((cmxyn[i][0][k][ind].times(cjzyn[i][0][k][ind].conj())).minus(cmzyn[i][0][k][ind].times(cjxyn[i][0][k][ind].conj()))));
		}
            }

            for (int j=0;j<uj-lj;j++){
		for (int k=0;k<uk-lk;k++){
                    sum_of_cxp[ind].set(sum_of_cxp[ind].plus((cmzxp[0][j][k][ind].times(cjyxp[0][j][k][ind].conj())).minus(cmyxp[0][j][k][ind].times(cjzxp[0][j][k][ind].conj()))));
                    sum_of_cxn[ind].set(sum_of_cxn[ind].plus((cmzxn[0][j][k][ind].times(cjyxn[0][j][k][ind].conj())).minus(cmyxn[0][j][k][ind].times(cjzxn[0][j][k][ind].conj()))));
		}
            }
            
            Complex temp1, temp2, temp3, temp4, temp5, temp6; 
            
            temp1 = sum_of_czp[ind].times(c.getDeltaX()*c.getDeltaY());
            temp2 = sum_of_czn[ind].times(c.getDeltaX()*c.getDeltaY());
            temp3 = sum_of_cyp[ind].times(c.getDeltaX()*c.getDeltaZ());
            temp4 = sum_of_cyn[ind].times(c.getDeltaX()*c.getDeltaZ());
            temp5 = sum_of_cxp[ind].times(c.getDeltaY()*c.getDeltaZ());
            temp6 = sum_of_cxn[ind].times(c.getDeltaY()*c.getDeltaZ());
            
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
    
    public void CalculateJandM(int timeIndex, EField E, HField H){
        if(this.farFieldDisable){
            return;
        }
        Complex zi = new Complex(0.0, 1.0);
	Complex[] exp_h =  Common.genComplex1DArray(frequencies.length, new Complex(0.0,0.0));
	Complex[] exp_e =  Common.genComplex1DArray(frequencies.length, new Complex(0.0,0.0));
        
	for (int mi=0; mi<frequencies.length; mi++){
	    exp_h[mi] = ( Common.complexExp(new Complex(0,-omegas[mi]*(timeIndex + 0.5)*c.getDeltaT()))).times(c.getDeltaT());
	    exp_e[mi] = ( Common.complexExp(new Complex(0,-omegas[mi]*(timeIndex + 1)*c.getDeltaT()))).times(c.getDeltaT());
	}
	
	if (frequencies.length > 0){
	    for(int i = li; i <= ui; i++){
		for (int j = lj; j <= uj; j++){
		    for (int k = lk; k <= uk; k++){
			if (i == ui && j != uj && k != uk){
			    tmyxp[0][0][j-lj][k-lk] =  0.5*(E.getEZ(ui,j,k) + E.getEZ(ui,j+1,k));
			    tmzxp[0][0][j-lj][k-lk] = -0.5*(E.getEY(ui,j,k) + E.getEY(ui,j,k+1));
			    
			    tjyxp[0][0][j-lj][k-lk] = -0.25*(H.getHZ(ui,j,k) + H.getHZ(ui,j,k+1) + H.getHZ(ui-1,j,k) + H.getHZ(ui-1,j,k+1));
			    tjzxp[0][0][j-lj][k-lk] =  0.25*(H.getHY(ui,j,k) + H.getHY(ui,j+1,k) + H.getHY(ui-1,j,k) + H.getHY(ui-1,j+1,k));
			    
			    tmyxn[0][0][j-lj][k-lk] = -0.5*(E.getEZ(li,j,k) + E.getEZ(li,j+1,k));
			    tmzxn[0][0][j-lj][k-lk] =  0.5*(E.getEY(li,j,k) + E.getEY(li,j,k+1));
			    
			    tjyxn[0][0][j-lj][k-lk] =  0.25*(H.getHZ(li,j,k) + H.getHZ(li,j,k+1) + H.getHZ(li-1,j,k) + H.getHZ(li-1,j,k+1));
			    tjzxn[0][0][j-lj][k-lk] = -0.25*(H.getHY(li,j,k) + H.getHY(li,j+1,k) + H.getHY(li-1,j,k) + H.getHY(li-1,j+1,k));
			    
			    //fourier transform
			    for (int mi=0; mi<frequencies.length; mi++){
				cjyxp[0][j-lj][k-lk][mi].set(cjyxp[0][j-lj][k-lk][mi].plus(exp_h[mi].times(tjyxp[0][0][j-lj][k-lk])));
				cjzxp[0][j-lj][k-lk][mi].set(cjzxp[0][j-lj][k-lk][mi].plus(exp_h[mi].times(tjzxp[0][0][j-lj][k-lk])));
				
				cjyxn[0][j-lj][k-lk][mi].set(cjyxn[0][j-lj][k-lk][mi].plus(exp_h[mi].times(tjyxn[0][0][j-lj][k-lk])));
				cjzxn[0][j-lj][k-lk][mi].set(cjzxn[0][j-lj][k-lk][mi].plus(exp_h[mi].times(tjzxn[0][0][j-lj][k-lk])));
				
				cmyxp[0][j-lj][k-lk][mi].set(cmyxp[0][j-lj][k-lk][mi].plus(exp_e[mi].times(tmyxp[0][0][j-lj][k-lk])));
				cmzxp[0][j-lj][k-lk][mi].set(cmzxp[0][j-lj][k-lk][mi].plus(exp_e[mi].times(tmzxp[0][0][j-lj][k-lk]))); 
				
				cmyxn[0][j-lj][k-lk][mi].set(cmyxn[0][j-lj][k-lk][mi].plus(exp_e[mi].times(tmyxn[0][0][j-lj][k-lk])));
				cmzxn[0][j-lj][k-lk][mi].set(cmzxn[0][j-lj][k-lk][mi].plus(exp_e[mi].times(tmzxn[0][0][j-lj][k-lk]))); 
			    }
			}
			
			if (i != ui && j == uj && k != uk){
			    tmxyp[0][i-li][0][k-lk] = -0.5*(E.getEZ(i,uj,k) + E.getEZ(i+1,uj,k));
			    tmzyp[0][i-li][0][k-lk] =  0.5*(E.getEX(i,uj,k) + E.getEX(i,uj,k+1));
			    
			    tjzyp[0][i-li][0][k-lk] = -0.25*(H.getHX(i,uj,k) + H.getHX(i+1,uj,k) + H.getHX(i,uj-1,k) + H.getHX(i+1,uj-1,k));
			    tjxyp[0][i-li][0][k-lk] =  0.25*(H.getHZ(i,uj,k) + H.getHZ(i,uj,k+1) + H.getHZ(i,uj-1,k) + H.getHZ(i,uj-1,k+1));
			    
			    tmxyn[0][i-li][0][k-lk] =  0.5*(E.getEZ(i,lj,k) + E.getEZ(i+1,lj,k));
			    tmzyn[0][i-li][0][k-lk] = -0.5*(E.getEX(i,lj,k) + E.getEX(i,lj,k+1));
			    
			    tjzyn[0][i-li][0][k-lk] =  0.25*(H.getHX(i,lj,k) + H.getHX(i+1,lj,k) + H.getHX(i,lj-1,k) + H.getHX(i+1,lj-1,k));
			    tjxyn[0][i-li][0][k-lk] = -0.25*(H.getHZ(i,lj,k) + H.getHZ(i,lj,k+1) + H.getHZ(i,lj-1,k) + H.getHZ(i,lj-1,k+1));
			    
			    //fourier transform
			    for (int mi=0; mi<frequencies.length; mi++){
				cjxyp[i-li][0][k-lk][mi].set(cjxyp[i-li][0][k-lk][mi].plus(exp_h[mi].times(tjxyp[0][i-li][0][k-lk])));
				cjzyp[i-li][0][k-lk][mi].set(cjzyp[i-li][0][k-lk][mi].plus(exp_h[mi].times(tjzyp[0][i-li][0][k-lk])));
				
				cjxyn[i-li][0][k-lk][mi].set(cjxyn[i-li][0][k-lk][mi].plus(exp_h[mi].times(tjxyn[0][i-li][0][k-lk])));
				cjzyn[i-li][0][k-lk][mi].set(cjzyn[i-li][0][k-lk][mi].plus(exp_h[mi].times(tjzyn[0][i-li][0][k-lk])));
				
				cmxyp[i-li][0][k-lk][mi].set(cmxyp[i-li][0][k-lk][mi].plus(exp_e[mi].times(tmxyp[0][i-li][0][k-lk]))); 
				cmzyp[i-li][0][k-lk][mi].set(cmzyp[i-li][0][k-lk][mi].plus(exp_e[mi].times(tmzyp[0][i-li][0][k-lk]))); 
				
				cmxyn[i-li][0][k-lk][mi].set(cmxyn[i-li][0][k-lk][mi].plus(exp_e[mi].times(tmxyn[0][i-li][0][k-lk]))); 
				cmzyn[i-li][0][k-lk][mi].set(cmzyn[i-li][0][k-lk][mi].plus(exp_e[mi].times(tmzyn[0][i-li][0][k-lk]))); 
			    }
			}
			
			if (i != ui && j != uj && k == uk){
			    tmxzp[0][i-li][j-lj][0] =  0.5*(E.getEY(i,j,uk) + E.getEY(i+1,j,uk));
			    tmyzp[0][i-li][j-lj][0] = -0.5*(E.getEX(i,j,uk) + E.getEX(i,j+1,uk));
			    
			    tjyzp[0][i-li][j-lj][0] =  0.25*(H.getHX(i,j,uk) + H.getHX(i+1,j,uk) + H.getHX(i,j,uk-1) + H.getHX(i+1,j,uk-1));
			    tjxzp[0][i-li][j-lj][0] = -0.25*(H.getHY(i,j,uk) + H.getHY(i,j+1,uk) + H.getHY(i,j,uk-1) + H.getHY(i,j+1,uk-1));
			    
			    tmxzn[0][i-li][j-lj][0] = -0.5*(E.getEY(i,j,lk) + E.getEY(i+1,j,lk));
			    tmyzn[0][i-li][j-lj][0] =  0.5*(E.getEX(i,j,lk) + E.getEX(i,j+1,lk));
			    
			    tjyzn[0][i-li][j-lj][0] = -0.25*(H.getHX(i,j,lk) + H.getHX(i+1,j,lk) + H.getHX(i,j,lk-1) + H.getHX(i+1,j,lk-1));
			    tjxzn[0][i-li][j-lj][0] =  0.25*(H.getHY(i,j,lk) + H.getHY(i,j+1,lk) + H.getHY(i,j,lk-1) + H.getHY(i,j+1,lk-1));
			    
			    //fourier transform
			    for (int mi=0; mi<frequencies.length; mi++){
				cjxzp[i-li][j-lj][0][mi].set(cjxzp[i-li][j-lj][0][mi].plus(exp_h[mi].times(tjxzp[0][i-li][j-lj][0])));
				cjyzp[i-li][j-lj][0][mi].set(cjyzp[i-li][j-lj][0][mi].plus(exp_h[mi].times(tjyzp[0][i-li][j-lj][0]))); 
				
				cjxzn[i-li][j-lj][0][mi].set(cjxzn[i-li][j-lj][0][mi].plus(exp_h[mi].times(tjxzn[0][i-li][j-lj][0])));
				cjyzn[i-li][j-lj][0][mi].set(cjyzn[i-li][j-lj][0][mi].plus(exp_h[mi].times(tjyzn[0][i-li][j-lj][0])));
				
				cmxzp[i-li][j-lj][0][mi].set(cmxzp[i-li][j-lj][0][mi].plus(exp_e[mi].times(tmxzp[0][i-li][j-lj][0])));
				cmyzp[i-li][j-lj][0][mi].set(cmyzp[i-li][j-lj][0][mi].plus(exp_e[mi].times(tmyzp[0][i-li][j-lj][0])));
				
				cmxzn[i-li][j-lj][0][mi].set(cmxzn[i-li][j-lj][0][mi].plus(exp_e[mi].times(tmxzn[0][i-li][j-lj][0]))); 
				cmyzn[i-li][j-lj][0][mi].set(cmyzn[i-li][j-lj][0][mi].plus(exp_e[mi].times(tmyzn[0][i-li][j-lj][0]))); 
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
    
    public void initFarFieldArrays(ProblemSpace ps,Boundary b){
        if(b.getCPMLxn())
            li = this.outerBoundaryCellCount;
        if(b.getCPMLyn())
            lj = this.outerBoundaryCellCount;
        if(b.getCPMLzn())
            lk = this.outerBoundaryCellCount;
        if(b.getCPMLxp())
            ui = ps.getNX() - this.outerBoundaryCellCount;
        if(b.getCPMLyp())
            uj = ps.getNY() - this.outerBoundaryCellCount;
        if(b.getCPMLzp())
            uk = ps.getNZ() - this.outerBoundaryCellCount;
        if(b.getPBCX()){
            li = 3;
            ui = ps.getNX()-3;
        }
        if(b.getPBCY()){
            lj = 3;
            uj = ps.getNY()-3;
        }
        if(b.getPBCZ()){
            lk = 3;
            uk = ps.getNZ()-3;
        }

        this.omegas =  Common.genDouble1DArray(this.frequencies.length, 0.0);
        
	for (int i = 0; i < frequencies.length ; i++){
            omegas[i] = 2*Math.PI*frequencies[i];
	}
        
	tjxyp =  Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tjxzp =  Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tjyxp =  Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tjyzp =  Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tjzxp =  Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tjzyp =  Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tjxyn =  Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tjxzn =  Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tjyxn =  Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tjyzn =  Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tjzxn =  Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tjzyn =  Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tmxyp =  Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tmxzp =  Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tmyxp =  Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tmyzp =  Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tmzxp =  Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tmzyp =  Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tmxyn =  Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
	tmxzn =  Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tmyxn =  Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tmyzn =  Common.genDouble4DArray(1,ui-li,uj-lj,1, 0.0);
	tmzxn =  Common.genDouble4DArray(1,1,uj-lj,uk-lk, 0.0);
	tmzyn =  Common.genDouble4DArray(1,ui-li,1,uk-lk, 0.0);
        
	cjxyp =  Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjxzp =  Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cjyxp =  Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjyzp =  Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cjzxp =  Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjzyp =  Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjxyn =  Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjxzn =  Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cjyxn =  Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjyzn =  Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cjzxn =  Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cjzyn =  Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmxyp =  Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmxzp =  Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cmyxp =  Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmyzp =  Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cmzxp =  Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmzyp =  Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmxyn =  Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmxzn =  Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cmyxn =  Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmyzn =  Common.genComplex4DArray(ui-li,uj-lj,1,frequencies.length, new Complex(0.0,0.0));
	cmzxn =  Common.genComplex4DArray(1,uj-lj,uk-lk,frequencies.length, new Complex(0.0,0.0));
	cmzyn =  Common.genComplex4DArray(ui-li,1,uk-lk,frequencies.length, new Complex(0.0,0.0));
    }
}
