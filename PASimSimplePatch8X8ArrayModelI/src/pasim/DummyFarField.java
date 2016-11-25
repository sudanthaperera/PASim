package pasim;

public class DummyFarField{
    private double[] frequencies, omegas;;
    private int outerBoundaryCellCount, angleStepsCount, thetaStepsCount, phiStepsCount;
    private double[] radiatedPower;
    private int li,lj,lk,ui,uj,uk;
    private double[] theta,phi;
    private double[][] farfield_dataTheta;
    private double[][] farfield_dataPhi;
    private double[][][] farfield_dataTheta3D;
    private double[][][] farfield_dataPhi3D;
    private double[][][] farfield_dataThetaPhase3D;
    private double[][][] farfield_dataPhiPhase3D;
    private Cell c;
    
    private Complex[][][][] cjxyp, cjxzp, cjyxp, cjyzp, cjzxp, cjzyp, cjxyn, cjxzn, cjyxn, cjyzn, cjzxn, cjzyn, cmxyp, cmxzp, cmyxp, cmyzp, cmzxp, cmzyp, cmxyn, cmxzn, cmyxn, cmyzn, cmzxn, cmzyn;
    
    public void printSize(String ArrayName){
        System.out.printf(ArrayName.concat(" ---- ").concat("(X , Y , Z) = (%d , %d , %d) \n"),(this.ui - this.li),(this.uj - this.lj),(this.uk - this.lk));
    }
    
    public void setCjxyp(Complex[][][][] val){
        cjxyp = val;
    }
    
    public void setCjxzp(Complex[][][][] val){
        cjxzp = val;
    }
    
    public void setCjyxp(Complex[][][][] val){
        cjyxp = val;
    }
    
    public void setCjyzp(Complex[][][][] val){
        cjyzp = val;
    }
    
    public void setCjzxp(Complex[][][][] val){
        cjzxp = val;
    }
    
    public void setCjzyp(Complex[][][][] val){
        cjzyp = val;
    }
    
    public void setCjxyn(Complex[][][][] val){
        cjxyn = val;
    }
    
    public void setCjxzn(Complex[][][][] val){
        cjxzn = val;
    }
    
    public void setCjyxn(Complex[][][][] val){
        cjyxn = val;
    }
    
    public void setCjyzn(Complex[][][][] val){
        cjyzn = val;
    }
    
    public void setCjzxn(Complex[][][][] val){
        cjzxn = val;
    }
    
    public void setCjzyn(Complex[][][][] val){
        cjzyn = val;
    }
    
    public void setCmxyp(Complex[][][][] val){
        cmxyp = val;
    }
    
    public void setCmxzp(Complex[][][][] val){
        cmxzp = val;
    }
    
    public void setCmyxp(Complex[][][][] val){
        cmyxp = val;
    }
    
    public void setCmyzp(Complex[][][][] val){
        cmyzp = val;
    }
    
    public void setCmzxp(Complex[][][][] val){
        cmzxp = val;
    }
    
    public void setCmzyp(Complex[][][][] val){
        cmzyp = val;
    }
    
    public void setCmxyn(Complex[][][][] val){
        cmxyn = val;
    }
    
    public void setCmxzn(Complex[][][][] val){
        cmxzn = val;
    }
    
    public void setCmyxn(Complex[][][][] val){
        cmyxn = val;
    }
    
    public void setCmyzn(Complex[][][][] val){
        cmyzn = val;
    }
    
    public void setCmzxn(Complex[][][][] val){
        cmzxn = val;
    }
    
    public void setCmzyn(Complex[][][][] val){
        cmzyn = val;
    }
    
    public DummyFarField(double singleFreq, Cell c){
        this.c = c;
        frequencies = new double[1];
        frequencies[0] = singleFreq;
        angleStepsCount = 360;
        theta =  Common.genDouble1DArray(angleStepsCount, 0.0);
        phi =  Common.genDouble1DArray(angleStepsCount, 0.0);
        
        for(int i=0;i<angleStepsCount;i++){
            theta[i]=i*Math.PI/angleStepsCount;
            phi[i]=i*(2*Math.PI/angleStepsCount)-2*Math.PI;
        }
    }
    
    public DummyFarField(double start,double stop,int count, Cell c){
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
    
    public void setupFarFieldOutPut(double[] theta, double[] phi){
        this.angleStepsCount = theta.length;
        this.theta = theta;
        this.phi = phi;
    }
    
    public void setupTheta(double start,double end,int stepsCount){
        start = (Math.PI/180)*start;
        end = (Math.PI/180)*end;
        double stepSize = (end-start)/((double)stepsCount-1);
        theta = Common.genDouble1DArray(stepsCount, 0.0);
        thetaStepsCount = stepsCount;
        for(int i=0;i<thetaStepsCount;i++){
            theta[i] = start + stepSize*i;
        }
    }
    
    public void setupPhi(double start,double end,int stepsCount){
        start = (Math.PI/180)*start;
        end = (Math.PI/180)*end;
        double stepSize = (end-start)/((double)stepsCount-1);
        phi = Common.genDouble1DArray(stepsCount, 0.0);
        phiStepsCount = stepsCount;
        for(int i=0;i<phiStepsCount;i++){
            phi[i] = start + stepSize*i;
        }
    }
    
    public void saveFarFieldThetaPhiPattern(){
         Common.save2DArray(farfield_dataPhi, "FarFieldPhiPattern");
         Common.save2DArray(farfield_dataTheta, "FarFieldThetaPattern");
         Common.save1DArray(phi, "FarFIeldPhi");
         Common.save1DArray(theta, "FarFieldTheta");
    }
    
    public void saveFarFieldThetaPhiPattern(String FarFieldPhiPattern,String FarFieldThetaPattern,String FarFIeldPhi,String FarFieldTheta){
        Common.save2DArray(farfield_dataPhi, FarFieldPhiPattern);
        Common.save2DArray(farfield_dataTheta, FarFieldThetaPattern);
        Common.save1DArray(phi, FarFIeldPhi);
        Common.save1DArray(theta, FarFieldTheta);
    }
    
    public void saveFarFieldThetaPhi3DPattern(String FarFieldPhiPattern,String FarFieldThetaPattern,String FarFIeldPhi,String FarFieldTheta){
        double[][] dataTheta = Common.genDouble2DArray(thetaStepsCount, phiStepsCount, 0.0);
        double[][] dataPhi = Common.genDouble2DArray(thetaStepsCount, phiStepsCount, 0.0);
        double[][] dataThetaPhase = Common.genDouble2DArray(thetaStepsCount, phiStepsCount, 0.0);
        double[][] dataPhiPhase = Common.genDouble2DArray(thetaStepsCount, phiStepsCount, 0.0);
        
        for (int freqIndex = 0;freqIndex<this.frequencies.length;freqIndex++){
            for (int i = 0;i<thetaStepsCount;i++){
                for (int j = 0;j<phiStepsCount;j++){
                    dataPhi[i][j] = farfield_dataPhi3D[freqIndex][i][j];
                    dataTheta[i][j] = farfield_dataTheta3D[freqIndex][i][j];
                    dataPhiPhase[i][j] = farfield_dataPhiPhase3D[freqIndex][i][j];
                    dataThetaPhase[i][j] = farfield_dataThetaPhase3D[freqIndex][i][j];
                }
            }
            Common.save2DArray(dataPhi, FarFieldPhiPattern.concat(String.valueOf(freqIndex)).concat("Magnitude"));
            Common.save2DArray(dataTheta, FarFieldThetaPattern.concat(String.valueOf(freqIndex)).concat("Magnitude"));
            Common.save2DArray(dataPhiPhase, FarFieldPhiPattern.concat(String.valueOf(freqIndex)).concat("Phase"));
            Common.save2DArray(dataThetaPhase, FarFieldThetaPattern.concat(String.valueOf(freqIndex)).concat("Phase"));
        }
        Common.save1DArray(phi, FarFIeldPhi);
        Common.save1DArray(theta, FarFieldTheta);
    }
    
    public void calFarField3D(){
        Complex[][] exp_jk_rpr = Common.genComplex2DArray(thetaStepsCount,phiStepsCount, new Complex());
        double[][] dx_sinth_cosphi = Common.genDouble2DArray(thetaStepsCount,phiStepsCount, 0.0);
        double[][] dy_sinth_sinphi = Common.genDouble2DArray(thetaStepsCount,phiStepsCount, 0.0);
        double[][] dz_costh = Common.genDouble2DArray(thetaStepsCount,phiStepsCount, 0.0);
        double[][] dy_dz_costh_sinphi = Common.genDouble2DArray(thetaStepsCount,phiStepsCount, 0.0);
        double[][] dy_dz_sinth = Common.genDouble2DArray(thetaStepsCount,phiStepsCount, 0.0);
        double[][] dy_dz_cosphi = Common.genDouble2DArray(thetaStepsCount,phiStepsCount, 0.0);
        double[][] dx_dz_costh_cosphi = Common.genDouble2DArray(thetaStepsCount,phiStepsCount, 0.0);
        double[][] dx_dz_sinth = Common.genDouble2DArray(thetaStepsCount,phiStepsCount, 0.0);
        double[][] dx_dz_sinphi = Common.genDouble2DArray(thetaStepsCount,phiStepsCount, 0.0);
        double[][] dx_dy_costh_cosphi = Common.genDouble2DArray(thetaStepsCount,phiStepsCount, 0.0);
        double[][] dx_dy_costh_sinphi = Common.genDouble2DArray(thetaStepsCount,phiStepsCount, 0.0);
        double[][] dx_dy_sinphi = Common.genDouble2DArray(thetaStepsCount,phiStepsCount, 0.0);
        double[][] dx_dy_cosphi = Common.genDouble2DArray(thetaStepsCount,phiStepsCount, 0.0);
        
        
        farfield_dataTheta3D = Common.genDouble3DArray(frequencies.length,thetaStepsCount,phiStepsCount,0.0);
        farfield_dataPhi3D = Common.genDouble3DArray(frequencies.length,thetaStepsCount,phiStepsCount,0.0);
        farfield_dataThetaPhase3D = Common.genDouble3DArray(frequencies.length,thetaStepsCount,phiStepsCount,0.0);
        farfield_dataPhiPhase3D = Common.genDouble3DArray(frequencies.length,thetaStepsCount,phiStepsCount,0.0);
                
        for (int i = 0;i<thetaStepsCount;i++){
            for (int j = 0;j<phiStepsCount;j++){
                dx_sinth_cosphi[i][j] = c.getDeltaX()*Math.sin(theta[i])*Math.cos(phi[j]);
                dy_sinth_sinphi[i][j] = c.getDeltaY()*Math.sin(theta[i])*Math.sin(phi[j]);
                dz_costh[i][j] = c.getDeltaZ()*Math.cos(theta[i]);
                dy_dz_costh_sinphi[i][j] = c.getDeltaY()*c.getDeltaZ()*Math.cos(theta[i])*Math.sin(phi[j]);
                dy_dz_sinth[i][j] = c.getDeltaY()*c.getDeltaZ()*Math.sin(theta[i]);
                dy_dz_cosphi[i][j] = c.getDeltaY()*c.getDeltaZ()*Math.cos(phi[j]);
                dx_dz_costh_cosphi[i][j] = c.getDeltaX()*c.getDeltaZ()*Math.cos(theta[i])*Math.cos(phi[j]);
                dx_dz_sinth[i][j] = c.getDeltaX()*c.getDeltaZ()*Math.sin(theta[i]);
                dx_dz_sinphi[i][j] = c.getDeltaX()*c.getDeltaZ()*Math.sin(phi[j]);
                dx_dy_costh_cosphi[i][j] = c.getDeltaX()*c.getDeltaY()*Math.cos(theta[i])*Math.cos(phi[j]);
                dx_dy_costh_sinphi[i][j] = c.getDeltaX()*c.getDeltaY()*Math.cos(theta[i])*Math.sin(phi[j]);
                dx_dy_sinphi[i][j] = c.getDeltaX()*c.getDeltaY()*Math.sin(phi[j]);
                dx_dy_cosphi[i][j] = c.getDeltaX()*c.getDeltaY()*Math.cos(phi[j]);
            }
        }
        
        double ci = 0.5*(ui+li);
        double cj = 0.5*(uj+lj);
        double ck = 0.5*(uk+lk);
        
        for (int mi=0; mi<this.frequencies.length; mi++){
                    
            Complex[][] Ntheta = Common.genComplex2DArray(thetaStepsCount,phiStepsCount, new Complex());
            Complex[][] Ltheta = Common.genComplex2DArray(thetaStepsCount,phiStepsCount, new Complex());
            Complex[][] Nphi = Common.genComplex2DArray(thetaStepsCount,phiStepsCount, new Complex());
            Complex[][] Lphi = Common.genComplex2DArray(thetaStepsCount,phiStepsCount, new Complex());
            double[][] rpr = Common.genDouble2DArray(thetaStepsCount,phiStepsCount,1);
        
            double k;
            
            System.out.println("Calculating directivity for ".concat(String.valueOf(this.frequencies[mi]).concat(" Hz")));
            k = 2*Math.PI*frequencies[mi]*Math.sqrt(Constants.MU0*Constants.EPS0);
            for (int nj=lj;nj<uj;nj++){
                for (int nk=lk;nk<uk;nk++){
                    for (int i=0;i<thetaStepsCount;i++){
                        for (int j = 0;j<phiStepsCount;j++){
                            //for +ax direction
                            rpr[i][j] = (ui - ci)*dx_sinth_cosphi[i][j] + (nj-cj+0.5)*dy_sinth_sinphi[i][j] + (nk-ck+0.5)*dz_costh[i][j];
                            exp_jk_rpr[i][j].set(Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i][j])));
                            Ntheta[i][j].set(Ntheta[i][j].plus((cjyxp[0][nj-lj][nk-lk][mi].times(dy_dz_costh_sinphi[i][j]).minus(cjzxp[0][nj-lj][nk-lk][mi].times(dy_dz_sinth[i][j]))).times(exp_jk_rpr[i][j])));
                            Ltheta[i][j].set(Ltheta[i][j].plus((cmyxp[0][nj-lj][nk-lk][mi].times(dy_dz_costh_sinphi[i][j]).minus(cmzxp[0][nj-lj][nk-lk][mi].times(dy_dz_sinth[i][j]))).times(exp_jk_rpr[i][j])));
                            Nphi[i][j].set(Nphi[i][j].plus((cjyxp[0][nj-lj][nk-lk][mi].times(dy_dz_cosphi[i][j])).times(exp_jk_rpr[i][j])));
                            Lphi[i][j].set(Lphi[i][j].plus((cmyxp[0][nj-lj][nk-lk][mi].times(dy_dz_cosphi[i][j])).times(exp_jk_rpr[i][j])));
                        
                            //for -ax direction
                            rpr[i][j] = (li - ci)*dx_sinth_cosphi[i][j] + (nj-cj+0.5)*dy_sinth_sinphi[i][j] + (nk-ck+0.5)*dz_costh[i][j];
                            exp_jk_rpr[i][j].set(Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i][j])));
                            Ntheta[i][j].set(Ntheta[i][j].plus((cjyxn[0][nj-lj][nk-lk][mi].times(dy_dz_costh_sinphi[i][j]).minus(cjzxn[0][nj-lj][nk-lk][mi].times(dy_dz_sinth[i][j]))).times(exp_jk_rpr[i][j])));
                            Ltheta[i][j].set(Ltheta[i][j].plus(cmyxn[0][nj-lj][nk-lk][mi].times(dy_dz_costh_sinphi[i][j]).minus(cmzxn[0][nj-lj][nk-lk][mi].times(dy_dz_sinth[i][j])).times(exp_jk_rpr[i][j])));
                            Nphi[i][j].set(Nphi[i][j].plus(cjyxn[0][nj-lj][nk-lk][mi].times(dy_dz_cosphi[i][j]).times(exp_jk_rpr[i][j])));
                            Lphi[i][j].set(Lphi[i][j].plus(cmyxn[0][nj-lj][nk-lk][mi].times(dy_dz_cosphi[i][j]).times(exp_jk_rpr[i][j])));
                        }
                    }
                }
            }
           
            for (int ni=li;ni<ui;ni++){
                for (int nk=lk;nk<uk;nk++){
                    for (int i=0;i<thetaStepsCount;i++){
                        for (int j = 0;j<phiStepsCount;j++){
                            //for +ay direction
                            rpr[i][j] = (ni - ci + 0.5)*dx_sinth_cosphi[i][j] + (uj-cj)*dy_sinth_sinphi[i][j] + (nk-ck+0.5)*dz_costh[i][j];
                            exp_jk_rpr[i][j].set(Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i][j])));
                            Ntheta[i][j].set( Ntheta[i][j].plus((cjxyp[ni-li][0][nk-lk][mi].times(dx_dz_costh_cosphi[i][j]).minus(cjzyp[ni-li][0][nk-lk][mi].times(dx_dz_sinth[i][j]))).times(exp_jk_rpr[i][j])));
                            Ltheta[i][j].set(Ltheta[i][j].plus((cmxyp[ni-li][0][nk-lk][mi].times(dx_dz_costh_cosphi[i][j]).minus(cmzyp[ni-li][0][nk-lk][mi].times(dx_dz_sinth[i][j]))).times(exp_jk_rpr[i][j])));
                            Nphi[i][j].set(Nphi[i][j].plus(((cjxyp[ni-li][0][nk-lk][mi].times(-1.0)).times(dx_dz_sinphi[i][j])).times(exp_jk_rpr[i][j])));
                            Lphi[i][j].set(Lphi[i][j].plus(((cmxyp[ni-li][0][nk-lk][mi].times(-1.0)).times(dx_dz_sinphi[i][j])).times(exp_jk_rpr[i][j])));

                            //for -ay direction
                            rpr[i][j] = (ni - ci + 0.5)*dx_sinth_cosphi[i][j] + (lj-cj)*dy_sinth_sinphi[i][j] + (nk-ck+0.5)*dz_costh[i][j];
                            exp_jk_rpr[i][j].set(Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i][j])));
                            Ntheta[i][j].set(Ntheta[i][j].plus(((cjxyn[ni-li][0][nk-lk][mi].times(dx_dz_costh_cosphi[i][j])).minus(cjzyn[ni-li][0][nk-lk][mi].times(dx_dz_sinth[i][j]))).times(exp_jk_rpr[i][j])));
                            Ltheta[i][j].set(Ltheta[i][j].plus(((cmxyn[ni-li][0][nk-lk][mi].times(dx_dz_costh_cosphi[i][j])).minus(cmzyn[ni-li][0][nk-lk][mi].times(dx_dz_sinth[i][j]))).times(exp_jk_rpr[i][j])));
                            Nphi[i][j].set(Nphi[i][j].plus(((cjxyn[ni-li][0][nk-lk][mi].times(-1.0)).times(dx_dz_sinphi[i][j])).times(exp_jk_rpr[i][j])));
                            Lphi[i][j].set(Lphi[i][j].plus(((cmxyn[ni-li][0][nk-lk][mi].times(-1.0)).times(dx_dz_sinphi[i][j])).times(exp_jk_rpr[i][j])));
                        }
                    }
                }
            }
            
            for (int ni=li;ni<ui;ni++){
                for (int nj=lj;nj<uj;nj++){
                    for (int i=0;i<thetaStepsCount;i++){
                        for (int j = 0;j<phiStepsCount;j++){
                            //for +az direction
                            
                            rpr[i][j] = (ni-ci+0.5)*dx_sinth_cosphi[i][j] + (nj - cj + 0.5)*dy_sinth_sinphi[i][j] + (uk-ck)*dz_costh[i][j];
                            exp_jk_rpr[i][j].set(Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i][j])));
                            Ntheta[i][j].set(Ntheta[i][j].plus(((cjxzp[ni-li][nj-lj][0][mi].times(dx_dy_costh_cosphi[i][j])).plus(cjyzp[ni-li][nj-lj][0][mi].times(dx_dy_costh_sinphi[i][j]))).times(exp_jk_rpr[i][j])));
                            Ltheta[i][j].set(Ltheta[i][j].plus((((cmxzp[ni-li][nj-lj][0][mi].times(dx_dy_costh_cosphi[i][j])).plus((cmyzp[ni-li][nj-lj][0][mi].times(dx_dy_costh_sinphi[i][j])))).times(exp_jk_rpr[i][j]))));
                            Nphi[i][j].set(Nphi[i][j].plus((((cjxzp[ni-li][nj-lj][0][mi].times(-1.0)).times(dx_dy_sinphi[i][j])).plus((cjyzp[ni-li][nj-lj][0][mi].times(dx_dy_cosphi[i][j])))).times(exp_jk_rpr[i][j])));
                            Lphi[i][j].set(Lphi[i][j].plus((((cmxzp[ni-li][nj-lj][0][mi].times(-1.0)).times(dx_dy_sinphi[i][j])).plus(cmyzp[ni-li][nj-lj][0][mi].times(dx_dy_cosphi[i][j]))).times(exp_jk_rpr[i][j])));
                            
                            //for -az direction
                            rpr[i][j] = (ni-ci+0.5)*dx_sinth_cosphi[i][j] + (nj - cj + 0.5)*dy_sinth_sinphi[i][j] + (lk-ck)*dz_costh[i][j];
                            exp_jk_rpr[i][j].set(Common.complexExp((new Complex(0.0,-1.0)).times(k*rpr[i][j])));
                            Ntheta[i][j].set(Ntheta[i][j].plus(((cjxzn[ni-li][nj-lj][0][mi].times(dx_dy_costh_cosphi[i][j])).plus(cjyzn[ni-li][nj-lj][0][mi].times(dx_dy_costh_sinphi[i][j]))).times(exp_jk_rpr[i][j])));
                            Ltheta[i][j].set(Ltheta[i][j].plus(((cmxzn[ni-li][nj-lj][0][mi].times(dx_dy_costh_cosphi[i][j])).plus(cmyzn[ni-li][nj-lj][0][mi].times(dx_dy_costh_sinphi[i][j]))).times(exp_jk_rpr[i][j])));
                            Nphi[i][j].set(Nphi[i][j].plus((((cjxzn[ni-li][nj-lj][0][mi].times(-1.0)).times(dx_dy_sinphi[i][j])).plus(cjyzn[ni-li][nj-lj][0][mi].times(dx_dy_cosphi[i][j]))).times(exp_jk_rpr[i][j])));
                            Lphi[i][j].set(Lphi[i][j].plus((((cmxzn[ni-li][nj-lj][0][mi].times(-1.0)).times(dx_dy_sinphi[i][j])).plus(cmyzn[ni-li][nj-lj][0][mi].times(dx_dy_cosphi[i][j]))).times(exp_jk_rpr[i][j])));                        
                        }
                    }
                }
            }
            
            for (int i=0;i<thetaStepsCount;i++){
                for (int j = 0;j<phiStepsCount;j++){
                    farfield_dataTheta3D[mi][i][j] = (Math.pow(k, 2)/(8*Math.PI*Constants.ETA0*radiatedPower[mi]))*Math.pow((Lphi[i][j].plus(Ntheta[i][j].times(Constants.ETA0))).abs(),2);
                    farfield_dataPhi3D[mi][i][j] = (Math.pow(k, 2)/(8*Math.PI*Constants.ETA0*radiatedPower[mi]))*Math.pow((Ltheta[i][j].minus(Nphi[i][j].times(Constants.ETA0))).abs(),2);
                    farfield_dataThetaPhase3D[mi][i][j] = 2*(Lphi[i][j].plus(Ntheta[i][j].times(Constants.ETA0))).angle();
                    farfield_dataPhiPhase3D[mi][i][j] = 2*(Ltheta[i][j].minus(Nphi[i][j].times(Constants.ETA0))).angle();
                }
            }
        }
    }
        
    
    public void calFarField(){
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
        Common.save1DArray(radiatedPower, "RadiatedPower");
    }
    
    public void calTotalRadiatedPower(){
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

        this.omegas =  Common.genDouble1DArray(this.frequencies.length, 0.0);
        
	for (int i = 0; i < frequencies.length ; i++){
            omegas[i] = 2*Math.PI*frequencies[i];
	}
        
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
