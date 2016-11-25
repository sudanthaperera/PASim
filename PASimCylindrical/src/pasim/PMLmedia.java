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
public class PMLmedia {
    private Boundary b;
    private Cell c;
    private ProblemSpace ps;
    private EMobject EM;
    
    private double[] sigma_per_rn, sigma_pmr_rn, kappa_er_rn, kappa_mr_rn, alpha_er_rn, alpha_mr_rn, cpml_b_er_rn, cpml_a_er_rn,cpml_b_mr_rn, cpml_a_mr_rn;
    
    private double[][][] Psi_ear_rn, Psi_ezr_rn, Psi_har_rn, Psi_hzr_rn; 
    
    private double[][][] CPsi_ear_rn, CPsi_ezr_rn, CPsi_har_rn, CPsi_hzr_rn;
    
    private double[] sigma_per_rp, sigma_pmr_rp, kappa_er_rp, kappa_mr_rp, alpha_er_rp, alpha_mr_rp, cpml_b_er_rp, cpml_a_er_rp, cpml_b_mr_rp, cpml_a_mr_rp;
    
    private double[][][] Psi_ear_rp, Psi_ezr_rp, Psi_har_rp, Psi_hzr_rp;
    
    private double[][][] CPsi_ear_rp, CPsi_ezr_rp, CPsi_har_rp, CPsi_hzr_rp;
    
    private double[] sigma_pea_an, sigma_pma_an, kappa_ea_an, kappa_ma_an, alpha_ea_an, alpha_ma_an, cpml_b_ea_an, cpml_a_ea_an, cpml_b_ma_an, cpml_a_ma_an;
    
    private double[][][] Psi_eza_an, Psi_era_an, Psi_hza_an, Psi_hra_an; 
    
    private double[][][] CPsi_eza_an, CPsi_era_an, CPsi_hza_an, CPsi_hra_an;
    
    private double[] sigma_pea_ap, sigma_pma_ap, kappa_ea_ap, kappa_ma_ap, alpha_ea_ap, alpha_ma_ap, cpml_b_ea_ap, cpml_a_ea_ap, cpml_b_ma_ap, cpml_a_ma_ap;
    
    private double[][][] Psi_era_ap, Psi_eza_ap, Psi_hra_ap, Psi_hza_ap; 
    
    private double[][][] CPsi_eza_ap, CPsi_era_ap, CPsi_hza_ap, CPsi_hra_ap;
    
    private double[] sigma_pez_zn, sigma_pmz_zn, kappa_ez_zn, kappa_mz_zn, alpha_ez_zn, alpha_mz_zn, cpml_b_ez_zn, cpml_a_ez_zn, cpml_b_mz_zn, cpml_a_mz_zn;
    
    private double[][][] Psi_eaz_zn, Psi_erz_zn, Psi_haz_zn, Psi_hrz_zn; 
    
    private double[][][] CPsi_eaz_zn, CPsi_erz_zn, CPsi_haz_zn, CPsi_hrz_zn;
    
    private double[] sigma_pez_zp, sigma_pmz_zp, kappa_ez_zp, kappa_mz_zp, alpha_ez_zp, alpha_mz_zp, cpml_b_ez_zp, cpml_a_ez_zp, cpml_b_mz_zp, cpml_a_mz_zp;
    
    private double[][][] Psi_erz_zp, Psi_eaz_zp, Psi_hrz_zp, Psi_haz_zp; 
    
    private double[][][] CPsi_eaz_zp, CPsi_erz_zp, CPsi_haz_zp, CPsi_hrz_zp;
    
    public void saveAllArrayRn(){
        Common.save1DArray(sigma_per_rn, "sigma_per_rn");
        Common.save1DArray(sigma_pmr_rn, "sigma_pmr_rn");
        Common.save1DArray(kappa_er_rn, "kappa_er_rn");
        Common.save1DArray(kappa_mr_rn, "kappa_mr_rn");
        Common.save1DArray(alpha_er_rn, "alpha_er_rn");
        Common.save1DArray(alpha_mr_rn, "alpha_mr_rn");
        Common.save1DArray(cpml_b_er_rn, "cpml_b_er_rn");
        Common.save1DArray(cpml_a_er_rn, "cpml_a_er_rn");
        Common.save1DArray(cpml_b_mr_rn, "cpml_b_mr_rn");
        Common.save1DArray(cpml_a_mr_rn, "cpml_a_mr_rn");
        Common.save3DArray(Psi_ear_rn, "Psi_ear_rn");
        Common.save3DArray(Psi_ezr_rn, "Psi_ezr_rn");
        Common.save3DArray(Psi_har_rn, "Psi_har_rn");
        Common.save3DArray(Psi_hzr_rn, "Psi_hzr_rn");
        Common.save3DArray(CPsi_ear_rn, "CPsi_ear_rn");
        Common.save3DArray(CPsi_ezr_rn, "CPsi_ezr_rn");
        Common.save3DArray(CPsi_har_rn, "CPsi_har_rn");
        Common.save3DArray(CPsi_hzr_rn, "CPsi_hzr_rn");
    }
    
    public void saveAllArrayRp(){
        Common.save1DArray(sigma_per_rp, "sigma_per_rp");
        Common.save1DArray(sigma_pmr_rp, "sigma_pmr_rp");
        Common.save1DArray(kappa_er_rp, "kappa_er_rp");
        Common.save1DArray(kappa_mr_rp, "kappa_mr_rp");
        Common.save1DArray(alpha_er_rp, "alpha_er_rp");
        Common.save1DArray(alpha_mr_rp, "alpha_mr_rp");
        Common.save1DArray(cpml_b_er_rp, "cpml_b_er_rp");
        Common.save1DArray(cpml_a_er_rp, "cpml_a_er_rp");
        Common.save1DArray(cpml_b_mr_rp, "cpml_b_mr_rp");
        Common.save1DArray(cpml_a_mr_rp, "cpml_a_mr_rp");
        Common.save3DArray(Psi_ear_rp, "Psi_ear_rp");
        Common.save3DArray(Psi_ezr_rp, "Psi_ezr_rp");
        Common.save3DArray(Psi_har_rp, "Psi_har_rp");
        Common.save3DArray(Psi_hzr_rp, "Psi_hzr_rp");
        Common.save3DArray(CPsi_ear_rp, "CPsi_ear_rp");
        Common.save3DArray(CPsi_ezr_rp, "CPsi_ezr_rp");
        Common.save3DArray(CPsi_har_rp, "CPsi_har_rp");
        Common.save3DArray(CPsi_hzr_rp, "CPsi_hzr_rp");
    }

    public void saveAllArrayAn(){
        Common.save1DArray(sigma_pea_an, "sigma_pea_an");
        Common.save1DArray(sigma_pma_an, "sigma_pma_an");
        Common.save1DArray(kappa_ea_an, "kappa_ea_an");
        Common.save1DArray(kappa_ma_an, "kappa_ma_an");
        Common.save1DArray(alpha_ea_an, "alpha_ea_an");
        Common.save1DArray(alpha_ma_an, "alpha_ma_an");
        Common.save1DArray(cpml_b_ea_an, "cpml_b_ea_an");
        Common.save1DArray(cpml_a_ea_an, "cpml_a_ea_an");
        Common.save1DArray(cpml_b_ma_an, "cpml_b_ma_an");
        Common.save1DArray(cpml_a_ma_an, "cpml_a_ma_an");
        Common.save3DArray(Psi_era_an, "Psi_era_an");
        Common.save3DArray(Psi_eza_an, "Psi_eza_an");
        Common.save3DArray(Psi_hra_an, "Psi_hra_an");
        Common.save3DArray(Psi_hza_an, "Psi_hza_an");
        Common.save3DArray(CPsi_hra_an, "CPsi_hra_an");
        Common.save3DArray(CPsi_eza_an, "CPsi_eza_an");
        Common.save3DArray(CPsi_era_an, "CPsi_era_an");
        Common.save3DArray(CPsi_hza_an, "CPsi_hza_an");
    }
    
    public void saveAllArraaAp(){    
        Common.save1DArray(sigma_pea_ap, "sigma_pea_ap");
        Common.save1DArray(sigma_pma_ap, "sigma_pma_ap");
        Common.save1DArray(kappa_ea_ap, "kappa_ea_ap");
        Common.save1DArray(kappa_ma_ap, "kappa_ma_ap");
        Common.save1DArray(alpha_ea_ap, "alpha_ea_ap");
        Common.save1DArray(alpha_ma_ap, "alpha_ma_ap");
        Common.save1DArray(cpml_b_ea_ap, "cpml_b_ea_ap");
        Common.save1DArray(cpml_a_ea_ap, "cpml_a_ea_ap");
        Common.save1DArray(cpml_b_ma_ap, "cpml_b_ma_ap");
        Common.save1DArray(cpml_a_ma_ap, "cpml_a_ma_ap");
        Common.save3DArray(Psi_era_ap, "Psi_era_ap");
        Common.save3DArray(Psi_eza_ap, "Psi_eza_ap");
        Common.save3DArray(Psi_hra_ap, "Psi_hra_ap");
        Common.save3DArray(Psi_hza_ap, "Psi_hza_ap");
        Common.save3DArray(CPsi_eza_ap, "CPsi_eza_ap");
        Common.save3DArray(CPsi_era_ap, "CPsi_era_ap");
        Common.save3DArray(CPsi_hza_ap, "CPsi_hza_ap");
        Common.save3DArray(CPsi_hra_ap, "CPsi_hra_ap");
    }
    
    public void saveAllArrayZn(){
        Common.save1DArray(sigma_pez_zn, "sigma_pez_zn");
        Common.save1DArray(sigma_pmz_zn, "sigma_pmz_zn");
        Common.save1DArray(kappa_ez_zn, "kappa_ez_zn");
        Common.save1DArray(kappa_mz_zn, "kappa_mz_zn");
        Common.save1DArray(alpha_ez_zn, "alpha_ez_zn");
        Common.save1DArray(alpha_mz_zn, "alpha_mz_zn");
        Common.save1DArray(cpml_b_ez_zn, "cpml_b_ez_zn");
        Common.save1DArray(cpml_a_ez_zn, "cpml_a_ez_zn");
        Common.save1DArray(cpml_b_mz_zn, "cpml_b_mz_zn");
        Common.save1DArray(cpml_a_mz_zn, "cpml_a_mz_zn");
        Common.save3DArray(Psi_eaz_zn, "Psi_eaz_zn");
        Common.save3DArray(Psi_erz_zn, "Psi_erz_zn");
        Common.save3DArray(Psi_haz_zn, "Psi_haz_zn");
        Common.save3DArray(Psi_hrz_zn, "Psi_hrz_zn");
        Common.save3DArray(CPsi_eaz_zn, "CPsi_eaz_zn");
        Common.save3DArray(CPsi_erz_zn, "CPsi_erz_zn");
        Common.save3DArray(CPsi_haz_zn, "CPsi_haz_zn");
        Common.save3DArray(CPsi_hrz_zn, "CPsi_hrz_zn");
    }
    
    public void saveAllArrayZp(){
        Common.save1DArray(sigma_pez_zp, "sigma_pez_zp");
        Common.save1DArray(sigma_pmz_zp, "sigma_pmz_zp");
        Common.save1DArray(kappa_ez_zp, "kappa_ez_zp");
        Common.save1DArray(kappa_mz_zp, "kappa_mz_zp");
        Common.save1DArray(alpha_ez_zp, "alpha_ez_zp");
        Common.save1DArray(alpha_mz_zp, "alpha_mz_zp");
        Common.save1DArray(cpml_b_ez_zp, "cpml_b_ez_zp");
        Common.save1DArray(cpml_a_ez_zp, "cpml_a_ez_zp");
        Common.save1DArray(cpml_b_mz_zp, "cpml_b_mz_zp");
        Common.save1DArray(cpml_a_mz_zp, "cpml_a_mz_zp");
        Common.save3DArray(Psi_erz_zp, "Psi_erz_zp");
        Common.save3DArray(Psi_eaz_zp, "Psi_eaz_zp");
        Common.save3DArray(Psi_hrz_zp, "Psi_hrz_zp");
        Common.save3DArray(Psi_haz_zp, "Psi_haz_zp");
        Common.save3DArray(CPsi_eaz_zp, "CPsi_eaz_zp");
        Common.save3DArray(CPsi_erz_zp, "CPsi_erz_zp");
        Common.save3DArray(CPsi_haz_zp, "CPsi_haz_zp");
        Common.save3DArray(CPsi_hrz_zp, "CPsi_hrz_zp");
    }

    
    public PMLmedia(ProblemSpace ps,Boundary b, Cell c){
        this.ps = ps;
        this.b = b;
        this.c = c;
        this.EM = new EMobject();
    }
    
    public void initCPMLboundaryRn(){        
        //Initialize cpml for rn region
        if (b.getCPMLrn()){   
            int order = b.getOrder(); //order of the polynomial distribution
            double sigmaFactor = b.getSigmaFactor();
            double kappaMax = b.getKappaMax();
            double alphaMin = b.getAlphaMin();
            double alphaMax = b.getAlphaMax();
            double sigmaMax;
            int ncells;
            double[] rho_e_rn, rho_m_rn;
            
            //define one-dimensional temporary cpml parameter arrays 
            sigmaMax = sigmaFactor*(order+1)/(150*Math.PI*c.getDeltaR());
            ncells = b.getCellCountRn();
	    
            rho_e_rn = Common.genDouble1DArray(ncells,0.0);
            rho_m_rn = Common.genDouble1DArray(ncells,0.0);
            sigma_per_rn = Common.genDouble1DArray(ncells,0.0);
            sigma_pmr_rn = Common.genDouble1DArray(ncells,0.0);
            kappa_er_rn = Common.genDouble1DArray(ncells,0.0);
            kappa_mr_rn = Common.genDouble1DArray(ncells,0.0);
            alpha_er_rn = Common.genDouble1DArray(ncells,0.0);
            alpha_mr_rn = Common.genDouble1DArray(ncells,0.0);
            cpml_b_er_rn = Common.genDouble1DArray(ncells,0.0);
            cpml_a_er_rn = Common.genDouble1DArray(ncells,0.0);
            cpml_b_mr_rn = Common.genDouble1DArray(ncells,0.0);
            cpml_a_mr_rn = Common.genDouble1DArray(ncells,0.0);
	    
            for(int ncells_ind = 0; ncells_ind < ncells; ncells_ind++){		
                rho_e_rn[ncells_ind] = (((double)ncells - (double)ncells_ind) - 0.75)/((double)ncells);
                rho_m_rn[ncells_ind] = (((double)ncells - (double)ncells_ind) - 0.25)/((double)ncells);
		
                sigma_per_rn[ncells_ind] = 1;
                sigma_pmr_rn[ncells_ind] = 1;
		
		sigma_per_rn[ncells_ind] = sigmaMax*Math.pow(rho_e_rn[ncells_ind], order);
                sigma_pmr_rn[ncells_ind] = sigmaMax*Math.pow(rho_m_rn[ncells_ind], order);
		
                kappa_er_rn[ncells_ind] = 1 + (kappaMax - 1)*Math.pow(rho_e_rn[ncells_ind], order);
                kappa_mr_rn[ncells_ind] = 1 + (kappaMax - 1)*Math.pow(rho_m_rn[ncells_ind], order);
		
                sigma_pmr_rn[ncells_ind] = (Constants.MU0/Constants.EPS0)*sigma_pmr_rn[ncells_ind];
                alpha_er_rn[ncells_ind] = alphaMin + (alphaMax - alphaMin)*(1 - rho_e_rn[ncells_ind]);
                alpha_mr_rn[ncells_ind] = alphaMin + (alphaMax - alphaMin)*(1 - rho_m_rn[ncells_ind]);
                alpha_mr_rn[ncells_ind] = (Constants.MU0/Constants.EPS0)*alpha_mr_rn[ncells_ind];
		
                //define one-dimensional cpml parameter arrays 
                cpml_b_er_rn[ncells_ind] = Math.exp((-c.getDeltaT()/Constants.EPS0)*((sigma_per_rn[ncells_ind]/kappa_er_rn[ncells_ind]) + alpha_er_rn[ncells_ind])); 
                cpml_a_er_rn[ncells_ind] = (1/c.getDeltaR())*(cpml_b_er_rn[ncells_ind] - 1.0)* sigma_per_rn[ncells_ind]/(kappa_er_rn[ncells_ind]*(sigma_per_rn[ncells_ind] + kappa_er_rn[ncells_ind]*alpha_er_rn[ncells_ind]));
                cpml_b_mr_rn[ncells_ind] = Math.exp((-c.getDeltaT()/Constants.MU0)*((sigma_pmr_rn[ncells_ind]/kappa_mr_rn[ncells_ind]) + alpha_mr_rn[ncells_ind])); 
                cpml_a_mr_rn[ncells_ind] = (1/c.getDeltaR())*(cpml_b_mr_rn[ncells_ind] - 1.0)*sigma_pmr_rn[ncells_ind]/(kappa_mr_rn[ncells_ind]*(sigma_pmr_rn[ncells_ind] + kappa_mr_rn[ncells_ind]*alpha_mr_rn[ncells_ind]));
            }
	    
            //Create and initialize 2D cpml convolution parameters 
            Psi_ear_rn = Common.genDouble3DArray(ncells,ps.getNA(),ps.getNZ()+1,0.0); 
            Psi_ezr_rn = Common.genDouble3DArray(ncells,ps.getNA()+1,ps.getNZ(),0.0); 
            Psi_har_rn = Common.genDouble3DArray(ncells,ps.getNA()+1,ps.getNZ(),0.0); 
            Psi_hzr_rn = Common.genDouble3DArray(ncells,ps.getNA(),ps.getNZ()+1,0.0);
	    
            //Create and initialize 2D cpml convolution coefficients
	    
            CPsi_ear_rn = Common.genDouble3DArray(ncells,ps.getNA(),ps.getNZ()+1,0.0); 
            CPsi_hzr_rn = Common.genDouble3DArray(ncells,ps.getNA(),ps.getNZ()+1,0.0); 
            CPsi_ezr_rn = Common.genDouble3DArray(ncells,ps.getNA()+1,ps.getNZ(),0.0); 
            CPsi_har_rn = Common.genDouble3DArray(ncells,ps.getNA()+1,ps.getNZ(),0.0);
	    
            for(int i=0;i<ncells;i++){
                for(int j=0;j<ps.getNA()+1;j++){
                    for(int k=0;k<ps.getNZ()+1;k++){ 
                        if(j<ps.getNA()){	
                            CPsi_ear_rn[i][j][k] = EM.getCeahz(i+1,j,k)*c.getDeltaR();
                            EM.setCeahz(i+1, j, k, EM.getCeahz(i+1,j,k)/kappa_er_rn[i]);
                            CPsi_hzr_rn[i][j][k] = EM.getChzea(i,j,k)*c.getDeltaR();
                            EM.setChzea(i, j, k, EM.getChzea(i,j,k)/kappa_mr_rn[i]);
                        }
                        if(k<ps.getNZ()){
                            CPsi_ezr_rn[i][j][k] = EM.getCezha(i+1, j, k)*c.getDeltaR();
                            EM.setCezha(i+1, j, k, EM.getCezha(i+1,j,k)/kappa_er_rn[i]);
                            CPsi_har_rn[i][j][k] = EM.getChaez(i,j,k)*c.getDeltaR();
                            EM.setChaez(i, j, k, EM.getChaez(i, j, k)/kappa_mr_rn[i]);
                        }
                    }
                }
            }
        }
    }
    
    public void initCPMLboundaryRp(){ 
        //Initialize cpml for rp region
	if (b.getCPMLrp()){
            int order = b.getOrder(); //order of the polynomial distribution
            double sigmaFactor = b.getSigmaFactor();
            double kappaMax = b.getKappaMax();
            double alphaMin = b.getAlphaMin();
            double alphaMax = b.getAlphaMax();
            double sigmaMax ;
            int ncells;
            double[] rho_e_rp, rho_m_rp;
            
            //define one-dimensional temporary cpml parameter arrays 
            sigmaMax = sigmaFactor*(order + 1)/(150*Math.PI*c.getDeltaR());
            ncells = b.getCellCountRp();
	    
            rho_e_rp = Common.genDouble1DArray(ncells,0.0);
            rho_m_rp = Common.genDouble1DArray(ncells,0.0);
            sigma_per_rp = Common.genDouble1DArray(ncells,0.0);
            sigma_pmr_rp = Common.genDouble1DArray(ncells,0.0);
            kappa_er_rp = Common.genDouble1DArray(ncells,0.0);
            kappa_mr_rp = Common.genDouble1DArray(ncells,0.0);
            alpha_er_rp = Common.genDouble1DArray(ncells,0.0);
            alpha_mr_rp = Common.genDouble1DArray(ncells,0.0);
            cpml_b_er_rp = Common.genDouble1DArray(ncells,0.0);
            cpml_a_er_rp = Common.genDouble1DArray(ncells,0.0);
            cpml_b_mr_rp = Common.genDouble1DArray(ncells,0.0);
            cpml_a_mr_rp = Common.genDouble1DArray(ncells,0.0);
	    
            for(int ncells_ind = 0;ncells_ind<ncells;ncells_ind++){
		rho_e_rp[ncells_ind] = (((double)ncells_ind + 1) - 0.75)/((double)ncells);
		rho_m_rp[ncells_ind] = (((double)ncells_ind + 1) - 0.25)/((double)ncells);
		
		sigma_per_rp[ncells_ind] = 1;
		sigma_pmr_rp[ncells_ind] = 1;
		
                sigma_per_rp[ncells_ind] = sigmaMax*Math.pow(rho_e_rp[ncells_ind], order);
                sigma_pmr_rp[ncells_ind] = sigmaMax*Math.pow(rho_m_rp[ncells_ind], order);
		
		kappa_er_rp[ncells_ind] = 1 + (kappaMax - 1)*Math.pow(rho_e_rp[ncells_ind], order);
		kappa_mr_rp[ncells_ind] = 1 + (kappaMax - 1)*Math.pow(rho_m_rp[ncells_ind], order);
		
		sigma_pmr_rp[ncells_ind] = (Constants.MU0/Constants.EPS0)*sigma_pmr_rp[ncells_ind];
		alpha_er_rp[ncells_ind] = alphaMin + (alphaMax - alphaMin)*(1 - rho_e_rp[ncells_ind]);
		alpha_mr_rp[ncells_ind] = alphaMin + (alphaMax - alphaMin)*(1 - rho_m_rp[ncells_ind]);
		alpha_mr_rp[ncells_ind] = (Constants.MU0/Constants.EPS0)*alpha_mr_rp[ncells_ind];
		
		//define one-dimensional cpml parameter arrays 
		cpml_b_er_rp[ncells_ind] = Math.exp((-c.getDeltaT()/Constants.EPS0)*((sigma_per_rp[ncells_ind]/kappa_er_rp[ncells_ind]) + alpha_er_rp[ncells_ind]));
		cpml_a_er_rp[ncells_ind] = (1/c.getDeltaR())*(cpml_b_er_rp[ncells_ind] - 1.0)* sigma_per_rp[ncells_ind]/(kappa_er_rp[ncells_ind]*(sigma_per_rp[ncells_ind] + kappa_er_rp[ncells_ind]*alpha_er_rp[ncells_ind]));
		cpml_b_mr_rp[ncells_ind] = Math.exp((-c.getDeltaT()/Constants.MU0)*((sigma_pmr_rp[ncells_ind]/kappa_mr_rp[ncells_ind]) + alpha_mr_rp[ncells_ind]));
		cpml_a_mr_rp[ncells_ind] = (1/c.getDeltaR())*(cpml_b_mr_rp[ncells_ind] - 1.0)*sigma_pmr_rp[ncells_ind]/(kappa_mr_rp[ncells_ind]*(sigma_pmr_rp[ncells_ind] + kappa_mr_rp[ncells_ind]*alpha_mr_rp[ncells_ind]));
            }
	    
            //Create and initialize 2D cpml convolution parameters 
            Psi_ear_rp = Common.genDouble3DArray(ncells,ps.getNA(),ps.getNZ()+1,0.0); 
            Psi_ezr_rp = Common.genDouble3DArray(ncells,ps.getNA()+1,ps.getNZ(),0.0); 
            Psi_har_rp = Common.genDouble3DArray(ncells,ps.getNA()+1,ps.getNZ(),0.0); 
            Psi_hzr_rp = Common.genDouble3DArray(ncells,ps.getNA(),ps.getNZ()+1,0.0); 
	    
            //Create and initialize 2D cpml convolution coefficients
	    
            CPsi_ear_rp = Common.genDouble3DArray(ncells,ps.getNA(),ps.getNZ()+1,0.0); 
            CPsi_hzr_rp = Common.genDouble3DArray(ncells,ps.getNA(),ps.getNZ()+1,0.0); 
            CPsi_ezr_rp = Common.genDouble3DArray(ncells,ps.getNA()+1,ps.getNZ(),0.0); 
            CPsi_har_rp = Common.genDouble3DArray(ncells,ps.getNA()+1,ps.getNZ(),0.0); 
	    
            for(int i=0;i<ncells;i++){
		for(int j=0;j<ps.getNA()+1;j++){
                    for(int k=0;k<ps.getNZ()+1;k++){ 
			if(j<ps.getNA()){	
                            CPsi_ear_rp[i][j][k] = EM.getCeahz(ps.getNR() - ncells + i, j, k)*c.getDeltaR();
                            EM.setCeahz(ps.getNR() - ncells + i, j, k, EM.getCeahz(ps.getNR() - ncells + i, j, k)/kappa_er_rp[i]);
                            CPsi_hzr_rp[i][j][k] = EM.getChzea(ps.getNR() - ncells + i, j, k)*c.getDeltaR();
                            EM.setChzea(ps.getNR() - ncells + i, j, k, EM.getChzea(ps.getNR() - ncells + i, j, k)/kappa_mr_rp[i]);
			}
			if(k<ps.getNZ()){
                            CPsi_ezr_rp[i][j][k] = EM.getCezha(ps.getNR() - ncells + i, j, k)*c.getDeltaR();
                            EM.setCezha(ps.getNR() - ncells + i, j, k, EM.getCezha(ps.getNR() - ncells + i, j, k)/kappa_er_rp[i]);
                            CPsi_har_rp[i][j][k] = EM.getChaez(ps.getNR() - ncells + i, j, k)*c.getDeltaR();
                            EM.setChaez(ps.getNR() - ncells + i, j, k, EM.getChaez(ps.getNR() - ncells + i, j, k)/kappa_mr_rp[i]);
                        }
                    }
		}
            }
  	}
    }
    
    public void initCPMLboundaryAn(){
	//Initialize cpml for yn region
	if (b.getCPMLan()){
            int order = b.getOrder(); //order of the polynomial distribution
            double sigmaFactor = b.getSigmaFactor();
            double kappaMax = b.getKappaMax();
            double alphaMin = b.getAlphaMin();
            double alphaMax = b.getAlphaMax();
            double sigmaMax ;
            double[] rho_e_an, rho_m_an;
            int ncells;
            //define one-dimensional temporary cpml parameter arrays 
            sigmaMax = sigmaFactor*(order+1)/(150*Math.PI*c.getDeltaA());
            ncells = b.getCellCountAn();
	    
            rho_e_an = Common.genDouble1DArray(ncells,0.0);
            rho_m_an = Common.genDouble1DArray(ncells,0.0);
            sigma_pea_an = Common.genDouble1DArray(ncells,0.0);
            sigma_pma_an = Common.genDouble1DArray(ncells,0.0);
            kappa_ea_an = Common.genDouble1DArray(ncells,0.0);
            kappa_ma_an = Common.genDouble1DArray(ncells,0.0);
            alpha_ea_an = Common.genDouble1DArray(ncells,0.0);
            alpha_ma_an = Common.genDouble1DArray(ncells,0.0);
            cpml_b_ea_an = Common.genDouble1DArray(ncells,0.0);
            cpml_a_ea_an = Common.genDouble1DArray(ncells,0.0);
            cpml_b_ma_an = Common.genDouble1DArray(ncells,0.0);
            cpml_a_ma_an = Common.genDouble1DArray(ncells,0.0);
	    
            for(int ncells_ind = 0;ncells_ind<ncells;ncells_ind++){
		rho_e_an[ncells_ind] = (((double)ncells - (double)ncells_ind) - 0.75)/((double)ncells);
		rho_m_an[ncells_ind] = (((double)ncells - (double)ncells_ind) - 0.25)/((double)ncells);
		
		sigma_pea_an[ncells_ind] = 1;
		sigma_pma_an[ncells_ind] = 1;
		
                sigma_pea_an[ncells_ind] = sigmaMax*Math.pow(rho_e_an[ncells_ind],order);
                sigma_pma_an[ncells_ind] = sigmaMax*Math.pow(rho_m_an[ncells_ind],order);
		
		kappa_ea_an[ncells_ind] = 1 + (kappaMax - 1)*Math.pow(rho_e_an[ncells_ind],order);
		kappa_ma_an[ncells_ind] = 1 + (kappaMax - 1)*Math.pow(rho_m_an[ncells_ind],order);
		
		sigma_pma_an[ncells_ind] = (Constants.MU0 / Constants.EPS0)*sigma_pma_an[ncells_ind];
		alpha_ea_an[ncells_ind] = alphaMin + (alphaMax - alphaMin)*(1 - rho_e_an[ncells_ind]);
		alpha_ma_an[ncells_ind] = alphaMin + (alphaMax - alphaMin)*(1 - rho_m_an[ncells_ind]);
		alpha_ma_an[ncells_ind] = (Constants.MU0/Constants.EPS0)*alpha_ma_an[ncells_ind];
		
		//define one-dimensional cpml parameter arrays 
		cpml_b_ea_an[ncells_ind] = Math.exp((-c.getDeltaT()/Constants.EPS0)*((sigma_pea_an[ncells_ind]/kappa_ea_an[ncells_ind]) + alpha_ea_an[ncells_ind])); 
		cpml_a_ea_an[ncells_ind] = (1/c.getDeltaA())*(cpml_b_ea_an[ncells_ind] - 1.0)* sigma_pea_an[ncells_ind]/(kappa_ea_an[ncells_ind]*(sigma_pea_an[ncells_ind] + kappa_ea_an[ncells_ind]*alpha_ea_an[ncells_ind]));
		cpml_b_ma_an[ncells_ind] = Math.exp((-c.getDeltaT()/Constants.MU0)*((sigma_pma_an[ncells_ind]/kappa_ma_an[ncells_ind]) + alpha_ma_an[ncells_ind])); 
		cpml_a_ma_an[ncells_ind] = (1/c.getDeltaA())*(cpml_b_ma_an[ncells_ind] - 1.0)*sigma_pma_an[ncells_ind]/(kappa_ma_an[ncells_ind]*(sigma_pma_an[ncells_ind] + kappa_ma_an[ncells_ind]*alpha_ma_an[ncells_ind]));
            }
	    
            //Create and initialize 2D cpml convolution parameters 
            Psi_eza_an = Common.genDouble3DArray(ps.getNR()+1,ncells,ps.getNZ(),0.0); 
            Psi_era_an = Common.genDouble3DArray(ps.getNR(),ncells,ps.getNZ()+1,0.0); 
            Psi_hza_an = Common.genDouble3DArray(ps.getNR(),ncells,ps.getNZ()+1,0.0); 
            Psi_hra_an = Common.genDouble3DArray(ps.getNR()+1,ncells,ps.getNZ(),0.0); 
	    
            //Create and initialize 2D cpml convolution coefficients
            CPsi_eza_an = Common.genDouble3DArray(ps.getNR()+1,ncells,ps.getNZ(),0.0); 
            CPsi_hra_an = Common.genDouble3DArray(ps.getNR()+1,ncells,ps.getNZ(),0.0); 
            CPsi_era_an = Common.genDouble3DArray(ps.getNR(),ncells,ps.getNZ()+1,0.0); 
            CPsi_hza_an = Common.genDouble3DArray(ps.getNR(),ncells,ps.getNZ()+1,0.0);
	    
            for(int i=0;i<ps.getNR()+1;i++){
		for(int j=0;j<ncells;j++){
		    for(int k=0;k<ps.getNZ()+1;k++){ 
			if(i<ps.getNR()){
			    CPsi_era_an[i][j][k] = EM.getCerhz(i, j+1, k)*c.getDeltaA();
			    EM.setCerhz(i, j+1, k, EM.getCerhz(i, j+1, k)/kappa_ea_an[j]);
			    CPsi_hza_an[i][j][k] = EM.getChzer(i, j, k)*c.getDeltaA();
			    EM.setChzer(i, j, k, EM.getChzer(i, j, k)/kappa_ma_an[j]);
			    
			}
			if(k<ps.getNZ()){
			    CPsi_hra_an[i][j][k] = EM.getChrez(i, j, k)*c.getDeltaA();
			    EM.setChrez(i,j,k, EM.getChrez(i, j, k)/kappa_ma_an[j]);
			    CPsi_eza_an[i][j][k] = EM.getCezhr(i, j+1, k)*c.getDeltaA();
			    EM.setCezhr(i, j+1, k, EM.getCezhr(i, j+1, k)/kappa_ea_an[j]);
			    
			}
		    }
		}
	    }
	}
    } 
    
    public void initCPMLboundaryAp(){
	//Initialize cpml for yp region
	if (b.getCPMLap()){
	    int order = b.getOrder(); //order of the polynomial distribution
	    double sigmaFactor = b.getSigmaFactor();
	    double kappaMax = b.getKappaMax();
	    double alphaMin = b.getAlphaMin();
	    double alphaMax = b.getAlphaMax();
	    double sigmaMax ;
	    int ncells;
	    double[] rho_e_ap, rho_m_ap;
	    
	    //define one-dimensional temporary cpml parameter arrays 
	    sigmaMax = sigmaFactor*(order + 1)/(150*Math.PI*c.getDeltaA());
	    ncells = b.getCellCountAp();
	    
	    rho_e_ap = Common.genDouble1DArray(ncells,0.0);
	    rho_m_ap = Common.genDouble1DArray(ncells,0.0);
	    sigma_pea_ap = Common.genDouble1DArray(ncells,0.0);
	    sigma_pma_ap = Common.genDouble1DArray(ncells,0.0);
	    kappa_ea_ap = Common.genDouble1DArray(ncells,0.0);
	    kappa_ma_ap = Common.genDouble1DArray(ncells,0.0);
	    alpha_ea_ap = Common.genDouble1DArray(ncells,0.0);
	    alpha_ma_ap = Common.genDouble1DArray(ncells,0.0);
	    cpml_b_ea_ap = Common.genDouble1DArray(ncells,0.0);
	    cpml_a_ea_ap = Common.genDouble1DArray(ncells,0.0);
	    cpml_b_ma_ap = Common.genDouble1DArray(ncells,0.0);
	    cpml_a_ma_ap = Common.genDouble1DArray(ncells,0.0);
	    
	    for(int ncells_ind = 0;ncells_ind<ncells;ncells_ind++){
		rho_e_ap[ncells_ind] = (((double)ncells_ind + 1) - 0.75)/((double)ncells);
		rho_m_ap[ncells_ind] = (((double)ncells_ind + 1) - 0.25)/((double)ncells);
		
		sigma_pea_ap[ncells_ind] = 1;
		sigma_pma_ap[ncells_ind] = 1;
		
		sigma_pea_ap[ncells_ind] = sigmaMax*Math.pow(rho_e_ap[ncells_ind],order);
		sigma_pma_ap[ncells_ind] = sigmaMax*Math.pow(rho_m_ap[ncells_ind],order);
		
		kappa_ea_ap[ncells_ind] = 1 + (kappaMax - 1)*Math.pow(rho_e_ap[ncells_ind],order);
		kappa_ma_ap[ncells_ind] = 1 + (kappaMax - 1)*Math.pow(rho_m_ap[ncells_ind],order);
		
		sigma_pma_ap[ncells_ind] = (Constants.MU0/Constants.EPS0)*sigma_pma_ap[ncells_ind];
		alpha_ea_ap[ncells_ind] = alphaMin + (alphaMax - alphaMin)*(1 - rho_e_ap[ncells_ind]);
		alpha_ma_ap[ncells_ind] = alphaMin + (alphaMax - alphaMin)*(1 - rho_m_ap[ncells_ind]);
		alpha_ma_ap[ncells_ind] = (Constants.MU0/Constants.EPS0)*alpha_ma_ap[ncells_ind];
		
		//define one-dimensional cpml parameter arrays 
		cpml_b_ea_ap[ncells_ind] = Math.exp((-c.getDeltaT()/Constants.EPS0)*((sigma_pea_ap[ncells_ind]/kappa_ea_ap[ncells_ind]) + alpha_ea_ap[ncells_ind]));
		cpml_a_ea_ap[ncells_ind] = (1/c.getDeltaA())*(cpml_b_ea_ap[ncells_ind] - 1.0)* sigma_pea_ap[ncells_ind]/(kappa_ea_ap[ncells_ind]*(sigma_pea_ap[ncells_ind] + kappa_ea_ap[ncells_ind]*alpha_ea_ap[ncells_ind]));
		cpml_b_ma_ap[ncells_ind] = Math.exp((-c.getDeltaT()/Constants.MU0)*((sigma_pma_ap[ncells_ind]/kappa_ma_ap[ncells_ind]) + alpha_ma_ap[ncells_ind]));
		cpml_a_ma_ap[ncells_ind] = (1/c.getDeltaA())*(cpml_b_ma_ap[ncells_ind] - 1.0)*sigma_pma_ap[ncells_ind]/(kappa_ma_ap[ncells_ind]*(sigma_pma_ap[ncells_ind] + kappa_ma_ap[ncells_ind]*alpha_ma_ap[ncells_ind]));
	    }
	    
	    //Create and initialize 2D cpml convolution parameters 
	    Psi_eza_ap = Common.genDouble3DArray(ps.getNR()+1,ncells,ps.getNZ(),0.0); 
	    Psi_era_ap = Common.genDouble3DArray(ps.getNR(),ncells,ps.getNZ()+1,0.0); 
	    Psi_hza_ap = Common.genDouble3DArray(ps.getNR(),ncells,ps.getNZ()+1,0.0); 
	    Psi_hra_ap = Common.genDouble3DArray(ps.getNR()+1,ncells,ps.getNZ(),0.0); 
	    
	    //Create and initialize 2D cpml convolution coefficients
	    
	    CPsi_era_ap = Common.genDouble3DArray(ps.getNR(),ncells,ps.getNZ()+1,0.0); 
	    CPsi_hza_ap = Common.genDouble3DArray(ps.getNR(),ncells,ps.getNZ()+1,0.0); 
	    CPsi_eza_ap = Common.genDouble3DArray(ps.getNR()+1,ncells,ps.getNZ(),0.0); 
	    CPsi_hra_ap = Common.genDouble3DArray(ps.getNR()+1,ncells,ps.getNZ(),0.0); 
	    
	    for(int i=0;i<ps.getNR()+1;i++){
		for(int j=0;j<ncells;j++){
		    for(int k=0;k<ps.getNZ()+1;k++){ 
			if(i<ps.getNR()){	
			    CPsi_era_ap[i][j][k] = EM.getCerhz(i, ps.getNA() - ncells + j, k)*c.getDeltaA();
			    EM.setCerhz(i, ps.getNA() - ncells + j, k, EM.getCerhz(i, ps.getNA() - ncells + j, k)/kappa_ea_ap[j]);
			    CPsi_hza_ap[i][j][k] = EM.getChzer(i, ps.getNA() - ncells + j, k)*c.getDeltaA();
			    EM.setChzer(i, ps.getNA() - ncells + j, k, EM.getChzer(i, ps.getNA() - ncells + j, k)/kappa_ma_ap[j]);
			}
			if(k<ps.getNZ()){
			    CPsi_eza_ap[i][j][k] = EM.getCezhr(i, ps.getNA() - ncells + j, k)*c.getDeltaA();
			    EM.setCezhr(i, ps.getNA() - ncells + j, k, EM.getCezhr(i, ps.getNA() - ncells + j, k)/kappa_ea_ap[j]);
			    CPsi_hra_ap[i][j][k] = EM.getChrez(i, ps.getNA() - ncells + j, k)*c.getDeltaA();
			    EM.setChrez(i, ps.getNA() - ncells + j, k, EM.getChrez(i, ps.getNA() - ncells + j, k)/kappa_ma_ap[j]);
			}
		    }
		}
	    }
	}        
    } 
    
    public void initCPMLboundaryZn(){ 
	//Initialize cpml for zn region
	if (b.getCPMLzn()){
	    int order = b.getOrder(); //order of the polynomial distribution
	    double sigmaFactor = b.getSigmaFactor();
	    double kappaMax = b.getKappaMax();
	    double alphaMin = b.getAlphaMin();
	    double alphaMax = b.getAlphaMax();
	    double sigmaMax ;
	    int ncells;
	    double[] rho_e_zn, rho_m_zn;
	    
	    //define one-dimensional temporary cpml parameter arrays 
	    sigmaMax = sigmaFactor*(order+1)/(150*Math.PI*c.getDeltaZ());
	    ncells = b.getCellCountZn();
	    
	    rho_e_zn = Common.genDouble1DArray(ncells,0.0);
	    rho_m_zn = Common.genDouble1DArray(ncells,0.0);
	    sigma_pez_zn = Common.genDouble1DArray(ncells,0.0);
	    sigma_pmz_zn = Common.genDouble1DArray(ncells,0.0);
	    kappa_ez_zn = Common.genDouble1DArray(ncells,0.0);
	    kappa_mz_zn = Common.genDouble1DArray(ncells,0.0);
	    alpha_ez_zn = Common.genDouble1DArray(ncells,0.0);
	    alpha_mz_zn = Common.genDouble1DArray(ncells,0.0);
	    cpml_b_ez_zn = Common.genDouble1DArray(ncells,0.0);
	    cpml_a_ez_zn = Common.genDouble1DArray(ncells,0.0);
	    cpml_b_mz_zn = Common.genDouble1DArray(ncells,0.0);
	    cpml_a_mz_zn = Common.genDouble1DArray(ncells,0.0);
	    
	    for(int ncells_ind = 0;ncells_ind<ncells;ncells_ind++){
		rho_e_zn[ncells_ind] = (((double)ncells - (double)ncells_ind) - 0.75)/((double)ncells);
		rho_m_zn[ncells_ind] = (((double)ncells - (double)ncells_ind) - 0.25)/((double)ncells);
		
		sigma_pez_zn[ncells_ind] = 1;
		sigma_pmz_zn[ncells_ind] = 1;
		
		sigma_pez_zn[ncells_ind] = sigmaMax*Math.pow(rho_e_zn[ncells_ind],order);
		sigma_pmz_zn[ncells_ind] = sigmaMax*Math.pow(rho_m_zn[ncells_ind],order);
		
		kappa_ez_zn[ncells_ind] = 1 + (kappaMax - 1)*Math.pow(rho_e_zn[ncells_ind],order);
		kappa_mz_zn[ncells_ind] = 1 + (kappaMax - 1)*Math.pow(rho_m_zn[ncells_ind],order);
		
		sigma_pmz_zn[ncells_ind] = (Constants.MU0 / Constants.EPS0)*sigma_pmz_zn[ncells_ind];
		alpha_ez_zn[ncells_ind] = alphaMin + (alphaMax - alphaMin)*(1 - rho_e_zn[ncells_ind]);
		alpha_mz_zn[ncells_ind] = alphaMin + (alphaMax - alphaMin)*(1 - rho_m_zn[ncells_ind]);
		alpha_mz_zn[ncells_ind] = (Constants.MU0/Constants.EPS0)*alpha_mz_zn[ncells_ind];
		
		//define one-dimensional cpml parameter arrays 
		cpml_b_ez_zn[ncells_ind] = Math.exp((-c.getDeltaT()/Constants.EPS0)*((sigma_pez_zn[ncells_ind]/kappa_ez_zn[ncells_ind]) + alpha_ez_zn[ncells_ind])); 
		cpml_a_ez_zn[ncells_ind] = (1/c.getDeltaZ())*(cpml_b_ez_zn[ncells_ind] - 1.0)* sigma_pez_zn[ncells_ind]/(kappa_ez_zn[ncells_ind]*(sigma_pez_zn[ncells_ind] + kappa_ez_zn[ncells_ind]*alpha_ez_zn[ncells_ind]));
		cpml_b_mz_zn[ncells_ind] = Math.exp((-c.getDeltaT()/Constants.MU0)*((sigma_pmz_zn[ncells_ind]/kappa_mz_zn[ncells_ind]) + alpha_mz_zn[ncells_ind])); 
		cpml_a_mz_zn[ncells_ind] = (1/c.getDeltaZ())*(cpml_b_mz_zn[ncells_ind] - 1.0)*sigma_pmz_zn[ncells_ind]/(kappa_mz_zn[ncells_ind]*(sigma_pmz_zn[ncells_ind] + kappa_mz_zn[ncells_ind]*alpha_mz_zn[ncells_ind]));
	    }
	    //Create and initialize 2D cpml convolution parameters 
	    Psi_eaz_zn = Common.genDouble3DArray(ps.getNR()+1,ps.getNA(),ncells,0.0); 
	    Psi_erz_zn = Common.genDouble3DArray(ps.getNR(),ps.getNA()+1,ncells,0.0); 
	    Psi_haz_zn = Common.genDouble3DArray(ps.getNR(),ps.getNA()+1,ncells,0.0); 
	    Psi_hrz_zn = Common.genDouble3DArray(ps.getNR()+1,ps.getNA(),ncells,0.0); 
	    
	    //Create and initialize 2D cpml convolution coefficients
	    
	    CPsi_eaz_zn = Common.genDouble3DArray(ps.getNR()+1,ps.getNA(),ncells,0.0); 
	    CPsi_hrz_zn = Common.genDouble3DArray(ps.getNR()+1,ps.getNA(),ncells,0.0); 
	    CPsi_erz_zn = Common.genDouble3DArray(ps.getNR(),ps.getNA()+1,ncells,0.0); 
	    CPsi_haz_zn = Common.genDouble3DArray(ps.getNR(),ps.getNA()+1,ncells,0.0); 
	    
	    for(int i=0;i<ps.getNR()+1;i++){
		for(int j=0;j<ps.getNA()+1;j++){
		    for(int k=0;k<ncells;k++){ 
			if(j<ps.getNA()){	
			    CPsi_eaz_zn[i][j][k] = EM.getCeahr(i, j, k+1)*c.getDeltaZ();
			    EM.setCeahr(i, j, k+1, EM.getCeahr(i, j, k+1)/kappa_ez_zn[k]);
			    CPsi_hrz_zn[i][j][k] = EM.getChrea(i, j, k)*c.getDeltaZ();
			    EM.setChrea(i, j, k, EM.getChrea(i, j, k)/kappa_mz_zn[k]);
			}
			if(i<ps.getNR()){
			    CPsi_erz_zn[i][j][k] = EM.getCerha(i, j, k+1)*c.getDeltaZ();
			    EM.setCerha(i, j, k+1, EM.getCerha(i, j, k+1)/kappa_ez_zn[k]);
			    CPsi_haz_zn[i][j][k] = EM.getChaer(i, j, k)*c.getDeltaZ();
			    EM.setChaer(i, j, k,  EM.getChaer(i, j, k)/kappa_mz_zn[k]);
			}
		    }
		}
	    }
	}
        
    }
    
    public void initCPMLboundaryZp(){ 
	//Initialize cpml for zp region
	if (b.getCPMLzp()){
	    int order = b.getOrder(); //order of the polynomial distribution
	    double sigmaFactor = b.getSigmaFactor();
	    double kappaMax = b.getKappaMax();
	    double alphaMin = b.getAlphaMin();
	    double alphaMax = b.getAlphaMax();
	    double sigmaMax;
	    int ncells;
	    double[] rho_e_zp, rho_m_zp;
	    
	    //define one-dimensional temporary cpml parameter arrays 
	    sigmaMax = sigmaFactor*(order + 1)/(150*Math.PI*c.getDeltaZ());
	    ncells = b.getCellCountZp();
	    
	    
	    rho_e_zp = Common.genDouble1DArray(ncells,0.0);
	    rho_m_zp = Common.genDouble1DArray(ncells,0.0);
	    sigma_pez_zp = Common.genDouble1DArray(ncells,0.0);
	    sigma_pmz_zp = Common.genDouble1DArray(ncells,0.0);
	    kappa_ez_zp = Common.genDouble1DArray(ncells,0.0);
	    kappa_mz_zp = Common.genDouble1DArray(ncells,0.0);
	    alpha_ez_zp = Common.genDouble1DArray(ncells,0.0);
	    alpha_mz_zp = Common.genDouble1DArray(ncells,0.0);
	    cpml_b_ez_zp = Common.genDouble1DArray(ncells,0.0);
	    cpml_a_ez_zp = Common.genDouble1DArray(ncells,0.0);
	    cpml_b_mz_zp = Common.genDouble1DArray(ncells,0.0);
	    cpml_a_mz_zp = Common.genDouble1DArray(ncells,0.0);
	    
	    for(int ncells_ind = 0;ncells_ind<ncells;ncells_ind++){
		rho_e_zp[ncells_ind] = (((double)ncells_ind + 1) - 0.75)/((double)ncells);
		rho_m_zp[ncells_ind] = (((double)ncells_ind + 1) - 0.25)/((double)ncells);
		
		sigma_pez_zp[ncells_ind] = 1;
		sigma_pmz_zp[ncells_ind] = 1;
		
		sigma_pez_zp[ncells_ind] = sigmaMax*Math.pow(rho_e_zp[ncells_ind],order);
		sigma_pmz_zp[ncells_ind] = sigmaMax*Math.pow(rho_m_zp[ncells_ind],order);
		
		kappa_ez_zp[ncells_ind] = 1 + (kappaMax - 1)*Math.pow(rho_e_zp[ncells_ind],order);
		kappa_mz_zp[ncells_ind] = 1 + (kappaMax - 1)*Math.pow(rho_m_zp[ncells_ind],order);
		
		sigma_pmz_zp[ncells_ind] = (Constants.MU0 / Constants.EPS0)*sigma_pmz_zp[ncells_ind];
		alpha_ez_zp[ncells_ind] = alphaMin + (alphaMax - alphaMin)*(1 - rho_e_zp[ncells_ind]);
		alpha_mz_zp[ncells_ind] = alphaMin + (alphaMax - alphaMin)*(1 - rho_m_zp[ncells_ind]);
		alpha_mz_zp[ncells_ind] = (Constants.MU0/Constants.EPS0)*alpha_mz_zp[ncells_ind];
		
		//define one-dimensional cpml parameter arrays 
		cpml_b_ez_zp[ncells_ind] = Math.exp((-c.getDeltaT()/Constants.EPS0)*((sigma_pez_zp[ncells_ind]/kappa_ez_zp[ncells_ind]) + alpha_ez_zp[ncells_ind]));
		cpml_a_ez_zp[ncells_ind] = (1/c.getDeltaZ())*(cpml_b_ez_zp[ncells_ind] - 1.0)* sigma_pez_zp[ncells_ind]/(kappa_ez_zp[ncells_ind]*(sigma_pez_zp[ncells_ind] + kappa_ez_zp[ncells_ind]*alpha_ez_zp[ncells_ind]));
		cpml_b_mz_zp[ncells_ind] = Math.exp((-c.getDeltaT()/Constants.MU0)*((sigma_pmz_zp[ncells_ind]/kappa_mz_zp[ncells_ind]) + alpha_mz_zp[ncells_ind]));
		cpml_a_mz_zp[ncells_ind] = (1/c.getDeltaZ())*(cpml_b_mz_zp[ncells_ind] - 1.0)*sigma_pmz_zp[ncells_ind]/(kappa_mz_zp[ncells_ind]*(sigma_pmz_zp[ncells_ind] + kappa_mz_zp[ncells_ind]*alpha_mz_zp[ncells_ind]));
	    }
	    
	    //Create and initialize 2D cpml convolution parameters 
	    Psi_erz_zp = Common.genDouble3DArray(ps.getNR(),ps.getNA()+1,ncells,0.0); 
	    Psi_eaz_zp = Common.genDouble3DArray(ps.getNR()+1,ps.getNA(),ncells,0.0); 
	    Psi_hrz_zp = Common.genDouble3DArray(ps.getNR()+1,ps.getNA(),ncells,0.0); 
	    Psi_haz_zp = Common.genDouble3DArray(ps.getNR(),ps.getNA()+1,ncells,0.0); 
	    
	    //Create and initialize 2D cpml convolution coefficients
	    
	    CPsi_eaz_zp = Common.genDouble3DArray(ps.getNR()+1,ps.getNA(),ncells,0.0); 
	    CPsi_hrz_zp = Common.genDouble3DArray(ps.getNR()+1,ps.getNA(),ncells,0.0); 
	    CPsi_erz_zp = Common.genDouble3DArray(ps.getNR(),ps.getNA()+1,ncells,0.0); 
	    CPsi_haz_zp = Common.genDouble3DArray(ps.getNR(),ps.getNA()+1,ncells,0.0); 
	    
	    for(int i=0;i<ps.getNR()+1;i++){
		for(int j=0;j<ps.getNA()+1;j++){
		    for(int k=0;k<ncells;k++){ 
			if(j<ps.getNA()){	
			    CPsi_eaz_zp[i][j][k] = EM.getCeahr(i, j, ps.getNZ() - ncells + k)*c.getDeltaZ();
			    EM.setCeahr(i, j, ps.getNZ() - ncells + k, EM.getCeahr(i, j, ps.getNZ() - ncells + k)/kappa_ez_zp[k]);
			    CPsi_hrz_zp[i][j][k] = EM.getChrea(i, j, ps.getNZ() - ncells + k)*c.getDeltaZ();
			    EM.setChrea(i, j, ps.getNZ() - ncells + k, EM.getChrea(i, j, ps.getNZ() - ncells + k)/kappa_mz_zp[k]);
			}
			if(i<ps.getNR()){
			    CPsi_erz_zp[i][j][k] = EM.getCerha(i, j, ps.getNZ() - ncells + k)*c.getDeltaZ();
			    EM.setCerha(i, j, ps.getNZ() - ncells + k, EM.getCerha(i, j, ps.getNZ() - ncells + k)/kappa_ez_zp[k]);
			    CPsi_haz_zp[i][j][k] = EM.getChaer(i, j, ps.getNZ() - ncells + k)*c.getDeltaZ();
			    EM.setChaer(i, j, ps.getNZ() - ncells + k, EM.getChaer(i, j, ps.getNZ() - ncells + k)/kappa_mz_zp[k]);
			}
		    }
		}
	    }
	}        
    }
    
    public void initAllCPMLboundary(){
	if(b.isAnySideCPML()) {
            //Initialize CPML boundary condition
            initCPMLboundaryRn();
            initCPMLboundaryRp();
            initCPMLboundaryAn();
            initCPMLboundaryAp();
            initCPMLboundaryZn();
            initCPMLboundaryZp();
	}        
    }
    
    public void applyCPML2Hfield(){
	for (int i = 0; i<b.getCellCountRn(); i++){
	    for(int j = 0; j<ps.getNA()+1; j++){
		for(int k = 0; k<ps.getNZ(); k++){
		    Psi_har_rn[i][j][k] = cpml_b_mr_rn[i]*Psi_har_rn[i][j][k] + cpml_a_mr_rn[i]*(Field.getEZ(i+1,j,k) - Field.getEZ(i,j,k)); 
		    Field.setHA(i,j,k, Field.getHA(i,j,k) + CPsi_har_rn[i][j][k]*Psi_har_rn[i][j][k]);
		}
	    }
	}
	
	//is_cpml_rn == TRUE && j<na
	for (int i = 0; i<b.getCellCountRn(); i++){
	    for(int j = 0; j<ps.getNA(); j++){
		for(int k = 0; k<ps.getNZ()+1; k++){
		    Psi_hzr_rn[i][j][k] = cpml_b_mr_rn[i]*Psi_hzr_rn[i][j][k] + cpml_a_mr_rn[i]*(Field.getEA(i+1,j,k) - Field.getEA(i,j,k)); 
		    Field.setHZ(i,j,k, Field.getHZ(i,j,k) + CPsi_hzr_rn[i][j][k]*Psi_hzr_rn[i][j][k]);    
		}
	    }
	}
	
	//is_cpml_rp == TRUE && k<nz
	for (int i = 0; i<b.getCellCountRp(); i++){
	    for(int j = 0; j<ps.getNA()+1; j++){
		for(int k = 0; k<ps.getNZ(); k++){
		    Psi_har_rp[i][j][k] = cpml_b_mr_rp[i]*Psi_har_rp[i][j][k] + cpml_a_mr_rp[i]*(Field.getEZ(i+(ps.getNR() - b.getCellCountRp())+1,j,k) - Field.getEZ(i+(ps.getNR() - b.getCellCountRp()),j,k)); 
		    Field.setHA((ps.getNR() - b.getCellCountRp())+i,j,k, Field.getHA((ps.getNR() - b.getCellCountRp())+i,j,k) + CPsi_har_rp[i][j][k]*Psi_har_rp[i][j][k]);
		}
	    }
	}
	
	//is_cpml_rp == TRUE && j<na
	for (int i = 0; i<b.getCellCountRp(); i++){
	    for(int j = 0; j<ps.getNA(); j++){
		for(int k = 0; k<ps.getNZ()+1; k++){
		    Psi_hzr_rp[i][j][k] = cpml_b_mr_rp[i]*Psi_hzr_rp[i][j][k] + cpml_a_mr_rp[i]*(Field.getEA(i+(ps.getNR() - b.getCellCountRp())+1,j,k) - Field.getEA(i+(ps.getNR() - b.getCellCountRp()),j,k)); 
		    Field.setHZ((ps.getNR() - b.getCellCountRp())+i,j,k, Field.getHZ((ps.getNR() - b.getCellCountRp())+i,j,k) + CPsi_hzr_rp[i][j][k]*Psi_hzr_rp[i][j][k]);    
		}
	    }
	}
	
	//is_cpml_an == TRUE && i<nr
	for(int i=0; i<ps.getNR(); i++){
	    for (int j = 0; j<b.getCellCountAn();j++){
		for (int k = 0; k<ps.getNZ()+1; k++){
		    Psi_hza_an[i][j][k] = cpml_b_ma_an[j]*Psi_hza_an[i][j][k] + cpml_a_ma_an[j]*(Field.getER(i,j+1,k) - Field.getER(i,j,k)); 
		    Field.setHZ(i,j,k, Field.getHZ(i,j,k) + CPsi_hza_an[i][j][k]*Psi_hza_an[i][j][k]);
		}
	    }
	}
	
	//is_cpml_an == TRUE && k<nz
	for(int i=0; i<ps.getNR()+1; i++){
	    for (int j = 0; j<b.getCellCountAn();j++){
		for (int k = 0; k<ps.getNZ(); k++){
		    Psi_hra_an[i][j][k] = cpml_b_ma_an[j]*Psi_hra_an[i][j][k] + cpml_a_ma_an[j]*(Field.getEZ(i,j+1,k) - Field.getEZ(i,j,k)); 
		    Field.setHR(i,j,k, Field.getHR(i,j,k) + CPsi_hra_an[i][j][k]*Psi_hra_an[i][j][k]);
		}
	    }
	}
	
	//is_cpml_ap == TRUE && i<nr
	for (int i=0;i<ps.getNR();i++){
	    for (int j = 0; j<b.getCellCountAp();j++){
		for (int k = 0;k<ps.getNZ()+1;k++){
		    Psi_hza_ap[i][j][k] = cpml_b_ma_ap[j] * Psi_hza_ap[i][j][k] + cpml_a_ma_ap[j]*(Field.getER(i,j+(ps.getNA() - b.getCellCountAp())+1,k) - Field.getER(i,j+(ps.getNA() - b.getCellCountAp()),k)); 
		    Field.setHZ(i,(ps.getNA() - b.getCellCountAp())+j,k, Field.getHZ(i,(ps.getNA() - b.getCellCountAp())+j,k) + CPsi_hza_ap[i][j][k]*Psi_hza_ap[i][j][k]);
		}
	    }
	}
	
	//is_cpml_ap == TRUE && k<nz
	for (int i=0;i<ps.getNR()+1;i++){
	    for (int j = 0; j<b.getCellCountAp();j++){
		for (int k = 0;k<ps.getNZ();k++){
		    Psi_hra_ap[i][j][k] = cpml_b_ma_ap[j] * Psi_hra_ap[i][j][k] + cpml_a_ma_ap[j]*(Field.getEZ(i,j+(ps.getNA() - b.getCellCountAp())+1,k) - Field.getEZ(i,j+(ps.getNA() - b.getCellCountAp()),k)); 
		    Field.setHR(i,(ps.getNA() - b.getCellCountAp())+j,k, Field.getHR(i,(ps.getNA() - b.getCellCountAp())+j,k) + CPsi_hra_ap[i][j][k]*Psi_hra_ap[i][j][k]);
		}
	    }
	}
	
	//is_cpml_zn == TRUE && j<na
	for (int i = 0; i<ps.getNR()+1;i++){
	    for (int j = 0; j<ps.getNA();j++){
		for (int k = 0; k<b.getCellCountZn(); k++){
		    Psi_hrz_zn[i][j][k] = cpml_b_mz_zn[k] * Psi_hrz_zn[i][j][k] + cpml_a_mz_zn[k]*(Field.getEA(i,j,k+1) - Field.getEA(i,j,k));
		    Field.setHR(i,j,k, Field.getHR(i,j,k) + CPsi_hrz_zn[i][j][k]*Psi_hrz_zn[i][j][k]);
		}
	    }
	}
	
	//is_cpml_zn == TRUE && i<nr
	for (int i = 0; i<ps.getNR();i++){
	    for (int j = 0; j<ps.getNA()+1;j++){
		for (int k = 0; k<b.getCellCountZn(); k++){
		    Psi_haz_zn[i][j][k] = cpml_b_mz_zn[k] * Psi_haz_zn[i][j][k] + cpml_a_mz_zn[k]*(Field.getER(i,j,k+1) - Field.getER(i,j,k));
		    Field.setHA(i,j,k, Field.getHA(i,j,k) + CPsi_haz_zn[i][j][k]*Psi_haz_zn[i][j][k]);
		}	
	    }
	}
	
	//is_cpml_zp ==TRUE && j<na
	for (int i = 0; i<ps.getNR()+1; i++){
	    for (int j = 0; j<ps.getNA(); j++){
		for (int k = 0; k<b.getCellCountZp(); k++){
		    Psi_hrz_zp[i][j][k] = cpml_b_mz_zp[k] * Psi_hrz_zp[i][j][k] + cpml_a_mz_zp[k]*(Field.getEA(i,j,k+(ps.getNZ() - b.getCellCountZp())+1) - Field.getEA(i,j,k+(ps.getNZ() - b.getCellCountZp()))); 
		    Field.setHR(i,j,(ps.getNZ() - b.getCellCountZp())+k, Field.getHR(i,j,(ps.getNZ() - b.getCellCountZp())+k) + CPsi_hrz_zp[i][j][k]*Psi_hrz_zp[i][j][k]);
		}
	    }
	}
	
	//is_cpml_zp ==TRUE && i<nr
	for (int i = 0; i<ps.getNR(); i++){
	    for (int j = 0; j<ps.getNA()+1; j++){
		for (int k = 0; k<b.getCellCountZp(); k++){					
		    Psi_haz_zp[i][j][k] = cpml_b_mz_zp[k] * Psi_haz_zp[i][j][k] + cpml_a_mz_zp[k]*(Field.getER(i,j,k+(ps.getNZ() - b.getCellCountZp())+1) - Field.getER(i,j,k+(ps.getNZ() - b.getCellCountZp()))); 
		    Field.setHA(i,j,(ps.getNZ() - b.getCellCountZp())+k, Field.getHA(i,j,(ps.getNZ() - b.getCellCountZp())+k) + CPsi_haz_zp[i][j][k]*Psi_haz_zp[i][j][k]);
		}    
	    }
	}
	
    }
    
    public void applyCPML2Efield(){
        int temp;
	for (int i = 0; i<b.getCellCountRn(); i++){
	    for (int j = 0; j<ps.getNA(); j++){
		for (int k = 0; k<ps.getNZ()+1;k++){
		    Psi_ear_rn[i][j][k] = cpml_b_er_rn[i]*Psi_ear_rn[i][j][k] + cpml_a_er_rn[i]*(Field.getHZ(i+1,j,k) - Field.getHZ(i,j,k));
		    Field.setEA(i+1,j,k, Field.getEA(i+1,j,k) + CPsi_ear_rn[i][j][k]*Psi_ear_rn[i][j][k]);
		}
	    }
	}
	
	//is_cpml_rn == TRUE && k<nz
	for (int i = 0; i<b.getCellCountRn(); i++){
	    for (int j = 0; j<ps.getNA()+1; j++){
		for (int k = 0; k<ps.getNZ();k++){ 
		    Psi_ezr_rn[i][j][k] = cpml_b_er_rn[i]*Psi_ezr_rn[i][j][k] + cpml_a_er_rn[i]*(Field.getHA(i+1,j,k) - Field.getHA(i,j,k)); 
		    Field.setEZ(i+1,j,k, Field.getEZ(i+1,j,k) + CPsi_ezr_rn[i][j][k]*Psi_ezr_rn[i][j][k]);
		}
	    }
	}
	//is_cpml_rp == TRUE && j<na
	temp = ps.getNR() - b.getCellCountRp();
	for (int i = 0; i<b.getCellCountRp(); i++){
	    for (int j = 0; j<ps.getNA(); j++){
		for (int k =0; k<ps.getNZ()+1; k++){
		    Psi_ear_rp[i][j][k] = cpml_b_er_rp[i]*Psi_ear_rp[i][j][k] + cpml_a_er_rp[i]*(Field.getHZ(i+temp,j,k) - Field.getHZ(i+temp-1,j,k));
		    Field.setEA(temp+i,j,k, Field.getEA(temp+i,j,k) + CPsi_ear_rp[i][j][k]*Psi_ear_rp[i][j][k]);
		}
	    }
	}
	
	//is_cpml_rp == TRUE && k<nz
	for (int i = 0; i<b.getCellCountRp(); i++){
	    for (int j = 0; j<ps.getNA()+1; j++){
		for (int k =0; k<ps.getNZ(); k++){
		    Psi_ezr_rp[i][j][k] = cpml_b_er_rp[i]*Psi_ezr_rp[i][j][k] + cpml_a_er_rp[i]*(Field.getHA(i+temp,j,k) - Field.getHA(i+temp-1,j,k));
		    Field.setEZ(temp+i,j,k, Field.getEZ(temp+i,j,k) + CPsi_ezr_rp[i][j][k]*Psi_ezr_rp[i][j][k]);
		}
	    }
	}
	
	//is_cpml_an == TRUE && k<nz
	for (int i = 0; i<ps.getNR()+1;i++){
	    for (int j = 0; j<b.getCellCountAn(); j++){
		for (int k = 0; k<ps.getNZ(); k++){
		    Psi_eza_an[i][j][k] = cpml_b_ea_an[j] * Psi_eza_an[i][j][k] + cpml_a_ea_an[j]*(Field.getHR(i,j+1,k) - Field.getHR(i,j,k)); 
		    Field.setEZ(i,j+1,k, Field.getEZ(i,j+1,k) + CPsi_eza_an[i][j][k]* Psi_eza_an[i][j][k]);
		}
	    }
	}
	
	//is_cpml_an == TRUE && i<nr
	for (int i = 0; i<ps.getNR();i++){
	    for (int j = 0; j<b.getCellCountAn(); j++){
		for (int k = 0; k<ps.getNZ()+1; k++){
		    Psi_era_an[i][j][k] = cpml_b_ea_an[j] * Psi_era_an[i][j][k] + cpml_a_ea_an[j]*(Field.getHZ(i,j+1,k) - Field.getHZ(i,j,k)); 
		    Field.setER(i,j+1,k, Field.getER(i,j+1,k) + CPsi_era_an[i][j][k]* Psi_era_an[i][j][k]);
		}    
	    }
	}
	
	//is_cpml_ap == TRUE && k<nz
	temp = ps.getNA() - b.getCellCountAp();
	for(int i = 0; i<ps.getNR()+1; i++){
	    for (int j = 0; j<b.getCellCountAp(); j++){
		for (int k = 0; k<ps.getNZ(); k++){
		    Psi_eza_ap[i][j][k] = cpml_b_ea_ap[j]*Psi_eza_ap[i][j][k] + cpml_a_ea_ap[j]*(Field.getHR(i,j+temp,k) - Field.getHR(i,j+temp-1,k)); 
		    Field.setEZ(i,temp+j,k, Field.getEZ(i,temp+j,k) + CPsi_eza_ap[i][j][k]*Psi_eza_ap[i][j][k]);
		}
	    }
	}
	
	//is_cpml_ap == TRUE && i<nr
	for(int i = 0; i<ps.getNR(); i++){
	    for (int j = 0; j<b.getCellCountAp(); j++){
		for (int k = 0; k<ps.getNZ()+1; k++){
		    Psi_era_ap[i][j][k] = cpml_b_ea_ap[j]*Psi_era_ap[i][j][k] + cpml_a_ea_ap[j]*(Field.getHZ(i,j+temp,k) - Field.getHZ(i,j+temp-1,k));
		    Field.setER(i,temp+j,k, Field.getER(i,temp+j,k) + CPsi_era_ap[i][j][k]*Psi_era_ap[i][j][k]);
		}
	    }
	}
	
	//is_cpml_zn == TRUE && i<nr
	for (int i =0; i<ps.getNR();i++){
	    for (int j=0;j<ps.getNA()+1;j++){
		for (int k = 0; k<b.getCellCountZn(); k++){
		    Psi_erz_zn[i][j][k] = cpml_b_ez_zn[k] * Psi_erz_zn[i][j][k] + cpml_a_ez_zn[k]*(Field.getHA(i,j,k+1) - Field.getHA(i,j,k)); 
		    Field.setER(i,j,k+1, Field.getER(i,j,k+1) + CPsi_erz_zn[i][j][k]*Psi_erz_zn[i][j][k]);
		}
	    }
	}
	
	//is_cpml_zn == TRUE && j<na
	for (int i =0; i<ps.getNR()+1;i++){
	    for (int j=0;j<ps.getNA();j++){
		for (int k = 0; k<b.getCellCountZn(); k++){
		    Psi_eaz_zn[i][j][k] = cpml_b_ez_zn[k] * Psi_eaz_zn[i][j][k] + cpml_a_ez_zn[k]*(Field.getHR(i,j,k+1) - Field.getHR(i,j,k)); 
		    Field.setEA(i,j,k+1, Field.getEA(i,j,k+1) + CPsi_eaz_zn[i][j][k]*Psi_eaz_zn[i][j][k]);
		}
	    }
	}
	
	//is_cpml_zp = TRUE && i<nr
	temp = ps.getNZ() - b.getCellCountZp();
	for (int i = 0; i<ps.getNR();i++){
	    for (int j = 0; j<ps.getNA()+1; j++){
		for (int k = 0; k<b.getCellCountZp(); k++){
		    Psi_erz_zp[i][j][k] = cpml_b_ez_zp[k]*Psi_erz_zp[i][j][k] + cpml_a_ez_zp[k]*(Field.getHA(i,j,k+temp)-Field.getHA(i,j,k+temp-1)); 	
		    Field.setER(i,j,temp+k, Field.getER(i,j,temp+k) + CPsi_erz_zp[i][j][k]*Psi_erz_zp[i][j][k]);
		}
	    }
	}
	
	//is_cpml_zp = TRUE && j<na
	for (int i = 0; i<ps.getNR()+1;i++){
	    for (int j = 0; j<ps.getNA(); j++){
		for (int k = 0; k<b.getCellCountZp(); k++){
		    Psi_eaz_zp[i][j][k] = cpml_b_ez_zp[k]*Psi_eaz_zp[i][j][k] + cpml_a_ez_zp[k]*(Field.getHR(i,j,k+temp)-Field.getHR(i,j,k+temp-1));
		    Field.setEA(i,j,temp+k, Field.getEA(i,j,temp+k) + CPsi_eaz_zp[i][j][k]*Psi_eaz_zp[i][j][k]);
		}
	    }
	}
    }
}
