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
public class ProblemSpace {
    private double Rmin=10000;
    private double Amin=10000;
    private double Zmin=10000;
    private double Rmax=-10000;
    private double Amax=-10000;
    private double Zmax=-10000;
    private double sizeR;
    private double sizeA;
    private double sizeZ;
    private final int nr,na,nz;
    private int[][][] ms;
    private final Brick[] bricks;
    private int brickCount;
    private final Cell c;
    private final Boundary b;
    private MaterialGrid mg;
    private Material[] m;
    private int blr,bla,blz,bur,bua,buz;
    
    public ProblemSpace(int brickCount, Brick[] bricks,Cell c, Boundary b, Material[] m){
        this.brickCount = brickCount;
        this.bricks = bricks;
        this.c = c;
        this.b = b;
        this.m = m;
        
        for(int i = 0; i < brickCount; i++)
	{
            if(this.Rmin > bricks[i].getRmin())
                this.Rmin = bricks[i].getRmin();
            if(this.Amin > bricks[i].getAmin())
                this.Amin = bricks[i].getAmin();
            if(this.Zmin > bricks[i].getZmin())
                this.Zmin = bricks[i].getZmin();
            if(this.Rmax < bricks[i].getRmax())
                this.Rmax = bricks[i].getRmax(); 
            if(this.Amax < bricks[i].getAmax())
                this.Amax = bricks[i].getAmax(); 
            if(this.Zmax < bricks[i].getZmax())
                this.Zmax = bricks[i].getZmax();
	}
        
        this.Rmin = this.Rmin - c.getDeltaR() * b.getAirBufferCellCountRn();
        this.Amin = this.Amin - c.getDeltaA() * b.getAirBufferCellCountAn();
        this.Zmin = this.Zmin - c.getDeltaZ() * b.getAirBufferCellCountZn();
        this.Rmax = this.Rmax + c.getDeltaR() * b.getAirBufferCellCountRp();
        this.Amax = this.Amax + c.getDeltaA() * b.getAirBufferCellCountAp();
        this.Zmax = this.Zmax + c.getDeltaZ() * b.getAirBufferCellCountZp();
        
	sizeR = this.Rmax - this.Rmin;
	sizeA = this.Amax - this.Amin;
	sizeZ = this.Zmax - this.Zmin;
        
	this.nr = (int)Math.round(sizeR/c.getDeltaR());
	this.na = (int)Math.round(sizeA/c.getDeltaA());
	this.nz = (int)Math.round(sizeZ/c.getDeltaZ());
        
        sizeR = nr*c.getDeltaR();
        sizeA = na*c.getDeltaA();
        sizeZ = nz*c.getDeltaZ();
        
	this.Rmax = this.Rmin + sizeR;
	this.Amax = this.Amin + sizeA;
	this.Zmax = this.Zmin + sizeZ;
    }

    public MaterialGrid getMaterialGrid(){
        return this.mg;
    }
    
    public double getRmax(){
        return this.Rmax;
    }
    
    public double getAmax(){
        return this.Amax;
    }

    public double getZmax(){
        return this.Zmax;
    }
        
    public double getRmin(){
        return this.Rmin;
    }

    public double getAmin(){
        return this.Amin;
    }

    public double getZmin(){
        return this.Zmin;
    }
    
    public int getNR(){
        return nr;
    }
    
    public int getNA(){
        return na;
    }
    
    public int getNZ(){
        return nz;
    }
    
    public void initMaterialGrid(){
        this.ms = Common.genInt3DArray(nr, na, nz, 1);
        double sigmaPEC;
        System.out.println("assign material type of the brick to the cells...");
	for(int index = 0; index < this.brickCount; index++){
            //convert brick end coordinates to node indices 
            blr = (int)Math.round((this.bricks[index].getRmin() - this.Rmin)/c.getDeltaR()); 
            bla = (int)Math.round((this.bricks[index].getAmin() - this.Amin)/c.getDeltaA()); 
            blz = (int)Math.round((this.bricks[index].getZmin() - this.Zmin)/c.getDeltaZ()); 

            bur = (int)Math.round((this.bricks[index].getRmax() - this.Rmin)/c.getDeltaR()); 
            bua = (int)Math.round((this.bricks[index].getAmax() - this.Amin)/c.getDeltaA()); 
            buz = (int)Math.round((this.bricks[index].getZmax() - this.Zmin)/c.getDeltaZ());
            
            for(int i=blr;i <= bur-1;i++){
		for(int j=bla;j <= bua-1;j++){
                    for(int k=blz;k <= buz-1;k++){
			this.ms[i][j][k] = this.bricks[index].getMaterialType();
                    }
		}
            }
	}
        
        mg = new MaterialGrid(this.nr,this.na,this.nz,c,m,this.ms);
        mg.setAll();
        mg.averageAll();
        
        System.out.println("find the zero thickness bricks...");        
        for(int index = 0; index < this.brickCount; index++){
            sigmaPEC = this.m[this.bricks[index].getMaterialType()].sigmaE();

            //convert brick end coordinates to node indices 
            blr = (int)Math.round((this.bricks[index].getRmin() - this.Rmin)/c.getDeltaR()); 
            bla = (int)Math.round((this.bricks[index].getAmin() - this.Amin)/c.getDeltaA()); 
            blz = (int)Math.round((this.bricks[index].getZmin() - this.Zmin)/c.getDeltaZ()); 

            bur = (int)Math.round((this.bricks[index].getRmax() - this.Rmin)/c.getDeltaR()); 
            bua = (int)Math.round((this.bricks[index].getAmax() - this.Amin)/c.getDeltaA()); 
            buz = (int)Math.round((this.bricks[index].getZmax() - this.Zmin)/c.getDeltaZ());
            
            if (blr == bur){
		for(int indexBa = bla; indexBa <= bua; indexBa++){
                    for(int indexBz = blz; indexBz <= buz; indexBz++){
			if(indexBa < bua){  
                            mg.setSigmaEA(blr, indexBa, indexBz, sigmaPEC);
			}
			if(indexBz < buz){
                            mg.setSigmaEZ(blr, indexBa, indexBz, sigmaPEC);
			}
                    }
		}
            }
		
            if (bla == bua){
		for(int indexBr = blr;indexBr <= bur;indexBr++){
                    for(int indexBz = blz;indexBz <= buz;indexBz++){
			if(indexBz < buz){
                            mg.setSigmaEZ(indexBr, bla, indexBz, sigmaPEC);
			}
			if(indexBr < bur){
                            mg.setSigmaER(indexBr, bla, indexBz, sigmaPEC);
			}
                    }
		}
            }
		
            if (blz == buz){
		for(int indexBr = blr;indexBr <= bur;indexBr++){
                    for(int indexBa = bla;indexBa <= bua;indexBa++){
			if(indexBr < bur){
                            mg.setSigmaER(indexBr, indexBa, blz, sigmaPEC);
			}
			if(indexBa < bua){	
                            mg.setSigmaEA(indexBr, indexBa, blz, sigmaPEC);
			}
                    }
		}
            }
        }
    }
}
