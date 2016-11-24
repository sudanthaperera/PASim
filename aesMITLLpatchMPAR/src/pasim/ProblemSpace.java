package pasim;

public class ProblemSpace {
    private double Xmin=10000;
    private double Ymin=10000;
    private double Zmin=10000;
    private double Xmax=-10000;
    private double Ymax=-10000;
    private double Zmax=-10000;
    private double sizeX;
    private double sizeY;
    private double sizeZ;
    private final int nx,ny,nz;
    private int[][][] ms;
    private final Brick[] bricks;
    private int brickCount;
    private final Cell c;
    private final Boundary b;
    private MaterialGrid mg;
    private Material[] m;
    private int blx,bly,blz,bux,buy,buz;
    
    public ProblemSpace(int brickCount, Brick[] bricks,Cell c, Boundary b, Material[] m){
        this.brickCount = brickCount;
        this.bricks = bricks;
        this.c = c;
        this.b = b;
        this.m = m;
        
        for(int i = 0; i < brickCount; i++)
	{
            if(this.Xmin > bricks[i].getXmin())
                this.Xmin = bricks[i].getXmin();
            if(this.Ymin > bricks[i].getYmin())
                this.Ymin = bricks[i].getYmin();
            if(this.Zmin > bricks[i].getZmin())
                this.Zmin = bricks[i].getZmin();
            if(this.Xmax < bricks[i].getXmax())
                this.Xmax = bricks[i].getXmax(); 
            if(this.Ymax < bricks[i].getYmax())
                this.Ymax = bricks[i].getYmax(); 
            if(this.Zmax < bricks[i].getZmax())
                this.Zmax = bricks[i].getZmax();
	}
        
        this.Xmin = this.Xmin - c.getDeltaX() * b.getAirBufferCellCountXn();
        this.Ymin = this.Ymin - c.getDeltaY() * b.getAirBufferCellCountYn();
        this.Zmin = this.Zmin - c.getDeltaZ() * b.getAirBufferCellCountZn();
        this.Xmax = this.Xmax + c.getDeltaX() * b.getAirBufferCellCountXp();
        this.Ymax = this.Ymax + c.getDeltaY() * b.getAirBufferCellCountYp();
        this.Zmax = this.Zmax + c.getDeltaZ() * b.getAirBufferCellCountZp();
        
	sizeX = this.Xmax - this.Xmin;
	sizeY = this.Ymax - this.Ymin;
	sizeZ = this.Zmax - this.Zmin;
        
	this.nx = (int)Math.round(sizeX/c.getDeltaX());
	this.ny = (int)Math.round(sizeY/c.getDeltaY());
	this.nz = (int)Math.round(sizeZ/c.getDeltaZ());
        
        sizeX = nx*c.getDeltaX();
        sizeY = ny*c.getDeltaY();
        sizeZ = nz*c.getDeltaZ();
        
	this.Xmax = this.Xmin + sizeX;
	this.Ymax = this.Ymin + sizeY;
	this.Zmax = this.Zmin + sizeZ;
    }

    public MaterialGrid getMaterialGrid(){
        return this.mg;
    }
    
    public double getXmax(){
        return this.Xmax;
    }
    
    public double getYmax(){
        return this.Ymax;
    }

    public double getZmax(){
        return this.Zmax;
    }
        
    public double getXmin(){
        return this.Xmin;
    }

    public double getYmin(){
        return this.Ymin;
    }

    public double getZmin(){
        return this.Zmin;
    }
    
    public int getNX(){
        return nx;
    }
    
    public int getNY(){
        return ny;
    }
    
    public int getNZ(){
        return nz;
    }
    
    public void initMaterialGrid(){
        this.ms = Common.genInt3DArray(nx, ny, nz, 1);
        double sigmaPEC;
	for(int index = 0; index < this.brickCount; index++){
            //convert brick end coordinates to node indices 
            blx = (int)Math.round((this.bricks[index].getXmin() - this.Xmin)/c.getDeltaX());
            bly = (int)Math.round((this.bricks[index].getYmin() - this.Ymin)/c.getDeltaY());
            blz = (int)Math.round((this.bricks[index].getZmin() - this.Zmin)/c.getDeltaZ());

            bux = (int)Math.round((this.bricks[index].getXmax() - this.Xmin)/c.getDeltaX());
            buy = (int)Math.round((this.bricks[index].getYmax() - this.Ymin)/c.getDeltaY());
            buz = (int)Math.round((this.bricks[index].getZmax() - this.Zmin)/c.getDeltaZ());
            
            for(int i=blx;i <= bux-1;i++){
		for(int j=bly;j <= buy-1;j++){
                    for(int k=blz;k <= buz-1;k++){
			this.ms[i][j][k] = this.bricks[index].getMaterialType();
                    }
		}
            }
	}
        
        mg = new MaterialGrid(this.nx,this.ny,this.nz,c,m,this.ms);
        mg.setAll();
        mg.averageAll();
        
        for(int index = 0; index < this.brickCount; index++){
            sigmaPEC = this.m[this.bricks[index].getMaterialType()].sigmaE();

            //convert brick end coordinates to node indices 
            blx = (int)Math.round((this.bricks[index].getXmin() - this.Xmin)/c.getDeltaX());
            bly = (int)Math.round((this.bricks[index].getYmin() - this.Ymin)/c.getDeltaY());
            blz = (int)Math.round((this.bricks[index].getZmin() - this.Zmin)/c.getDeltaZ());

            bux = (int)Math.round((this.bricks[index].getXmax() - this.Xmin)/c.getDeltaX());
            buy = (int)Math.round((this.bricks[index].getYmax() - this.Ymin)/c.getDeltaY());
            buz = (int)Math.round((this.bricks[index].getZmax() - this.Zmin)/c.getDeltaZ());
            
            if (blx == bux){
		for(int indexBy = bly; indexBy <= buy; indexBy++){
                    for(int indexBz = blz; indexBz <= buz; indexBz++){
			if(indexBy < buy){  
                            mg.setSigmaEY(blx, indexBy, indexBz, sigmaPEC);
			}
			if(indexBz < buz){
                            mg.setSigmaEZ(blx, indexBy, indexBz, sigmaPEC);
			}
                    }
		}
            }
		
            if (bly == buy){
		for(int indexBx = blx;indexBx <= bux;indexBx++){
                    for(int indexBz = blz;indexBz <= buz;indexBz++){
			if(indexBz < buz){
                            mg.setSigmaEZ(indexBx, bly, indexBz, sigmaPEC);
			}
			if(indexBx < bux){
                            mg.setSigmaEX(indexBx, bly, indexBz, sigmaPEC);
			}
                    }
		}
            }
		
            if (blz == buz){
		for(int indexBx = blx;indexBx <= bux;indexBx++){
                    for(int indexBy = bly;indexBy <= buy;indexBy++){
			if(indexBx < bux){
                            mg.setSigmaEX(indexBx, indexBy, blz, sigmaPEC);
			}
			if(indexBy < buy){	
                            mg.setSigmaEY(indexBx, indexBy, blz, sigmaPEC);
			}
                    }
		}
            }
        }
    }
}
