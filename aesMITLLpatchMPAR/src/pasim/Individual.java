package pasim;

public class Individual {
    private int[] chromosome;
    private double fitness = -1;

    public Individual(int[] chromosome) {
        this.chromosome = chromosome;
    }

    public Individual(int chromosomeLength) {
        Chromosome c = new Chromosome();
	this.chromosome = new int[chromosomeLength];
        
        for (int gene = 0; gene < chromosomeLength; gene++) {
            if (0.5 < Math.random()) {
		chromosome[gene] = 1;
            } else {
		chromosome[gene] = 0;
            }
	}
        
        this.makeSymmetric();
    }
    
    public void makeSymmetric(){
        Chromosome c = new Chromosome();
        int[][] patch = new int[c.xDim][c.yDim];
        int[] port = new int[c.portDim];
        int[] portLocationH = {c.xDim/2, c.yDim/2};
        int[] portLocationV = {c.xDim/2, c.yDim/2};
        
        for (int patchGeneX = 0; patchGeneX < c.xDim; patchGeneX++) {
            for (int patchGeneY = 0; patchGeneY < c.yDim; patchGeneY++) {
                patch[patchGeneX][patchGeneY] = this.getGene(patchGeneX*c.yDim+patchGeneY);
            }
        }
        
        for (int portGene = 0; portGene < c.portDim; portGene++) {
            port[portGene] = this.getGene(portGene + c.xDim*c.yDim);
            portLocationH[0] = portLocationH[0] - port[portGene]*(int)Math.pow(2, portGene);
            portLocationV[1] = portLocationV[1] + port[portGene]*(int)Math.pow(2, portGene);
        }
        
        
        for (int patchGeneX = 0; patchGeneX < c.xDim - 1; patchGeneX++) {
            for (int patchGeneY = 0; patchGeneY < c.yDim - patchGeneX - 1; patchGeneY++) {
                patch[c.xDim - (patchGeneY + 1)][c.yDim - (patchGeneX + 1)] = patch[patchGeneX][patchGeneY];
            }
        }
        
        portLocationFilling(patch, portLocationH, portLocationV);
        
        for (int patchGeneX = 0; patchGeneX < c.xDim; patchGeneX++) {
            for (int patchGeneY = 0; patchGeneY < c.yDim; patchGeneY++) {
                this.setGene(patchGeneX*c.yDim+patchGeneY, patch[patchGeneX][patchGeneY]);
            }
        }
        
	for (int gene = c.xDim*c.yDim; gene < chromosome.length; gene++) {
            this.setGene(gene,port[gene - c.xDim*c.yDim]);
	}
    }
    
    private void portLocationFilling(int[][] patch, int[]  portLocationH, int[] portLocationV){
        patch[portLocationH[0]-1][portLocationH[1]-1] = 1;
        patch[portLocationH[0]-1][portLocationH[1]-0] = 1;
        patch[portLocationH[0]-0][portLocationH[1]-1] = 1;
        patch[portLocationH[0]+1][portLocationH[1]+1] = 1;
        patch[portLocationH[0]+1][portLocationH[1]+0] = 1;
        patch[portLocationH[0]+0][portLocationH[1]+1] = 1;
        patch[portLocationH[0]+1][portLocationH[1]-1] = 1;
        patch[portLocationH[0]-1][portLocationH[1]+1] = 1;
        patch[portLocationH[0]+0][portLocationH[1]+0] = 1;
        
        patch[portLocationV[0]-1][portLocationV[1]-1] = 1;
        patch[portLocationV[0]-1][portLocationV[1]-0] = 1;
        patch[portLocationV[0]-0][portLocationV[1]-1] = 1;
        patch[portLocationV[0]+1][portLocationV[1]+1] = 1;
        patch[portLocationV[0]+1][portLocationV[1]+0] = 1;
        patch[portLocationV[0]+0][portLocationV[1]+1] = 1;
        patch[portLocationV[0]+1][portLocationV[1]-1] = 1;
        patch[portLocationV[0]-1][portLocationV[1]+1] = 1;
        patch[portLocationV[0]+0][portLocationV[1]+0] = 1;
    }

    public int[] getChromosome() {
	return this.chromosome;
    }

    public int getChromosomeLength() {
        return this.chromosome.length;
    }

    public void setGene(int offset, int gene) {
	this.chromosome[offset] = gene;
    }

    public int getGene(int offset) {
	return this.chromosome[offset];
    }

    public void setFitness(double fitness) {
	this.fitness = fitness;
    }

    public double getFitness() {
	return this.fitness;
    }

    public String toString() {
	String output = "";
	for (int gene = 0; gene < this.chromosome.length; gene++) {
            output += this.chromosome[gene];
	}
        return output;
    }
    
    public void save(int genIndex){
        double[] a = new double[this.chromosome.length];
        
        for(int index=0; index<a.length;index++){
            a[index] = this.chromosome[index];
        }
        Common.save1DArray(a, "BestSolustion".concat(String.valueOf(genIndex)));
    }
}
