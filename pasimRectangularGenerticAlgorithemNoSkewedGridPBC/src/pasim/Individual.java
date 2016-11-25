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
        int[][] patch = new int[c.xDim][c.yDim];
        int[] port = new int[c.portDim];
        
        for (int portGene = 0; portGene < c.portDim; portGene++) {
            if (0.5 < Math.random()) {
		port[portGene] = 1;
            } else {
		port[portGene] = 0;
            }
	}
        
        for (int patchGeneX = 0; patchGeneX < c.xDim; patchGeneX++) {
            for (int patchGeneY = 0; patchGeneY < c.yDim; patchGeneY++) {
                if (0.5 < Math.random()) {
                    patch[patchGeneX][patchGeneY] = 1;
                } else {
                    patch[patchGeneX][patchGeneY] = 0;
                }
            }
	}
        
        for (int patchGeneX = 0; patchGeneX < c.xDim; patchGeneX++) {
            for (int patchGeneY = 0; patchGeneY < c.yDim; patchGeneY++) {
                patch[c.xDim - patchGeneY-1][c.yDim - patchGeneX-1] = patch[patchGeneX][patchGeneY];
            }
        }
        
        for (int patchGeneX = 0; patchGeneX < c.xDim; patchGeneX++) {
            for (int patchGeneY = 0; patchGeneY < c.yDim; patchGeneY++) {
                this.setGene(patchGeneX*c.yDim+patchGeneY, patch[patchGeneX][patchGeneY]);
            }
        }
        
	for (int gene = c.xDim*c.yDim; gene < chromosomeLength; gene++) {
            this.setGene(gene,port[gene - c.xDim*c.yDim]);
	}
    }
    
    public void makeSymmetric(){
        Chromosome c = new Chromosome();
        int[][] patch = new int[c.xDim][c.yDim];
        int[] port = new int[c.portDim];
        
        for (int patchGeneX = 0; patchGeneX < c.xDim; patchGeneX++) {
            for (int patchGeneY = 0; patchGeneY < c.yDim; patchGeneY++) {
                patch[patchGeneX][patchGeneY] = this.getGene(patchGeneX*c.yDim+patchGeneY);
            }
        }
        
        for (int portGene = 0; portGene < c.portDim; portGene++) {
            port[portGene] = this.getGene(portGene + c.xDim*c.yDim);
        }
        
        for (int patchGeneX = 0; patchGeneX < c.xDim; patchGeneX++) {
            for (int patchGeneY = 0; patchGeneY < c.yDim; patchGeneY++) {
                patch[c.xDim - patchGeneY][c.yDim - patchGeneX] = patch[patchGeneX][patchGeneY];
            }
        }
        
        for (int patchGeneX = 0; patchGeneX < c.xDim; patchGeneX++) {
            for (int patchGeneY = 0; patchGeneY < c.yDim; patchGeneY++) {
                this.setGene(patchGeneX*c.yDim+patchGeneY, patch[patchGeneX][patchGeneY]);
            }
        }
        
	for (int gene = c.xDim*c.yDim; gene < chromosome.length; gene++) {
            this.setGene(gene,port[gene - c.xDim*c.yDim]);
	}
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
}
