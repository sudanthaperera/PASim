package pasim;

/**
 *
 * @author Sudantha
 */
public class Chromosome {
    public static int xDim;
    public static int yDim;
    public static int layerCount;
    public static int portCount;
    public static int portDim;
    
    public Chromosome(){
        
    }
    
    public Chromosome(int xDim,int yDim,int layerCount,int portCount,int portDim){
        Chromosome.xDim = xDim;
        Chromosome.yDim = yDim;
        Chromosome.portCount = portCount;
        Chromosome.layerCount = layerCount;
        Chromosome.portDim = portDim;
    }
    
    public static int xDim(){
        return Chromosome.xDim;
    }
    
    public static int yDim(){
        return Chromosome.xDim;
    }
    
    public static int layerCount(){
        return Chromosome.layerCount;
    }
    
    public static int portCount(){
        return Chromosome.portCount;
    }
    
    public static int portDim(){
        return Chromosome.portDim;
    }
    
    public static void xDim(int xDim){
        Chromosome.xDim = xDim;
    }
    
    public static void yDim(int yDim){
        Chromosome.yDim = yDim;
    }
    
    public static void layerCount(int layerCount){
        Chromosome.layerCount = layerCount;
    }
    
    public static void portCount(int portCount){
        Chromosome.portCount = portCount;
    }
    
    public static void portDim(int portDim){
        Chromosome.portDim = portDim;
    }
}
