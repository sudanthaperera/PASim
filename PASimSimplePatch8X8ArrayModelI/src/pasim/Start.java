package pasim;

public class Start {
    public static void main(String[] args) {
        boolean optimization = false;
        boolean finiteArraySimulation = true;
        boolean isolatedElement = false;
        
        if(isolatedElement){
            IsolatedElement IE = new IsolatedElement();
            IE.start();
        }
        
        if(finiteArraySimulation){
            DummySim ds = new DummySim();
            ds.run();
        }
        else if(optimization){
            ElementOptimization EO = new ElementOptimization();
        }
    }
}
