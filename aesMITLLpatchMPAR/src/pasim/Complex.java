package pasim;

public class Complex {
    private double real;
    private double imag;

    public Complex() {
        this.real = 0;
        this.imag = 0;
    }
    
    public Complex(double angle){
        this.real = Math.cos(angle);
        this.imag = Math.sin(angle);
    }

    public Complex( double real, double imag ) {
        this.real      = real;
        this.imag = imag;
    }
   
    public Complex( Complex val ) {
        this.real      = val.real;
        this.imag = val.imag;
    }
    
    public void set(Complex val){
        this.real = val.real();
        this.imag = val.imag();
    }
    
    public double real(){
        return real;
    }

    public double imag(){
        return imag;
    }
    
    public void setReal(double real){
        this.real = real;
    }

    public void setImag(double imag){
        this.imag = imag;
    } 
    
    public String toString() {
        if (imag >= 0)
            return real + "+" +  imag + "i";
        else
            return real + "-" + -imag + "i";
    }

    public Complex conj() {
        Complex val = new Complex();

        val.real      = this.real;
        val.imag = - this.imag;

        return val;
    }

    public Complex plus( Complex num ) {
        Complex val = new Complex();

        val.real      = real + num.real;
        val.imag = imag + num.imag;

        return val;
    }

    public Complex minus( Complex num ) {
        Complex val = new Complex();

        val.real      = real - num.real;
        val.imag = imag - num.imag;

        return val;
    }

    public Complex times( Complex num ) {
        Complex val = new Complex();

        val.real      = real*num.real - imag*num.imag;
        val.imag = imag*num.real + real*num.imag;

        return val;
    }
   
    public Complex times(double num) {
        Complex val = new Complex();

        val.real      = real*num;
        val.imag = imag*num;

        return val;
    }

    public Complex divide( Complex num ) {
        double denominator = Math.pow(num.real(),2) + Math.pow(num.imag(), 2);
        double numeratorReal = this.real*num.real()+this.imag*num.imag();
        double numeratorImag = this.imag*num.real()-this.real*num.imag();
        
        if(denominator ==0){
            return (new Complex(0.0,0.0));
        }
        else{
            return (new Complex(numeratorReal/denominator,numeratorImag/denominator));
        }
    }
   
    public Complex divide( double num ) {
        Complex val = new Complex(this.real/num, this.imag/num);
        return val;
    }

    public double abs() {
        return Math.sqrt((Math.pow(real, 2) + Math.pow(imag, 2)));
    }
    
    public double angle(){
        return Math.atan(imag/real);
    }

    public Complex sqrt() { 
        double r, sinTheta, cosTheta;
        double sqrtR, sinHalfOfTheta, cosHalfOfTheta;
        
        r = Math.sqrt(Math.pow(real, 2) + Math.pow(imag, 2));
        cosTheta = real/r;
        sinTheta = imag/r;
        
        sqrtR = Math.sqrt(r);
        sinHalfOfTheta = Math.sin(Math.asin(sinTheta)/2);
        cosHalfOfTheta = Math.cos(Math.acos(cosTheta)/2);

        return (new Complex(sqrtR*cosHalfOfTheta,sqrtR*sinHalfOfTheta));
    }
   
    public boolean isEqual(Complex that){
        return this.real == that.real() && this.imag == that.imag();
    }
}
