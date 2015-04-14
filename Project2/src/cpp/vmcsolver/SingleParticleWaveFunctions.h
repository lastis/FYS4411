namespace wave_functions{
    static double alpha;

    static double phi1s(double r){
        return exp(-alpha*r);
    }

    static double phi2s(double r){
        return (1-alpha*r/2)*exp(-alpha*r/2);
    }

    static double phi2p(double r){
        return alpha*r*exp(-alpha*r/2);
    }
}
