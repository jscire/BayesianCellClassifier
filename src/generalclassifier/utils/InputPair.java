package generalclassifier.utils;

import beast.core.Input;
import beast.core.parameter.RealParameter;

public class InputPair {
    Input<RealParameter> mean;
    Input<RealParameter> standardDev;

    public InputPair(Input<RealParameter> m, Input<RealParameter> sd) {
        this.mean = m;
        this.standardDev = sd;
    }

    public Input<RealParameter> getMean(){
        return this.mean;
    }

    public Input<RealParameter> getStandardDev(){
        return this.standardDev;
    }
}
