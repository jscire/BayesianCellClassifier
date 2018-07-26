package generalclassifier.utils;

import beast.core.Input;
import beast.core.parameter.RealParameter;

public class InputPair {

    public enum DistributionType {
        NORMAL,
        LOGNORMAL
    }

    Input<RealParameter> mean;
    Input<RealParameter> standardDev;
    DistributionType distributionType;

    public InputPair(Input<RealParameter> m, Input<RealParameter> sd, DistributionType type) {
        this.mean = m;
        this.standardDev = sd;
        this.distributionType = type;
    }

    public Input<RealParameter> getMean(){
        return this.mean;
    }

    public Input<RealParameter> getStandardDev(){
        return this.standardDev;
    }

    public DistributionType getDistributionType() {
        return distributionType;
    }
}


