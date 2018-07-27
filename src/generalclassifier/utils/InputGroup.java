package generalclassifier.utils;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import com.sun.org.apache.bcel.internal.generic.BREAKPOINT;

public class InputGroup {

    public enum DistributionType {
        NORMAL,
        LOGNORMAL
    }

    Input<RealParameter> mean;
    Input<RealParameter> standardDev;
    Input<RealParameter> zeroFraction;
    DistributionType distributionType;

    public InputGroup(Input<RealParameter> m, Input<RealParameter> sd, DistributionType type) {
        this.mean = m;
        this.standardDev = sd;
        this.distributionType = type;
    }

    public InputGroup(Input<RealParameter> m, Input<RealParameter> sd, Input<RealParameter> zeroFraction, DistributionType type) {
        this.mean = m;
        this.standardDev = sd;
        this.zeroFraction = zeroFraction;
        this.distributionType = type;
    }

    public Input<RealParameter> getMean(){
        return this.mean;
    }

    public Input<RealParameter> getStandardDev(){
        return this.standardDev;
    }

    public Input<RealParameter> getZeroFraction() {
        return zeroFraction;
    }

    public DistributionType getDistributionType() {
        return distributionType;
    }
}


