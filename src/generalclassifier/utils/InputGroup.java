package generalclassifier.utils;

import beast.core.Input;
import beast.core.parameter.RealParameter;

public class InputGroup {

    public enum DistributionType {
        NORMAL,
        LOGNORMAL,
        GAMMA
    }

    Input<RealParameter> parm1;
    Input<RealParameter> parm2;
    Input<RealParameter> zeroFraction;
    DistributionType distributionType;
    boolean hasZeroFraction = false;

    public InputGroup(Input<RealParameter> parm1, Input<RealParameter> parm2, DistributionType type) {
        this.parm1 = parm1;
        this.parm2 = parm2;
        this.distributionType = type;
    }

    public InputGroup(Input<RealParameter> parm1, Input<RealParameter> parm2, Input<RealParameter> zeroFraction, DistributionType type) {
        this.parm1 = parm1;
        this.parm2 = parm2;
        this.zeroFraction = zeroFraction;
        this.distributionType = type;
        hasZeroFraction = true;
    }

    public Input<RealParameter> getParm1(){
        return this.parm1;
    }

    public Input<RealParameter> getMean(){
        if (this.distributionType == DistributionType.NORMAL || this.distributionType == DistributionType.LOGNORMAL)
            return this.parm1;
        else if (this.distributionType == DistributionType.GAMMA)
            return this.parm2;
            throw new IllegalArgumentException("This distribution is not parametrized with a mean.");
    }

    public Input<RealParameter> getShape(){
        if (this.distributionType == DistributionType.GAMMA)
            return this.parm1;
        else
            throw new IllegalArgumentException("This distribution is not parametrized with a shape.");
    }


    public Input<RealParameter> getStandardDev(){
        if (this.distributionType == DistributionType.NORMAL || this.distributionType == DistributionType.LOGNORMAL)
            return this.parm2;
        else
            throw new IllegalArgumentException("This distribution is not parametrized with a standard deviation.");
    }

    public Input<RealParameter> getZeroFraction() {
        return zeroFraction;
    }

    public DistributionType getDistributionType() {
        return distributionType;
    }


    public double getDensity(int cellType, double measuredValue){

        double result = 1;
        if (measuredValue == 0 && hasZeroFraction && this.zeroFraction.get() != null)
            return this.zeroFraction.get().getArrayValue(cellType);
        else {
            switch (distributionType) {
                case NORMAL:
                    result *= Utils.getNormalDensity(measuredValue, this.getMean().get().getArrayValue(cellType), this.getStandardDev().get().getArrayValue(cellType));
                    break;
                case LOGNORMAL:
                    result *= Utils.getLogNormalDensity(measuredValue, this.getMean().get().getArrayValue(cellType), this.getStandardDev().get().getArrayValue(cellType));
                    break;
                case GAMMA:
                    result *= Utils.getGammaDensityShapeMeanParam(measuredValue, this.getShape().get().getArrayValue(cellType), this.getMean().get().getArrayValue(cellType));
                    //TODO remove, to debug
//                    result *= Utils.getGammaDensityShapeRateParam(measuredValue, this.getShape().get().getArrayValue(cellType), this.getShape().get().getArrayValue(cellType)/this.getMean().get().getArrayValue(cellType));
//                    if(result > 10E30)
//                        System.out.printf("Approaching big values");
                    break;
                default:
                    throw new IllegalStateException("Unknown distribution type.");
            }

            if(hasZeroFraction && this.zeroFraction.get() != null && measuredValue != 0)
                result *= (1-this.zeroFraction.get().getArrayValue(cellType));
        }
        return result;
    }
}


