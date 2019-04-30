package generalclassifier.parametrization;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import generalclassifier.lineagetree.Cell;
import generalclassifier.utils.Utils;

//TODO change name of class
public class DistributionForMeasurement extends CalculationNode {

    public enum DistributionType{
        WEIBULL_MEDIAN_SHAPE,
        WEIBULL_SCALE_SHAPE,
        BETA,
        GAMMA_MEAN_SHAPE,
        LOGNORMAL,
        NORMAL
    }

    public enum EstimateType {
        MIN,
        MAX,
        MEAN,
    }


    public Input<RealParameter> parm1DistributionInput = new Input<>("parm1Distribution",
            "First parameter of the distribution.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> parm2DistributionInput = new Input<>("parm2Distribution",
            "Second parameter of the distribution.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> zeroFractionInput = new Input<>("zeroFraction", "Optional. " +
            "Probability mass associated with point value 0. " +
            "Useful when using distributions that exclude negative values, but measurements show negative values. " +
            "These negative values can be approximated as being 0, and a zeroFraction specified.");

    public Input<String> measurementTagInput = new Input<>("measurementTag",
            "Tag used for tracking what phenotype measure we are dealing with. " +
            "Logger will use this measurementTag as well.",
            Input.Validate.REQUIRED);

    public Input<String> distributionTypeInput = new Input<>("distributionType", "One of: " +
            "'gamma', 'lognormal', 'normal', 'beta', 'weibull_scale_shape', 'weibull_median_shape'." +
            "Default: normal",
            "normal",
            new String[]{"gamma", "lognormal", "normal", "beta", "weibull_scale_shape", "weibull_median_shape"});

    public Input<String> estimateTypeInput = new Input<>("estimateType", "One of 'mean', 'max', 'min'. " +
            "Type of summary taken on the time series of measurements. " +
            "Needed to know how to deal with measures on cells whose life has not entirely been recorded. " +
            "Useful for root cells for instance. Default: 'mean'",
            "mean",
            new String[]{"mean", "max", "min"});

    public Input<Boolean> isAppliedOnRootCellsInput = new Input<>("isAppliedOnRootCells",
            "Default: true",
            true);

    DistributionType distributionType;

    EstimateType estimateType;

    String measurementTag;

    boolean hasZeroFraction;

    boolean isAppliedOnRootCells;

    public int numberOfCellTypes;

    @Override
    public void initAndValidate() {
        if(parm1DistributionInput.get().getDimension() != parm2DistributionInput.get().getDimension())
            throw new IllegalArgumentException("Invalid number of dimensions in parm1Distribution or parm2Distribution." +
                    "The number of dimensions should be equal for these two input and be equal to the number of cell types.");

        numberOfCellTypes = parm1DistributionInput.get().getDimension();

        hasZeroFraction = zeroFractionInput.get() != null;

        if(hasZeroFraction && zeroFractionInput.get().getDimension() != numberOfCellTypes)
            throw new IllegalArgumentException("Invalid number of dimemsions in zero fraction." +
                    " It should correspond to the number of cell types.");

        switch(estimateTypeInput.get()) {
            case "mean":
                estimateType = EstimateType.MEAN;
                break;
            case "max":
                estimateType = EstimateType.MAX;
                break;
            case "min":
                estimateType = EstimateType.MIN;
                break;
            default:
                throw new IllegalArgumentException("Unknown estimate type. Input must be 'mean', 'max' or 'min'.");
        }

        switch(distributionTypeInput.get()){
            case "gamma":
                distributionType = DistributionType.GAMMA_MEAN_SHAPE;
                break;
            case "lognormal":
                distributionType = DistributionType.LOGNORMAL;
                break;
            case "normal":
                distributionType = DistributionType.NORMAL;
                break;
            case "beta":
                distributionType = DistributionType.BETA;
                break;
            case "weibull_scale_shape":
                distributionType = DistributionType.WEIBULL_SCALE_SHAPE;
                break;
            case "weibull_median_shape":
                distributionType = DistributionType.WEIBULL_MEDIAN_SHAPE;
                break;
            default:
                throw new IllegalArgumentException("Unknown distribution type." +
                        "Input must be one of: 'gamma', 'lognormal', 'normal', 'beta', " +
                        "'weibull_scale_shape', 'weibull_median_shape'");
        }

        measurementTag = measurementTagInput.get();

        isAppliedOnRootCells = isAppliedOnRootCellsInput.get();
    }


    /**
     * For now, value is independent of cell fate (does not matter if divides or dies, even for lifetime)
     * However, we do take into account the added uncertainty when the end of the cell life is not observed (Fate.U)
     * @param measuredValue
     * @param cellType
     * @param isIncompleteObservation
     * @param cellFate
     * @return
     */
    public double getProbability(double measuredValue, int cellType, boolean isIncompleteObservation, Cell.Fate cellFate) {

        if(cellType >= numberOfCellTypes)
            throw new IllegalArgumentException("Incorrect cell type. Cell type must be >= 0 and < numberOfCellTypes");

        switch(estimateType) {
            case MEAN:
                return getProbabilityDensity(measuredValue, cellType);
            case MAX:
                if(isIncompleteObservation || cellFate == Cell.Fate.U)
                    return getOppositeCumulativeDistribution(measuredValue, cellType);
                else
                    return getProbabilityDensity(measuredValue, cellType);
            case MIN:
                if(isIncompleteObservation || cellFate == Cell.Fate.U)
                    return getCumulativeDistribution(measuredValue, cellType);
                else
                    return getProbabilityDensity(measuredValue, cellType);
            default:
                throw new IllegalArgumentException("Unknown estimate type.");
        }

    }

    private double getProbabilityDensity(double measuredValue, int cellType){

        double p = 1;

        if(measuredValue == 0 && hasZeroFraction) return zeroFractionInput.get().getArrayValue(cellType);

        switch(distributionType) {
            case WEIBULL_MEDIAN_SHAPE:
                p *= Utils.getWeibullDensityMedianShapeParam(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case WEIBULL_SCALE_SHAPE:
                p *= Utils.getWeibullDensity(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case BETA:
                p *= Utils.getBetaDensity(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case GAMMA_MEAN_SHAPE:
                p *= Utils.getGammaDensityMeanShapeParam(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case LOGNORMAL:
                p *= Utils.getLogNormalDensity(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case NORMAL:
                p *= Utils.getNormalDensity(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            default:
                throw new IllegalArgumentException("Unknown distribution type");
        }

        if(hasZeroFraction) p *= (1 - zeroFractionInput.get().getArrayValue(cellType));

        return p;
    }

    private double getCumulativeDistribution(double measuredValue, int cellType) {

        double p = 0;

        switch(distributionType) {
            case WEIBULL_MEDIAN_SHAPE:
                p += Utils.getWeibullCumulativeDistributionMedianShapeParam(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case WEIBULL_SCALE_SHAPE:
                p += Utils.getWeibullCumulativeDistribution(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case BETA:
                p += Utils.getBetaCumulativeDistribution(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case GAMMA_MEAN_SHAPE:
                p += Utils.getGammaCumulativeDistributionMeanShapeParam(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case LOGNORMAL:
                p += Utils.getLogNormalCumulativeDistribution(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case NORMAL:
                p += Utils.getNormalCumulativeDistribution(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            default:
                throw new IllegalArgumentException("Unknown distribution type");

        }

        if(hasZeroFraction) {
            p /= 1 + zeroFractionInput.get().getArrayValue(cellType); // normalize to 1.

            if (measuredValue >= 0)
                p += zeroFractionInput.get().getArrayValue(cellType);
        }

        return p;
    }

    private double getOppositeCumulativeDistribution(double measuredValue, int cellType) {

        double p = 0;

        switch(distributionType) {
            case WEIBULL_MEDIAN_SHAPE:
                p += 1- Utils.getWeibullCumulativeDistributionMedianShapeParam(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case WEIBULL_SCALE_SHAPE:
                p += 1- Utils.getWeibullCumulativeDistribution(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case BETA:
                p += 1- Utils.getBetaCumulativeDistribution(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case GAMMA_MEAN_SHAPE:
                p += 1- Utils.getGammaCumulativeDistributionMeanShapeParam(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case LOGNORMAL:
                p += 1- Utils.getLogNormalCumulativeDistribution(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            case NORMAL:
                p += 1- Utils.getNormalCumulativeDistribution(measuredValue,
                        parm1DistributionInput.get().getArrayValue(cellType),
                        parm2DistributionInput.get().getArrayValue(cellType));
                break;
            default:
                throw new IllegalArgumentException("Unknown distribution type");

        }

        if(hasZeroFraction) {
            p /= 1 + zeroFractionInput.get().getArrayValue(cellType); // normalize to 1.

            if (measuredValue == 0)
                p += zeroFractionInput.get().getArrayValue(cellType);
        }

        return p;
    }

    public String getMeasurementTag(){
        return measurementTag;
    }

    public int getNumberOfCellTypes(){
        return numberOfCellTypes;
    }

    public boolean isAppliedOnRootCells(){
        return isAppliedOnRootCells;
    }

}
