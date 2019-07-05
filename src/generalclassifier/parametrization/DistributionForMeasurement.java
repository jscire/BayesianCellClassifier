package generalclassifier.parametrization;

import beast.core.CalculationNode;
import beast.core.Distribution;
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

    public Input<Boolean> isAppliedOnRootCellsInput = new Input<>("isAppliedToRootCells",
            "Default: true",
            true);


    DistributionType distributionType;

    EstimateType estimateType;

    String measurementTag;

    boolean hasZeroFraction = false;

    boolean isAppliedToRootCells;

    public int numberOfCellTypes;

    @Override
    public void initAndValidate() {
        if(parm1DistributionInput.get().getDimension() != parm2DistributionInput.get().getDimension())
            throw new IllegalArgumentException("Invalid number of dimensions in parm1Distribution or parm2Distribution." +
                    "The number of dimensions should be equal for these two input and be equal to the number of cell types.");

        numberOfCellTypes = parm1DistributionInput.get().getDimension();


        // set hasZeroFraction to true if input zeroFraction is not null and if at least one of the values is not zero.
        // we assume that only scaling operators and that the zero fractions will stay 0 if they start at 0.
        if(zeroFractionInput.get() != null) {
            for (int i = 0; i < zeroFractionInput.get().getDimension(); i++) {
                if(zeroFractionInput.get().getArrayValue(i) != 0){
                    hasZeroFraction = true;
                    break;
                }
            }
        }

        if(hasZeroFraction && zeroFractionInput.get().getDimension() != numberOfCellTypes)
            throw new IllegalArgumentException("Invalid number of dimensions in zero fraction." +
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

        isAppliedToRootCells = isAppliedOnRootCellsInput.get();
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

            if (measuredValue >= 0)
                p += zeroFractionInput.get().getArrayValue(cellType);

            p /= (1 + zeroFractionInput.get().getArrayValue(cellType)); // normalize to 1.
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
            if (measuredValue == 0)
                p += zeroFractionInput.get().getArrayValue(cellType);

            p /= (1 + zeroFractionInput.get().getArrayValue(cellType)); // normalize to 1.
        }

        return p;
    }

    public String getMeasurementTag(){
        return measurementTag;
    }

    public DistributionType getDistributionType() {
        return distributionType;
    }

    public boolean getHasZeroFraction(){
        return hasZeroFraction;
    }

    public double getZeroFraction(){
        if (!hasZeroFraction) return 0.0;
        return zeroFractionInput.get().getValue();
    }

    public double getParm1(int cellType){
        if(cellType >= numberOfCellTypes)
            throw new IllegalArgumentException("Invalid cell type (higher than number of cell types.");
        return parm1DistributionInput.get().getArrayValue(cellType);
    }

    public double getParm2(int cellType){
        if(cellType >= numberOfCellTypes)
            throw new IllegalArgumentException("Invalid cell type (higher than number of cell types.");
        return parm2DistributionInput.get().getArrayValue(cellType);
    }

    public int getNumberOfCellTypes(){
        return numberOfCellTypes;
    }

    public boolean isAppliedToRootCells(){
        return isAppliedToRootCells;
    }

    public static void main(String[] args) {
        String tag = "lifetime";

        DistributionForMeasurement distr = new DistributionForMeasurement();

        distr.initByName("measurementTag", tag,
        "parm1Distribution", new RealParameter("1.0 2.3"),
        "parm2Distribution", new RealParameter("0.5 0.3"),
        "zeroFraction", new RealParameter("0.01 0.01"),
        "distributionType", "lognormal",
        "estimateType", "max",
        "isAppliedToRootCells", true);

        distr.initAndValidate();

        double p = distr.getProbability(3.8, 0, true, Cell.Fate.U);
        System.out.println("Prob: " + p);
    }

}
