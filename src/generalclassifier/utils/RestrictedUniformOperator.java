package generalclassifier.utils;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.util.Randomizer;


/**
 * Code adapted from UniformOperator and ScaleOperator classes from BEAST 2 .
 * Copied UniformOperator instead of extended because of accessibility restrictions on class fields.
 */

@Description("Assign one or more parameter values to a uniformly selected value in its range. " +
        "Lets you exclude some dimensions from the selection.")
public class RestrictedUniformOperator extends Operator {

    final public Input<IntegerParameter> parameterInput = new Input<>("parameter", "an integer parameter to sample individual values for", Input.Validate.REQUIRED);
    final public Input<BooleanParameter> indicatorInput = new Input<>("indicator",
            "indicates which of the dimensions of the parameters can be changed.",
            Input.Validate.REQUIRED);

    IntegerParameter parameter;
    int lowerIndex, upperIndex;

    @Override
    public void initAndValidate() {
        parameter = parameterInput.get();

        lowerIndex = parameter.getLower();
        upperIndex = parameter.getUpper();

        final BooleanParameter indicators = indicatorInput.get();
        if (indicators != null) {
            final int dataDim = parameterInput.get().getDimension();
            final int indsDim = indicators.getDimension();
            if (!(indsDim == dataDim || indsDim + 1 == dataDim)) {
                throw new IllegalArgumentException("indicator dimension not compatible from parameter dimension");
            }
        }
    }

    @Override
    public double proposal() {

        // which position to scale
        final int index;
        final int dim = parameter.getDimension();

        final BooleanParameter indicators = indicatorInput.get();
        if (indicators != null) {
            final int dimCount = indicators.getDimension();
            final Boolean[] indicator = indicators.getValues();
            final boolean impliedOne = dimCount == (dim - 1);

            // available bit locations. there can be hundreds of them. scan list only once.
            final int[] loc = new int[dimCount + 1];
            int locIndex = 0;

            if (impliedOne) {
                loc[locIndex] = 0;
                ++locIndex;
            }
            for (int i = 0; i < dimCount; i++) {
                if (indicator[i]) {
                    loc[locIndex] = i + (impliedOne ? 1 : 0);
                    ++locIndex;
                }
            }

            if (locIndex > 0) {
                final int rand = Randomizer.nextInt(locIndex);
                index = loc[rand];
            } else {
                return Double.NEGATIVE_INFINITY; // no active indicators
            }

        } else {
            index = Randomizer.nextInt(parameter.getDimension());
        }

        int newValue = Randomizer.nextInt(upperIndex - lowerIndex + 1) + lowerIndex; // from 0 to n-1, n must > 0,
        parameter.setValue(index, newValue);

        return 0.0;
    }

}
