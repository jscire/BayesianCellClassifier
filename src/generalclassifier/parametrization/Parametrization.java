package generalclassifier.parametrization;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import generalclassifier.lineagetree.Cell;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

//TODO rename class
public class Parametrization extends CalculationNode {

    public Input<List<DistributionForMeasurement>> distributionsInput = new Input<>("distributions",
            "List of DistributionForMeasurement," +
                    "containing all the types of phenotypic measurements taken into account " +
                    "and the characteristics of the distribution they are modeled by.",
            new ArrayList<DistributionForMeasurement>(), Input.Validate.REQUIRED);

    public Input<List<RealParameter>> transitionUponDivisionProbsInput = new Input<>("transitionUponDivisionProbs",
            "List of flattened 2-dimensional arrays containing the probabilites of transition between states upon division." +
                    "If mkj in array i represents the probability that a cell of type i divides into cells of type k and j (k<= j)," +
                    "the element mkj has index n*k + j - k(k+1)/2 in the input vector." +
                    "If j<k, invert k and j." +
                    "We assume here that the daughter cells are not ordered and " +
                    "so that the arrays are symmetric before being truncated (the triangle below the diagonal is removed)." +
                    "This renumbering is equivalent to flattening this truncated array by rows." +
                    "Note that all vectors in this list must have elements which sum up to 1.",
            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);

    public Input<List<RealParameter>> fateProbabilitiesInput = new Input<>("fateProbabilities",
            "List of 2-element vectors containing probabilities " +
                    "for each cell type of the cell being a divider or an apoptosing cell." +
                    "First element: prob of dividing, Second element: Prob of dying (apoptosis)" +
                    "Default: 1 for dividing, 0 for apoptosing for all types",
            new ArrayList<RealParameter>());

    public Input<RealParameter> lossProbInput = new Input<>("lossProb",
            "Probability of losing any given cell. Default: 0",
            new RealParameter("0"));

    public int numberOfCellTypes;

    HashSet<String> uniqueMeasurementTags;

    @Override
    public void initAndValidate() {

        String tag;

        numberOfCellTypes = -1;

        for(DistributionForMeasurement distr : distributionsInput.get()) {

            tag = distr.getMeasurementTag();

            //check for duplication in measurement tags.
            if(uniqueMeasurementTags.contains(tag))
                throw new IllegalArgumentException("Duplicated tag: " + tag + " in list of DistributionForMeasurement.");
            else
                uniqueMeasurementTags.add(tag);

            //check that all distribution inputs specify the same number of cell types.
            if(numberOfCellTypes == -1)
                numberOfCellTypes = distr.getNumberOfCellTypes();
            else if(numberOfCellTypes != distr.getNumberOfCellTypes())
                throw new IllegalArgumentException("All distributions must have the same number of cell types");
        }
        
        if(transitionUponDivisionProbsInput.get().size() != numberOfCellTypes)
            throw new IllegalArgumentException("Size of list of transitionUponDivisionProbs should " +
                    "correspond to the number of cell types");

        for (RealParameter param : transitionUponDivisionProbsInput.get()) {
            if(param.getDimension() != (numberOfCellTypes - 1) * numberOfCellTypes / 2.0)
                throw new IllegalArgumentException("Wrong number of elements in an element of transitionUponDivisionProbs." +
                        "There should be ((numberOfCellTypes - 1) * numberOfCellTypes / 2) elements");
        }


        if(fateProbabilitiesInput.get() == null) {
            List<RealParameter> fateProbabilities = new ArrayList<RealParameter>();
            for (int i = 0; i < numberOfCellTypes; i++) {
                fateProbabilities.add(new RealParameter("1.0 0.0"));
            }
            fateProbabilitiesInput.set(fateProbabilities);
        }
        else {
            for(RealParameter param : fateProbabilitiesInput.get()) {
                if(param.getDimension() != 2)
                    throw new IllegalArgumentException("Wrong number of elements in element of fateProbabilitiesInput." +
                            "There should be exactly two elements.");
            }
        }

        if(lossProbInput.get() != null && lossProbInput.get().getDimension() != 1)
            throw new IllegalArgumentException("lossProb must be of dimension 1.");
    }

    public List<DistributionForMeasurement> getDistributions() {
        return distributionsInput.get();
    }

    public double getTransitionProbability(int typeParent, int typeChild1, int typeChild2){

        if(typeChild1 > typeChild2) {
            int temp = typeChild2;
            typeChild2 = typeChild1;
            typeChild1 = temp;
        }

        int indexInFlattenedArray = numberOfCellTypes * typeChild1 + typeChild2 - typeChild1 * (typeChild1+1)/2;

        return transitionUponDivisionProbsInput.get().get(typeParent).getArrayValue(indexInFlattenedArray);
    }

    public double getLossProbability() {
        return lossProbInput.get().getValue();
    }

    public double getFateProbability(Cell.Fate fate, int cellType) {

        double lossProbability = lossProbInput.get().getValue();

        if(fate == Cell.Fate.L)
            return lossProbability;

        if(fate == Cell.Fate.U)
            return 1 - lossProbability;

        int idxFate = fate == Cell.Fate.D ? 0 : 1;

        return fateProbabilitiesInput.get().get(cellType).getArrayValue(idxFate) * (1 - lossProbability);
    }
}
