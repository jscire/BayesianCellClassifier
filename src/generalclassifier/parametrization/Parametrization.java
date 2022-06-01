package generalclassifier.parametrization;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import generalclassifier.lineagetree.Cell;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

//TODO rename class
public class Parametrization extends CalculationNode {

    public Input<List<DistributionForMeasurement>> distributionsInput = new Input<>("distribution",
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
                    "Note that all vectors in this list must have elements which sum up to 1." +
                    "If generation-specific transition probs are allowed " +
                    "(haveGenerationSpecificTransitionProbs set to TRUE)" +
                    "all first-generation transition probs are input first and then second-generation transition probs and so on.",
            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);

    public Input<List<RealParameter>> fateProbabilitiesInput = new Input<>("fateProbabilities",
            "List of 2-element vectors containing probabilities " +
                    "for each cell type of the cell being a divider or an apoptosing cell." +
                    "First element: prob of dividing, Second element: Prob of dying (apoptosis)" +
                    "Default: 1 for dividing, 0 for apoptosing for all types",
            new ArrayList<RealParameter>());

    public Input<RealParameter> typeFrequenciesInput = new Input<>("typeFrequencies",
            "Default: 1/numberOfCellTypes for each type.");

    public Input<RealParameter> lossProbInput = new Input<>("lossProb",
            "Probability of losing any given cell. Default: 0",
            new RealParameter("0"));

    public Input<Boolean> haveGenerationSpecificTransitionProbsInput = new Input<>("haveGenerationSpecificTransitionProbs",
            "Default: false", false);

    public Input<BooleanParameter> ignoreKinshipInfoInput = new Input<>("ignoreKinshipInfo",
            "If set to true, transition probabilities correspond to the probabilities of random draws" +
                    "for states distributed following type frequencies. Default: false", new BooleanParameter("false"));

    public int numberOfCellTypes;

    public final double ERRORMARGIN = 1e-3;

    SortedSet<String> uniqueMeasurementTags;

    @Override
    public void initAndValidate() {

        String tag;
        uniqueMeasurementTags = new TreeSet<>();

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
        // if generation-specific transition probs are allowed, check that number of vectors of transition probs is a multiple of the number of types
        if(haveGenerationSpecificTransitionProbsInput.get()) {
            if(transitionUponDivisionProbsInput.get().size() % numberOfCellTypes != 0)
                throw new IllegalArgumentException("Size of list of transitionUponDivisionProbs should " +
                        "correspond to a multiple of the number of cell types");
        } else {
            // check that number of vectors of transition probs is equal to the number of cell types
            if(transitionUponDivisionProbsInput.get().size() != numberOfCellTypes)
                throw new IllegalArgumentException("Size of list of transitionUponDivisionProbs should " +
                        "correspond to the number of cell types");
        }


        //TODO add check that all transitionUponDivisions elements sum up to 1 with an error margin.
        for (RealParameter param : transitionUponDivisionProbsInput.get()) {
            if(param.getDimension() != (numberOfCellTypes + 1) * numberOfCellTypes / 2.0)
                throw new IllegalArgumentException("Wrong number of elements in an element of transitionUponDivisionProbs." +
                        "There should be ((numberOfCellTypes + 1) * numberOfCellTypes / 2) elements");
        }


        if(fateProbabilitiesInput.get().size() == 0) {
            for (int i = 0; i < numberOfCellTypes; i++) {
                this.setInputValue("fateProbabilities", new RealParameter("1.0 0.0"));
            }
        }
        else {
            //TODO add check that all fateProbabilities elements sum up to 1 with an error margin.
            for(RealParameter param : fateProbabilitiesInput.get()) {
                if(param.getDimension() != 2)
                    throw new IllegalArgumentException("Wrong number of elements in element of fateProbabilitiesInput." +
                            "There should be exactly two elements.");
            }
        }

        if(typeFrequenciesInput.get() == null) {
            double defaultFrequency = 1.0/numberOfCellTypes;
            String freqInputVal = "";
            for (int i = 0; i < numberOfCellTypes; i++) {
                freqInputVal += (defaultFrequency + " ");
            }
            this.setInputValue("typeFrequencies", new RealParameter(freqInputVal));
        }
        else if (typeFrequenciesInput.get().getDimension() != numberOfCellTypes) {
            throw new IllegalArgumentException("Wrong number of dimensions in typeFrequenciesInput.");
        }
        else {
            double totalFreq = 0;
            for (int i = 0; i < numberOfCellTypes; i++) {
                totalFreq += typeFrequenciesInput.get().getArrayValue(i);
            }
            if(Math.abs(1-totalFreq) > ERRORMARGIN)
                throw new IllegalArgumentException("Type frequencies do not sum to 1.");
        }

        if(lossProbInput.get() != null && lossProbInput.get().getDimension() != 1)
            throw new IllegalArgumentException("lossProb must be of dimension 1.");
    }

    public List<DistributionForMeasurement> getDistributions() {
        return distributionsInput.get();
    }

    /**
     * Returns the value in transitionUponDivisionProbsInput corresponding to the type triplet.
     * Adds a 0.5 factor if type of child 1 is different to child 2
     * @param typeMother
     * @param typeChild1
     * @param typeChild2
     * @return
     */
    public double getTransitionProbability(int typeMother, int typeChild1, int typeChild2, int generationMother, boolean isTreeOfKnowType){

        if(ignoreKinshipInfoInput.get().getValue()) {
            // if kinship is ignored the transition probability is just the probability of randomly drawing
            // two cells of type typeChild1 and typeChild2 from the overall population

            // in the case where the entire tree is of know type, this transition prob is fixed to 1.
            // (The frequency of the known type is 1 in that case)
            if(isTreeOfKnowType)
                return 1;
            else
                return(typeFrequenciesInput.get().getArrayValue(typeChild1) * typeFrequenciesInput.get().getArrayValue(typeChild2));
        }

        if(typeChild1 > typeChild2) {
            int temp = typeChild2;
            typeChild2 = typeChild1;
            typeChild1 = temp;
        }

        int indexInFlattenedArray = numberOfCellTypes * typeChild1 + typeChild2 - typeChild1 * (typeChild1+1)/2;

        double valueInputArray;
        if(haveGenerationSpecificTransitionProbsInput.get()) {
            // we get the (typeMother + numberOfCellTypes * (generationMother-1))th vector of transition probabilities.
            // we take (generationMother-1) because the root cell is, by convention, at generation number 1.
            valueInputArray = transitionUponDivisionProbsInput.get().get(typeMother + numberOfCellTypes * (generationMother-1) ).getArrayValue(indexInFlattenedArray);
        } else {
            valueInputArray = transitionUponDivisionProbsInput.get().get(typeMother).getArrayValue(indexInFlattenedArray);
        }

        if(typeChild1 != typeChild2)
            return valueInputArray * 0.5; // 1/2 factor because 2 possibilities with child1 type j and child1 type k.
        else
            return valueInputArray;

    }
    
    public double[] getTransitionProbabilitiesForMotherType(int typeMother){

        double[] probs = new double[transitionUponDivisionProbsInput.get().get(typeMother).getDimension()];
        for (int i = 0; i < probs.length; i++) {
            probs[i] = transitionUponDivisionProbsInput.get().get(typeMother).getArrayValue(i);
        }
        return probs;
    }

    public SortedSet<String> getMeasurementTags(){
        return uniqueMeasurementTags;
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

    public double[] getFateProbability(int cellType) {
        double[] probs = new double[fateProbabilitiesInput.get().get(cellType).getDimension()];
        for (int i = 0; i < probs.length; i++) {
            probs[i] = fateProbabilitiesInput.get().get(cellType).getArrayValue(i);
        }
        return probs;
    }

    public double getTypeFreq(int type){
        return typeFrequenciesInput.get().getArrayValue(type);
    }

    public int getNumberOfCellTypes() {
        return numberOfCellTypes;
    }
    
    public static void main(String[] args){
        String tag1 = "lifetime";

        DistributionForMeasurement distr1 = new DistributionForMeasurement();

        distr1.initByName("measurementTag", tag1,
                "parm1Distribution", new RealParameter("1.0 2.3"),
                "parm2Distribution", new RealParameter("0.5 0.3"),
                "zeroFraction", new RealParameter("0.01 0.01"),
                "distributionType", "lognormal",
                "estimateType", "max",
                "isAppliedToRootCells", true);

        distr1.initAndValidate();

        String tag2 = "Sca1";

        DistributionForMeasurement distr2 = new DistributionForMeasurement();

        distr2.initByName("measurementTag", tag2,
                "parm1Distribution", new RealParameter("1.0 2.3"),
                "parm2Distribution", new RealParameter("0.5 0.3"),
                "zeroFraction", new RealParameter("0.01 0.01"),
                "distributionType", "lognormal",
                "estimateType", "max",
                "isAppliedToRootCells", true);

        distr2.initAndValidate();

        Parametrization parametrization = new Parametrization();

        parametrization.initByName("distribution", distr1,
                "distribution", distr2,
                "transitionUponDivisionProbs", new RealParameter("0.2 0.25 0.55"),
                "transitionUponDivisionProbs", new RealParameter("0.3 0.1 0.6"));

        System.out.println("Done");


    }
}
