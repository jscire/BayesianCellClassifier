package generalclassifier.core;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

import generalclassifier.lineagetree.Cell;
import generalclassifier.lineagetree.LineageTree;
import generalclassifier.lineagetree.MeasureType;
import generalclassifier.utils.InputPair;
import generalclassifier.utils.Utils;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

import java.io.*;
import java.util.*;


public class  LineageTreeProb extends Distribution {

    public Input<LineageTree> lineageTreeInput = new Input<>("tree",
            "Lineage tree.",
            Input.Validate.REQUIRED);

    public Input<List<RealParameter>> fateProbabilitiesInput = new Input<>("fateProbabilities",
            "List of vectors containing probabilities for each cell type of the cell being a divider, apoptoser, non-divider or type-transitioner (if defined).",
            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);

    public Input<List<RealParameter>> scaleWeibullInput = new Input<>("scaleWeibull",
            "List of vectors containing Weibull scale parameters for each type for division, apoptosis and type-transition (if defined) times.",
            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);

    public Input<List<RealParameter>> shapeWeibullInput = new Input<>("shapeWeibull",
            "List of vectors containing Weibull shape parameters for each type for division, apoptosis and type-transition (if defined) times.",
            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);

    public Input<List<RealParameter>> transitionUponDivisionProbsInput = new Input<>("transitionUponDivisionProbs",
            "List of flattened 2-dimensional arrays containing the probabilites of transition between states upon division." +
                    "If mkj in array i represents the probability that a cell of type i divides into cells of type k and j (k<= j)," +
                    "the element of index n*k + j - k(k+1)/2 gives the index of element mkj in the input vector." +
                    "If j<k, invert k and j. We assume here that the daughter cells are not ordered and so that the arrays are symmetric before being truncated." +
                    "This renumbering is equivalent to flattening the array by rows." +
                    "Note that all vectors in this list must have elements which sum up to 1.",
            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);


    public Input<RealParameter> lossProbInput = new Input<>("lossProb",
            "Probability of losing any given cell",
            Input.Validate.REQUIRED);

    public Input<BooleanParameter> rootIsHSCInput = new Input<>("rootIsHSC",
            "True if root is an HSC cell, false otherwise.",
            Input.Validate.REQUIRED);

    public Input<Integer> rootIsHSCidxInput = new Input<>("rootIsHSCidx",
            "If provided, this is an index into the rootIsHSC boolean " +
                    "parameter specifying which element corresponds to the root " +
                    "HSC status.");

    public Input<BooleanParameter> matrixOfAllowedTransitionsInput = new Input<> ("allowedTransitions",
            "Flattened matrix of allowed transitions between cell types. " +
                    "If mij is the boolean that characterizes transitions from type i to j, the element of index i*(n-1)+j if i>j, or i*(n-1) + j-1 otherwise, gives the index of element mij. n is the number of types.");

    // TODO possibly rework on the way inputs are written in. This is a first draft
    public Input<RealParameter> meanNormalAreaGrowthRateInput = new Input<>("meanAreaGrowthRate",
            "");
    public Input<RealParameter> sdNormalAreaGrowthRateInput = new Input<>("sdAreaGrowthRate",
            "");
    public Input<RealParameter> meanNormalPerimeterGrowthRateInput = new Input<>("meanPerimeterGrowthRate",
            "");
    public Input<RealParameter> sdNormalPerimeterGrowthRateInput = new Input<>("sdPerimeterGrowthRate",
            "");
    public Input<RealParameter> meanNormalEccentricityInput = new Input<>("meanEccentricity",
            "");
    public Input<RealParameter> sdNormalEccentricityInput = new Input<>("sdEccentricity",
            "");
    public Input<RealParameter> meanNormalInstantSpeedInput = new Input<>("meanInstantSpeed",
            "");
    public Input<RealParameter> sdNormalInstantSpeedInput = new Input<>("sdInstantSpeed",
            "");
    public Input<RealParameter> meanNormalCD41ProductionRateInput = new Input<>("meanCD41ProductionRate",
            "");
    public Input<RealParameter> sdNormalCD41ProductionRateInput = new Input<>("sdCD41ProductionRate",
            "");
    public Input<RealParameter> meanNormalFcgRIIIProductionRateInput = new Input<>("meanFcgRIIIProductionRate",
            "");
    public Input<RealParameter> sdNormalFcgRIIIProductionRateInput = new Input<>("sdFcgRIIIProductionRate",
            "");
    public Input<RealParameter> meanNormalROSProductionRateInput = new Input<>("meanROSProductionRate",
            "");
    public Input<RealParameter> sdNormalROSProductionRateInput = new Input<>("sdROSProductionRate",
            "");
    public Input<RealParameter> meanNormalTMRMProductionRateInput = new Input<>("meanTMRMProductionRate",
            "");
    public Input<RealParameter> sdNormalTMRMProductionRateInput = new Input<>("sdTMRMProductionRate",
            "");
    public Input<RealParameter> meanNormalTMRMMaxRateInput = new Input<>("meanTMRMMaxRate",
            "");
    public Input<RealParameter> sdNormalTMRMMaxRateInput = new Input<>("sdTMRMMaxRate",
            "");
    public Input<RealParameter> meanNormalSca1ProductionRateInput = new Input<>("meanSca1ProductionRate",
            "");
    public Input<RealParameter> sdNormalSca1ProductionRateInput = new Input<>("sdSca1ProductionRate",
            "");
    public Input<RealParameter> meanNormalIg2afcProductionRateInput = new Input<>("meanIg2afcProductionRate",
            "");
    public Input<RealParameter> sdNormalIg2afcProductionRateInput = new Input<>("sdIg2afcProductionRate",
            "");
    public Input<RealParameter> meanNormalCD71APCProductionRateInput = new Input<>("meanCD71APCProductionRate",
            "");
    public Input<RealParameter> sdNormalCD71APCProductionRateInput = new Input<>("sdCD71APCProductionRate",
            "");
    public Input<RealParameter> meanNormalCD71PEProductionRateInput = new Input<>("meanCD71PEProductionRate",
            "");
    public Input<RealParameter> sdNormalCD71PEProductionRateInput = new Input<>("sdCD71PEProductionRate",
            "");
    public Input<RealParameter> meanNormalcMycGFPMaxRateInput = new Input<>("meancMycGFPMaxRate",
            "");
    public Input<RealParameter> sdNormalcMycGFPMaxRateInput = new Input<>("sdcMycGFPMaxRate",
            "");


    Map<InputPair, MeasureType> mapMeasureTypeToInput = createMapOfMeasureInputs();

    private Map<InputPair, MeasureType> createMapOfMeasureInputs(){
        Map<InputPair, MeasureType> result = new HashMap<>();

        InputPair areaGrowthRateInput = new InputPair(meanNormalAreaGrowthRateInput, sdNormalAreaGrowthRateInput);
        InputPair perimeterGrowthRateInput = new InputPair(meanNormalPerimeterGrowthRateInput, sdNormalPerimeterGrowthRateInput);
        InputPair eccentricityInput = new InputPair(meanNormalEccentricityInput, sdNormalEccentricityInput);
        InputPair instantSpeedInput = new InputPair(meanNormalInstantSpeedInput, sdNormalInstantSpeedInput);
        InputPair CD41ProductionRateInput = new InputPair(meanNormalCD41ProductionRateInput, sdNormalCD41ProductionRateInput);
        InputPair FcgRIIIProductionRateInput = new InputPair(meanNormalFcgRIIIProductionRateInput, sdNormalFcgRIIIProductionRateInput);
        InputPair ROSProductionRateInput = new InputPair(meanNormalROSProductionRateInput, sdNormalROSProductionRateInput);
        InputPair TMRMProductionRateInput = new InputPair(meanNormalTMRMProductionRateInput, sdNormalTMRMProductionRateInput);
        InputPair TMRMMaxRateInput = new InputPair(meanNormalTMRMMaxRateInput, sdNormalTMRMMaxRateInput);
        InputPair Sca1ProductionRateInput = new InputPair(meanNormalSca1ProductionRateInput, sdNormalSca1ProductionRateInput);
        InputPair Ig2afcProductionRateInput = new InputPair(meanNormalIg2afcProductionRateInput, sdNormalIg2afcProductionRateInput);
        InputPair CD71APCProductionRateInput = new InputPair(meanNormalCD71APCProductionRateInput, sdNormalCD71APCProductionRateInput);
        InputPair CD71PEProductionRateInput = new InputPair(meanNormalCD71PEProductionRateInput, sdNormalCD71PEProductionRateInput);
        InputPair cMycGFPMaxRateInput = new InputPair(meanNormalcMycGFPMaxRateInput, sdNormalcMycGFPMaxRateInput);
        

        result.put(areaGrowthRateInput, MeasureType.Area);
        result.put(perimeterGrowthRateInput, MeasureType.Perimeter);
        result.put(eccentricityInput, MeasureType.Eccentricity);
        result.put(instantSpeedInput, MeasureType.InstantSpeed);
        result.put(CD41ProductionRateInput, MeasureType.CD41);
        result.put(FcgRIIIProductionRateInput, MeasureType.FcgRIII);
        result.put(ROSProductionRateInput, MeasureType.ROS);
        result.put(TMRMProductionRateInput, MeasureType.TMRMmean);
        result.put(TMRMMaxRateInput, MeasureType.TMRMmax);
        result.put(Sca1ProductionRateInput, MeasureType.Sca1);
        result.put(Ig2afcProductionRateInput, MeasureType.Ig2afc);
        result.put(CD71APCProductionRateInput, MeasureType.CD71APC);
        result.put(CD71PEProductionRateInput, MeasureType.CD71PE);
        result.put(cMycGFPMaxRateInput, MeasureType.cMycGFP);

        return Collections.unmodifiableMap(result);
    }

    LineageTree tree;

    boolean rootIsHSC;
    //TODO allow for specifying the types of leaves or even other cells of the tree.

    boolean transitionUponDivisionIsAllowed;
    boolean transitionDuringLifetimeIsAllowed;

    int numberOfTypes;

    // TrapezoidIntegrator numericalIntegrator = new TrapezoidIntegrator(1e-6, 1e-6, 1, 60);
    IterativeLegendreGaussIntegrator numericalIntegrator = new IterativeLegendreGaussIntegrator(2, 1e-6, 1e-6);

    int maxEval = 10000;


    //TODO turn every step in the calculation into log calculations
    //TODO change code so that there is a matrix of transition probabilities for transitions on branches

    public LineageTreeProb() {
    }

    @Override
    public void initAndValidate() {

        tree = lineageTreeInput.get();

        numberOfTypes = fateProbabilitiesInput.get().size();

        //TODO remove when tested with more than 2 types.
        if(numberOfTypes != 2)
            throw new IllegalStateException("Number of types is not two. Such a configuration is not tested yet.");

        if(shapeWeibullInput.get().size() != numberOfTypes || scaleWeibullInput.get().size() != numberOfTypes)
            throw new IllegalStateException("Inputs fateProbabilities, shapeWeibull and scaleWeibull should have the same number of elements: the number of cell types.");

        if(matrixOfAllowedTransitionsInput.get() != null) {
            transitionDuringLifetimeIsAllowed = true;
            throw new IllegalStateException("Transitions on branches are not fully implemented yet");

            //TODO uncomment check below when exception above is not thrown anymore (when transitions on branches are implemented)
//            if(matrixOfAllowedTransitionsInput.get().getDimension() != numberOfTypes * (numberOfTypes -1))
//                throw new IllegalStateException("Incorrect size of matrix of allowed transitions.");
        }
        else {
            transitionDuringLifetimeIsAllowed = false;
        }

        if(transitionUponDivisionProbsInput.get().get(0) != null) {
            transitionUponDivisionIsAllowed=true;

            if(transitionUponDivisionProbsInput.get().size() != numberOfTypes)
                throw new IllegalStateException("Incorrect size of list of probabilities of transition upon division.");
        }
        else {
            transitionUponDivisionIsAllowed = false;
        }

        if (lossProbInput.get() == null)
            throw new IllegalStateException("lossProb input is missing. Set to zero if cells are never lost.");

        //TODO add checks for parameter input format. In particular, check that fateProbabilities, transitionProbsUponDivision, shapeWeibull and scaleWeibull have correct dimension for each cell type.

    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        if (rootIsHSCidxInput.get() != null)
            rootIsHSC = rootIsHSCInput.get().getValue(rootIsHSCidxInput.get());
        else
            rootIsHSC = rootIsHSCInput.get().getValue();

        //TODO rework on this when adding more types
        int rootNodeType = rootIsHSC? 0 : 1;

        //for now prior is that frequency is equal for all types TODO allow for more flexibility
        try {
            logP = Math.log(calculatePruningProb((Cell) tree.getRoot(), rootNodeType)[rootNodeType] * 1.0/numberOfTypes);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return logP;
    }

    /**
     * Calculate the conditional probability for current node using pruning algorithm
     * Assumes the node is not a leaf
     * @param node
     * @param nodeType, -1 if unspecified
     * @return
     */
    public double[] calculatePruningProb(Cell node, int nodeType) { // here nodeType refers to the type of the node at the beginning of the branch

        if(node.isLeaf()) { // if a node is a leaf here then it is the root. otherwise leaves are caught below in a previous recursive call of calculatePruningProb
            double[] result = new double[numberOfTypes];
            result[nodeType] = getProbaAtLeaves(node, -1)[nodeType];
            return result;
        }


        //TODO: potential for parallelization here
        int childIndex = 0;
        double[] pruningProbaFirstChild;

        if (node.getChild(childIndex).isLeaf())
            pruningProbaFirstChild = getProbaAtLeaves((Cell) node.getChild(childIndex), -1);
        else
            pruningProbaFirstChild = calculatePruningProb((Cell) node.getChild(childIndex), -1);

        childIndex = 1;
        double[] pruningProbaSecondChild;

        if (node.getChild(childIndex).isLeaf())
            pruningProbaSecondChild = getProbaAtLeaves((Cell) node.getChild(childIndex), -1);
        else
            pruningProbaSecondChild = calculatePruningProb( (Cell) node.getChild(childIndex), -1);


        double[] probsEndBranch = new double[numberOfTypes];
        for (int i = 0; i < numberOfTypes; i++) {
            for (int j = 0; j < numberOfTypes; j++) {
                for (int k = 0; k < numberOfTypes; k++) {
                    // TODO there is room for improvement here, at least for getProbabilityAtDivision, because calculations are repeated since symmetric on j and k.
                    // for now getProbabilityAtDivision calculations are very lightweight, only two types, so it's not an issue, but it may become one if really n types
                    probsEndBranch[i] += getProbabilityAtDivisionNode(i, j, k) * pruningProbaFirstChild[j] * pruningProbaSecondChild[k];
                }
            }
        }

        double[] probsStartBranch = new double[numberOfTypes];

        if(nodeType == -1) { // type of node is not fixed
            for (int i = 0; i < numberOfTypes; i++) {
                for (int j = 0; j < numberOfTypes; j++) {
                    probsStartBranch[i] += getProbabilityCellBranch(node, i, j) * probsEndBranch[j];
                }
            }
        }
        else if (nodeType > -1 && nodeType < numberOfTypes) { // type of this node is fixed
            for (int i = 0; i < numberOfTypes; i++) {
                if (i != nodeType)
                    probsStartBranch[i] = 0;
                else {
                    for (int j = 0; j < numberOfTypes; j++) {
                        probsStartBranch[i] += getProbabilityCellBranch(node, i, j) * probsEndBranch[j];
                    }
                }
            }
        }
        else {
            throw new IllegalStateException("Undefined nodeType. Value must be -1 if unspecified, and -1 < i < numberofTypes otherwise");
        }

        return probsStartBranch;
    }

    public double[] getProbaAtLeaves(Cell leaf, int leafType) { // here leafType refers to the type of the cell at the end of branch
        //TODO improve: for now, we're just starting with an equal proba for each type.
        double[] startingProbas  = new double[numberOfTypes];
        for (int i = 0; i < numberOfTypes; i++) {
            startingProbas[i] = 1.0/numberOfTypes;
        }

        double[] resultProbs = new double[numberOfTypes];


        if(leafType == -1) { // type of leaf is not fixed
            for (int i = 0; i < numberOfTypes; i++) {
                resultProbs[i] = 0;
                for (int j = 0; j < numberOfTypes; j++) {
                    resultProbs[i] += getProbabilityCellBranch(leaf, i, j) * startingProbas[j];
                }
            }
        }
        else if (leafType > -1 && leafType < numberOfTypes) { // type of this leaf is fixed
            for (int i = 0; i < numberOfTypes; i++) {
                resultProbs[i] = getProbabilityCellBranch(leaf, i, leafType) * startingProbas[leafType];
            }
        }
        else {
            throw new IllegalStateException("Undefined leafType. Value must be -1 if unspecified, and -1 < i < numberOfTypes otherwise");
        }

        return resultProbs;
    }

    public double getProbabilityCellBranch(Cell node, int typeStartBranch, int typeEndBranch) {

        double branchProb=0;

        if(node.getFate() == Cell.Fate.L) return lossProbInput.get().getValue();

        if(typeStartBranch == typeEndBranch) { // no type transition

            if (node.getFate() == Cell.Fate.D || node.getFate() == Cell.Fate.A) {

                int cellFate = node.getFate() == Cell.Fate.D ? 0 : 1; // get the index of the cellFate of the cell (0 for dividers, 1 for apoptosers)
                branchProb = fateProbabilitiesInput.get().get(typeEndBranch).getArrayValue(cellFate)
                        * Utils.getWeibullDensity(node.getEdgeLength(),
                        scaleWeibullInput.get().get(typeEndBranch).getArrayValue(cellFate),
                        shapeWeibullInput.get().get(typeEndBranch).getArrayValue(cellFate));

            }
            else if (node.getFate() == Cell.Fate.N) { //TODO potentially remove this fate (and have only unobserved instead)

                branchProb = fateProbabilitiesInput.get().get(typeEndBranch).getValue(2)
                        + fateProbabilitiesInput.get().get(typeEndBranch).getValue(0) // cell divides after end of branch
                        * Math.exp(-Math.pow(node.getEdgeLength() / scaleWeibullInput.get().get(typeEndBranch).getValue(0), shapeWeibullInput.get().get(typeEndBranch).getValue(0)))
                        + fateProbabilitiesInput.get().get(typeEndBranch).getValue(1) // cell dies after end of branch
                        * Math.exp(-Math.pow(node.getEdgeLength() / scaleWeibullInput.get().get(typeEndBranch).getValue(1), shapeWeibullInput.get().get(typeEndBranch).getValue(1)));

                if (transitionDuringLifetimeIsAllowed && fateProbabilitiesInput.get().get(typeEndBranch).getDimension() > 3) // check if this state can transition
                    //TODO if more than 1 fate that the cell can transition to, take it into account, either by summing the different probas or by having one big proba for all
                    branchProb += fateProbabilitiesInput.get().get(typeEndBranch).getValue(3) // cell transitions after end of branch
                            * Math.exp(-Math.pow(node.getEdgeLength() / scaleWeibullInput.get().get(typeEndBranch).getValue(2), shapeWeibullInput.get().get(typeEndBranch).getValue(2)));

            }
            else if (node.getFate() == Cell.Fate.U) {

                branchProb = fateProbabilitiesInput.get().get(typeEndBranch).getValue(0) // cell divides after end of branch
                        * Math.exp(-Math.pow(node.getEdgeLength() / scaleWeibullInput.get().get(typeEndBranch).getValue(0), shapeWeibullInput.get().get(typeEndBranch).getValue(0)))
                        + fateProbabilitiesInput.get().get(typeEndBranch).getValue(1) // cell dies after end of branch
                        * Math.exp(-Math.pow(node.getEdgeLength() / scaleWeibullInput.get().get(typeEndBranch).getValue(1), shapeWeibullInput.get().get(typeEndBranch).getValue(1)));

                if (transitionDuringLifetimeIsAllowed && fateProbabilitiesInput.get().get(typeEndBranch).getDimension() > 3) // check if this state can transition
                    branchProb += fateProbabilitiesInput.get().get(typeEndBranch).getValue(3) // cell transitions after end of branch
                            * Math.exp(-Math.pow(node.getEdgeLength() / scaleWeibullInput.get().get(typeEndBranch).getValue(2), shapeWeibullInput.get().get(typeEndBranch).getValue(2)));

            }
            else {
                throw new IllegalStateException("Invalid cell fate.");
            }

            // account for the experimental measures performed on the cells
            branchProb *= getProbabilityCellMeasures(node, typeStartBranch, typeEndBranch);

        }
        else { // typeStartBranch != typeEndBranch

            //TODO transition on branches are not properly implemented yet.

            //get index of transition in the matrixOfAllowedTransitions
            int indexTransition = typeStartBranch > typeEndBranch ? typeStartBranch * (numberOfTypes - 1) + typeEndBranch : typeStartBranch * (numberOfTypes - 1) + typeEndBranch - 1;

            if (!transitionDuringLifetimeIsAllowed || !matrixOfAllowedTransitionsInput.get().getValue(indexTransition))
                return 0; // configuration has probability zero (impossible type transition)
            else {
                if (node.getFate() == Cell.Fate.D || node.getFate() == Cell.Fate.A) {
                    int cellEndFate = node.getFate() == Cell.Fate.D ? 0 : 1; // get the index of the cellFate of the cell (0 for dividers, 1 for apoptosers)

                    // integration over transition point
                    // the "2" below refers to the fate transitioner. TODO refactor this "2", make cleaner which fate is which, when generalising
                    DensityProbOfStateTransition f = new DensityProbOfStateTransition(
                            node.getEdgeLength(),
                            shapeWeibullInput.get().get(typeStartBranch).getValue(2),
                            scaleWeibullInput.get().get(typeStartBranch).getValue(2),
                            shapeWeibullInput.get().get(typeEndBranch).getValue(cellEndFate),
                            scaleWeibullInput.get().get(typeEndBranch).getValue(cellEndFate));

                    try {
                        double integrationRes = numericalIntegrator.integrate(maxEval, f, 0, 1);
                        // the "3" below refers to the fate transitioner. TODO refactor this "3", make cleaner which fate is which, when generalising
                        branchProb = integrationRes * fateProbabilitiesInput.get().get(typeStartBranch).getValue(3) * fateProbabilitiesInput.get().get(typeEndBranch).getValue(cellEndFate);
                    }
                    catch (Exception e) { //TODO solve the pb with some integrations that don't happen: it's a pb with properly exploring the state space
                        branchProb = 0;
                    }

                }
                else if (node.getFate() == Cell.Fate.N) { //TODO potentially remove this fate (and have only unobserved instead)

                    DensityProbTransitionBeforeEndAndDivisionDeathAfter f = new DensityProbTransitionBeforeEndAndDivisionDeathAfter(
                            node.getEdgeLength(),
                            shapeWeibullInput.get().get(typeStartBranch).getValue(2),
                            scaleWeibullInput.get().get(typeStartBranch).getValue(2),
                            shapeWeibullInput.get().get(typeEndBranch).getValue(0),
                            scaleWeibullInput.get().get(typeEndBranch).getValue(0),
                            shapeWeibullInput.get().get(typeEndBranch).getValue(1),
                            scaleWeibullInput.get().get(typeEndBranch).getValue(1),
                            fateProbabilitiesInput.get().get(typeEndBranch).getValue(0),
                            fateProbabilitiesInput.get().get(typeEndBranch).getValue(1),
                            fateProbabilitiesInput.get().get(typeEndBranch).getValue(2));

                    try {
                        double integrationRes = numericalIntegrator.integrate(maxEval, f, 0, 1);
                        branchProb = fateProbabilitiesInput.get().get(typeStartBranch).getValue(3) * integrationRes; // cell transitions before end of branch, resulting cell does not do anything before end of branch

                    }
                    catch (Exception e) {
                        branchProb = 0; //TODO deal with that integration problem
                    }
                }
                else if (node.getFate() == Cell.Fate.U) {
                    //TODO implement.
                    // will be similar to fate N but removing the possibility that fate N appears.
                }
                else {
                    throw new IllegalStateException("Invalid cell fate.");
                }
            }


        }
        branchProb *= (1 - lossProbInput.get().getValue());
        return branchProb;
    }

    public double getProbabilityCellMeasures(Cell node, int typeStartBranch, int typeEndBranch) {
        double branchProb = 1;

        //TODO implement transitions on branches
        if(typeEndBranch != typeStartBranch) {
            return 0;
        }

        for (InputPair input : mapMeasureTypeToInput.keySet()) {
            MeasureType measureType = mapMeasureTypeToInput.get(input);

            //TODO remove (for debugging)
//            boolean isMeasureMeasured = node.hasMeasure(measureType);
//            boolean isInputHere = (input.getMean().get() != null);

            // if input is provided by user and the cell actually contains info about the measure
            if(input.getMean().get() != null && node.hasMeasure(measureType)) {
                double x = node.getSummaryValue(mapMeasureTypeToInput.get(input));
                if(Double.isNaN(x)) continue; // skip this measure if there is no summary value
                //TODO remove, for debugging
                if(Utils.getTruncatedNormalDensity(x, typeEndBranch, input, measureType) == 0) {
                    System.out.println("Problem here, branchProb is zero.");
                }

                //TODO remove, for debugging
                if( Utils.getTruncatedNormalDensity(x, typeEndBranch, input, measureType) < 1e-80)
                    System.out.println("Problem here, branchProb is going to be very low.");

                branchProb *= Utils.getTruncatedNormalDensity(x, typeEndBranch, input, measureType);
            }
        }
        return branchProb;
    }

    public double getProbabilityAtDivisionNode(int typeEndParentBranch, int typeStartChildBranch1, int typeStartChildBranch2) {

        if (!transitionUponDivisionIsAllowed) {
            if (typeEndParentBranch == typeStartChildBranch1 && typeStartChildBranch1 == typeStartChildBranch2)
                return 1;
            else return 0; // if one of the three states is different, local configuration is impossible without transition upon division
        }
        else {

            if(typeStartChildBranch1 > typeStartChildBranch2) {
                int temp = typeStartChildBranch2;
                typeStartChildBranch2 = typeStartChildBranch1;
                typeStartChildBranch1 = temp;
            }

            int indexInFlattenedArray = numberOfTypes * typeStartChildBranch1 + typeStartChildBranch2 - typeStartChildBranch1 * (typeStartChildBranch1+1)/2;

            return transitionUponDivisionProbsInput.get().get(typeEndParentBranch).getValue(indexInFlattenedArray);
        }
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        throw new RuntimeException("Not implemented.");
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    public String getTreeID() {
        if (rootIsHSCidxInput.get() != null)
            return rootIsHSCidxInput.get() + "";
        else return "";
    }

    public boolean isHSCtree() {
        return rootIsHSC;
    }

    public String getCorrespondanceLabelsToNodeArray() {

        String result = "";
        // for debugging
        int count = 0;

        for (Node node : tree.getNodesAsArray()) {
            int nodeNumber = node.getNr();
            if (nodeNumber != count) { // TODO remove ,for debugging
                System.out.println("Problem");
            }
            int newNodeNumber = Integer.parseInt(node.getID().substring(0, 1));
            count++;
            result += newNodeNumber;
        }
        return result;
    }

    private class DensityProbOfStateTransition implements UnivariateFunction {

        double cellLifetime;

        double shapeLifetimeOriginState;
        double shapeLifetimeEndState;

        double scaleLifetimeOriginState;
        double scaleLifetimeEndState;

        DensityProbOfStateTransition(double cellLifetime, double shapeWeibullHSC, double scaleWeibullHSC, double shapeWeibullMPP, double scaleWeibullMPP) {
            this.cellLifetime = cellLifetime;

            this.shapeLifetimeOriginState = shapeWeibullHSC;
            this.scaleLifetimeOriginState = scaleWeibullHSC;

            this.shapeLifetimeEndState = shapeWeibullMPP;
            this.scaleLifetimeEndState = scaleWeibullMPP;
        }

        /**
         * @param x the fraction of the branch length for which the cell-state transition happens
         * @return the probability density value of the cell transitioning at proportion 'x' of the branch length
         */
        @Override
        public double value(double x) {

            double timeUntilTransition = x * cellLifetime;
            double beforeTransition = Utils.getWeibullDensity(timeUntilTransition, scaleLifetimeOriginState, shapeLifetimeOriginState);
            double afterTransition = Utils.getConditionalWeibullDensity(cellLifetime, timeUntilTransition, scaleLifetimeEndState, shapeLifetimeEndState);

            return beforeTransition * afterTransition;
        }
    }

    private class DensityProbTransitionBeforeEndAndDivisionDeathAfter implements UnivariateFunction {

        double cellLifetime;

        double shapeLifetimeOriginState;
        double shapeLifetimeEndState1;
        double shapeLifetimeEndState2;

        double scaleLifetimeOriginState;
        double scaleLifetimeEndState1;
        double scaleLifetimeEndState2;

        double probaEndState1;
        double probaEndState2;
        double probaEndState3;


        DensityProbTransitionBeforeEndAndDivisionDeathAfter(double cellLifetime,
                                                            double shapeWeibullTrans, double scaleWeibullTrans,
                                                            double shapeWeibullDiv, double scaleWeibullDiv,
                                                            double shapeWeibullDeath, double scaleWeibullDeath,
                                                            double probaDiv, double probaDeath, double probaNothing) {
            this.cellLifetime = cellLifetime;

            this.shapeLifetimeOriginState = shapeWeibullTrans;
            this.scaleLifetimeOriginState = scaleWeibullTrans;

            this.shapeLifetimeEndState1 = shapeWeibullDiv;
            this.scaleLifetimeEndState1 = scaleWeibullDiv;

            this.shapeLifetimeEndState2 = shapeWeibullDeath;
            this.scaleLifetimeEndState2 = scaleWeibullDeath;

            this.probaEndState1 = probaDiv;
            this.probaEndState2 = probaDeath;
            this.probaEndState3 = probaNothing;
        }

        /**
         * @param x the fraction of the branch length for which the cell-state transition happens
         * @return the probability density value of the cell transitioning at proportion 'x' of the branch length
         */
        @Override
        public double value(double x) {

            double timeUntilTransition = x * cellLifetime;
            double beforeTransition = Utils.getWeibullDensity(timeUntilTransition, scaleLifetimeOriginState, shapeLifetimeOriginState);
            double afterTransition1 = probaEndState1 * Utils.getConditionalWeibullProbaAboveT(cellLifetime, timeUntilTransition, scaleLifetimeEndState1, shapeLifetimeEndState1);
            double afterTransition2 = probaEndState2 * Utils.getConditionalWeibullProbaAboveT(cellLifetime, timeUntilTransition, scaleLifetimeEndState2, shapeLifetimeEndState2);
            double afterTransition3 = probaEndState3; // HSC transitions and MPP does nothing

            return beforeTransition * (afterTransition1 + afterTransition2 + afterTransition3);
        }
    }

    public static void main(String[] args) throws IOException {
        String fileName = "../Data/Examples/MinusInfinityLik4.csv";
        LineageTree tree =  new LineageTree();
        tree.setInputValue("measuresCSVFile", fileName);

        tree.initAndValidate();

        System.out.println(tree.toString());

        LineageTreeProb probTree = new LineageTreeProb();
        probTree.setInputValue("tree", tree);

        List<RealParameter> fateProbabilities = new ArrayList<>();
        fateProbabilities.add(new RealParameter("1.0 0. 0."));
        fateProbabilities.add(new RealParameter("1.0 0. 0."));

        List<RealParameter> scaleWeibull = new ArrayList<>();
        scaleWeibull.add(new RealParameter("10 10"));
        scaleWeibull.add(new RealParameter("10 10"));

        List<RealParameter> shapeWeibull = new ArrayList<>();
        shapeWeibull.add(new RealParameter("1 1"));
        shapeWeibull.add(new RealParameter("1 1"));

        List<RealParameter> transitionUponDivisionProbs = new ArrayList<>();
        transitionUponDivisionProbs.add(new RealParameter("0.3 0.3 0.4"));
        transitionUponDivisionProbs.add(new RealParameter("0 0 1.0"));

        RealParameter meanNormalAreaGrowthRateInput = new RealParameter("20.0 20.0");
        RealParameter sdNormalAreaGrowthRateInput = new RealParameter("20.0 20.0");

        RealParameter meanNormalEccentricityInput = new RealParameter("0.5 .5");
        RealParameter sdNormalEccentricityInput = new RealParameter("0.3 0.3");

        RealParameter meanNormalInstantSpeedInput = new RealParameter("10. 10.");
        RealParameter sdNormalInstantSpeedInput = new RealParameter("10.0 10.0");

        RealParameter meanNormalTMRMProductionRateInput = new RealParameter("5.0 5.0");
        RealParameter sdNormalTMRMProductionRateInput = new RealParameter("2.0 2.0");

        RealParameter meanNormalTMRMMaxRateInput = new RealParameter("20.0 20.0");
        RealParameter sdNormalTMRMMaxRateInput = new RealParameter("40.0 40.0");

        RealParameter meanNormalROSProductionRateInput = new RealParameter("1.0 1.0");
        RealParameter sdNormalROSProductionRateInput = new RealParameter("2.0 2.0");

        RealParameter meanNormalCD71APCProductionRateInput = new RealParameter("5.0 5.0");
        RealParameter sdNormalCD71APCProductionRateInput = new RealParameter("2.0 2.0");

        RealParameter meanNormalCD71PEProductionRateInput = new RealParameter("5.0 5.0");
        RealParameter sdNormalCD71PEProductionRateInput = new RealParameter("10.0 10.0");

        RealParameter meanNormalcMycGFPMaxRateInput = new RealParameter("1.0 1.0");
        RealParameter sdNormalcMycGFPMaxRateInput = new RealParameter("10.0 10.0");

//        RealParameter meanNormalPerimeterGrowthRateInput = new RealParameter("0. 1.");
//        RealParameter sdNormalPerimeterGrowthRateInput = new RealParameter("5.0 10.0");

        probTree.setInputValue("transitionUponDivisionProbs", transitionUponDivisionProbs);
        probTree.setInputValue("fateProbabilities", fateProbabilities);
        probTree.setInputValue("scaleWeibull",scaleWeibull);
        probTree.setInputValue("shapeWeibull",shapeWeibull);
        probTree.setInputValue("lossProb", new RealParameter("0."));
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));

        probTree.setInputValue("meanAreaGrowthRate", meanNormalAreaGrowthRateInput);
        probTree.setInputValue("sdAreaGrowthRate", sdNormalAreaGrowthRateInput);

        probTree.setInputValue("meanEccentricity", meanNormalEccentricityInput);
        probTree.setInputValue("sdEccentricity", sdNormalEccentricityInput);

        probTree.setInputValue("meanInstantSpeed", meanNormalInstantSpeedInput);
        probTree.setInputValue("sdInstantSpeed", sdNormalInstantSpeedInput);

        probTree.setInputValue("meanTMRMProductionRate", meanNormalTMRMProductionRateInput);
        probTree.setInputValue("sdTMRMProductionRate", sdNormalTMRMProductionRateInput);


        probTree.setInputValue("meanTMRMMaxRate", meanNormalTMRMMaxRateInput);
        probTree.setInputValue("sdTMRMMaxRate", sdNormalTMRMMaxRateInput);


        probTree.setInputValue("meanROSProductionRate", meanNormalROSProductionRateInput);
        probTree.setInputValue("sdROSProductionRate", sdNormalROSProductionRateInput);

        probTree.setInputValue("meanCD71APCProductionRate", meanNormalCD71APCProductionRateInput);
        probTree.setInputValue("sdCD71APCProductionRate", sdNormalCD71APCProductionRateInput);

        probTree.setInputValue("meanCD71PEProductionRate", meanNormalCD71PEProductionRateInput);
        probTree.setInputValue("sdCD71PEProductionRate", sdNormalCD71PEProductionRateInput);

        probTree.setInputValue("meancMycGFPMaxRate", meanNormalcMycGFPMaxRateInput);
        probTree.setInputValue("sdcMycGFPMaxRate", sdNormalcMycGFPMaxRateInput);

//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);


        //probTree.setInputValue("meanPerimeterGrowthRate", meanNormalPerimeterGrowthRateInput);
//        probTree.setInputValue("sdPerimeterGrowthRate", sdNormalPerimeterGrowthRateInput);
//        probTree.setInputValue("meanCD41ProductionRate", meanNormalCD41ProductionRateInput);
//        probTree.setInputValue("sdCD41ProductionRate", sdNormalCD41ProductionRateInput);
//
//        probTree.setInputValue("meanFcgRIIIProductionRate", meanNormalFcgRIIIProductionRateInput);
//        probTree.setInputValue("sdFcgRIIIProductionRate", sdNormalFcgRIIIProductionRateInput);
        
        probTree.initAndValidate();
        double logP;
        logP = probTree.calculateLogP();

        System.out.println("logP = " + logP);

        System.out.println(Math.log(Double.MAX_VALUE));
    }


}



