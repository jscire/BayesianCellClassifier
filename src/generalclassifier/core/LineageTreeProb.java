package generalclassifier.core;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

import generalclassifier.lineagetree.LineageTree;
import static generalclassifier.lineagetree.LineageTree.Fate.*;
import generalclassifier.utils.Utils;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

import java.io.*;
import java.util.*;


public class LineageTreeProb extends Distribution {

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

    //TODO set up as list to generalize
    public Input<RealParameter> transitionUponDivisionProbsInput = new Input<>("transitionProbs",
            "Vector containing transition-upon-division probabilities (q0, q1, q2). These correspond to transitions upon division."); // TODO rename input to be more explicit when this set-up works. For now I leave 'transitionProbs' to not mess up with what we already have that works.


    public Input<BooleanParameter> matrixOfAllowedTransitionsInput = new Input<> ("allowedTransitions",
            "Flattened matrix of allowed transitions between cell types. " +
                    "If mij is the boolean that characterizes transitions from type i to j, the element of index i*(n-1)+j if i>j, or i*(n-1) + j-1 otherwise, gives the index of element mij. n is the number of types.");


    LineageTree tree;

    boolean rootIsHSC;

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

        //TODO remove when generalized
        if(numberOfTypes != 2)
            throw new IllegalStateException("Number of types is not two. Such a configuration is not allowed by the implementation yet.");

        if(shapeWeibullInput.get().size() != numberOfTypes || scaleWeibullInput.get().size() != numberOfTypes)
            throw new IllegalStateException("Inputs fateProbabilities, shapeWeibull and scaleWeibull should have the same number of elements: the number of cell types.");

        if(matrixOfAllowedTransitionsInput.get() != null) {
            transitionDuringLifetimeIsAllowed = true;

            if(matrixOfAllowedTransitionsInput.get().getDimension() != numberOfTypes * (numberOfTypes -1))
                throw new IllegalStateException("Incorrect size of matrix of allowed transitions.");
        }
        else {
            transitionDuringLifetimeIsAllowed = false;
        }

        transitionUponDivisionIsAllowed = (transitionUponDivisionProbsInput.get() != null);

        if (lossProbInput.get() == null)
            throw new IllegalStateException("lossProb input is missing. Set to zero if cells are never lost.");

        //TODO add checks for parameter input format. In particular, check that fateProbabilities, shapeWeibull and scaleWeibull have correct dimension for each cell type.

    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        if (rootIsHSCidxInput.get() != null)
            rootIsHSC = rootIsHSCInput.get().getValue(rootIsHSCidxInput.get());
        else
            rootIsHSC = rootIsHSCInput.get().getValue();

        int rootNodeType = rootIsHSC? 0 : 1;

        //for now prior is that frequency is equal for all types TODO allow for more flexibility
        try {
            logP = Math.log(calculatePruningProb(tree.getRoot(), rootNodeType)[rootNodeType] * 1.0/numberOfTypes);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return logP;
    }


    public double getProbabilityCellBranch(Node node, int typeStartBranch, int typeEndBranch) {

        int nodeIdx = node.getNr();
        double branchProb=0;

        if(tree.getNodeFate(nodeIdx) == L) return lossProbInput.get().getValue();

        if(typeStartBranch == typeEndBranch) { // no type transition

            if (tree.getNodeFate(nodeIdx) == D || tree.getNodeFate(nodeIdx) == A) {

                int cellFate = tree.getNodeFate(nodeIdx) == D ? 0 : 1; // get the index of the cellFate of the cell (0 for dividers, 1 for apoptosers)
                branchProb = fateProbabilitiesInput.get().get(typeEndBranch).getArrayValue(cellFate)
                        * Utils.getWeibullDensity(tree.getEdgeLength(nodeIdx),
                        scaleWeibullInput.get().get(typeEndBranch).getArrayValue(cellFate),
                        shapeWeibullInput.get().get(typeEndBranch).getArrayValue(cellFate));

            }
            else if (tree.getNodeFate(nodeIdx) == N) { //TODO potentially remove this fate (and have only unobserved instead)

                branchProb = fateProbabilitiesInput.get().get(typeEndBranch).getValue(2)
                        + fateProbabilitiesInput.get().get(typeEndBranch).getValue(0) // cell divides after end of branch
                        * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx) / scaleWeibullInput.get().get(typeEndBranch).getValue(0), shapeWeibullInput.get().get(typeEndBranch).getValue(0)))
                        + fateProbabilitiesInput.get().get(typeEndBranch).getValue(1) // cell dies after end of branch
                        * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx) / scaleWeibullInput.get().get(typeEndBranch).getValue(1), shapeWeibullInput.get().get(typeEndBranch).getValue(1)));

                if (transitionDuringLifetimeIsAllowed && fateProbabilitiesInput.get().get(typeEndBranch).getDimension() > 3) // check if this state can transition
                    //TODO if more than 1 fate that the cell can transition to, take it into account, either by summing the different probas or by having one big proba for all
                    branchProb += fateProbabilitiesInput.get().get(typeEndBranch).getValue(3) // cell transitions after end of branch
                        * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx) / scaleWeibullInput.get().get(typeEndBranch).getValue(2), shapeWeibullInput.get().get(typeEndBranch).getValue(2)));

            }
            else if (tree.getNodeFate(nodeIdx) == U) {

                branchProb = fateProbabilitiesInput.get().get(typeEndBranch).getValue(0) // cell divides after end of branch
                        * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx) / scaleWeibullInput.get().get(typeEndBranch).getValue(0), shapeWeibullInput.get().get(typeEndBranch).getValue(0)))
                        + fateProbabilitiesInput.get().get(typeEndBranch).getValue(1) // cell dies after end of branch
                        * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx) / scaleWeibullInput.get().get(typeEndBranch).getValue(1), shapeWeibullInput.get().get(typeEndBranch).getValue(1)));

                if (transitionDuringLifetimeIsAllowed && fateProbabilitiesInput.get().get(typeEndBranch).getDimension() > 3) // check if this state can transition
                    branchProb += fateProbabilitiesInput.get().get(typeEndBranch).getValue(3) // cell transitions after end of branch
                            * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx) / scaleWeibullInput.get().get(typeEndBranch).getValue(2), shapeWeibullInput.get().get(typeEndBranch).getValue(2)));

            }
            else {
                throw new IllegalStateException("Invalid cell fate.");
            }

        }
        else { // typeStartBranch != typeEndBranch

            //get index of transition in the matrixOfAllowedTransitions
            int indexTransition = typeStartBranch > typeEndBranch ? typeStartBranch * (numberOfTypes - 1) + typeEndBranch : typeStartBranch * (numberOfTypes - 1) + typeEndBranch - 1;

            if (!transitionDuringLifetimeIsAllowed || !matrixOfAllowedTransitionsInput.get().getValue(indexTransition))
                return 0; // configuration has probability zero (impossible type transition)
            else {
                if (tree.getNodeFate(nodeIdx) == D || tree.getNodeFate(nodeIdx) == A) {
                    int cellEndFate = tree.getNodeFate(nodeIdx) == D ? 0 : 1; // get the index of the cellFate of the cell (0 for dividers, 1 for apoptosers)

                    // integration over transition point
                    // the "2" below refers to the fate transitioner. TODO refactor this "2", make cleaner which fate is which, when generalising
                    DensityProbOfStateTransition f = new DensityProbOfStateTransition(
                            tree.getEdgeLength(nodeIdx),
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
                else if (tree.getNodeFate(nodeIdx) == N) { //TODO potentially remove this fate (and have only unobserved instead)

                    DensityProbTransitionBeforeEndAndDivisionDeathAfter f = new DensityProbTransitionBeforeEndAndDivisionDeathAfter(
                            tree.getEdgeLength(nodeIdx),
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
                else if (tree.getNodeFate(nodeIdx) == U) {
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

    public double getProbabilityAtDivisionNode(Node node, int typeEndParentBranch, int typeStartChildBranch1, int typeStartChildBranch2) {

        if (!transitionUponDivisionIsAllowed) {
            if (typeEndParentBranch == typeStartChildBranch1 && typeStartChildBranch1 == typeStartChildBranch2)
                return 1;
            else return 0; // if one of the three states is different, local configuration is impossible without transition upon division

        }
        else {
            //TODO figure out a way to generalize transitions upon division to n states
            // for now only two types (0 and 1)

            if(typeEndParentBranch == 0) {
                switch (typeStartChildBranch1 + typeStartChildBranch2) {
                    case 0: // both daughters are of the type of the mother cell: 0 -> 0 + 0
                        return transitionUponDivisionProbsInput.get().getValue(2);

                    case 1: // asymmetric division 0 -> 0 + 1
                        return transitionUponDivisionProbsInput.get().getValue(1);

                    case 2: // 0 -> 1 + 1
                        return transitionUponDivisionProbsInput.get().getValue(0);

                    default:
                        throw new IllegalStateException("Invalid cell types. For now max n is 2");
                }
            }
            else if (typeStartChildBranch1 == 1 && typeStartChildBranch2 == 1) {
                return 1;
            }
            else { // TODO change when generalized
                return 0;
            }
        }
    }

    /**
     * Calculate the conditional probability for current node using pruning algorithm
     * Assumes the node is not a leaf
     * @param node
     * @param nodeType, -1 if unspecified
     * @return
     */
    public double[] calculatePruningProb(Node node, int nodeType) { // here nodeType referes to the type of the node at the beginning of the branch

        if(node.isLeaf()) { // if a node is a leaf here then it is the root. otherwise leaves are caught below in a previous recursive call of calculatePruningProb
            double[] result = new double[numberOfTypes];
            result[nodeType] = getProbaAtLeaves(node, -1)[nodeType];
            return result;
        }


        //TODO: potential for parallelization here
        int childIndex = 0;
        double[] pruningProbaFirstChild;

        if (node.getChild(childIndex).isLeaf())
            pruningProbaFirstChild = getProbaAtLeaves(node.getChild(childIndex), -1);
        else
            pruningProbaFirstChild = calculatePruningProb(node.getChild(childIndex), -1);

        childIndex = 1;
        double[] pruningProbaSecondChild;

        if (node.getChild(childIndex).isLeaf())
            pruningProbaSecondChild = getProbaAtLeaves(node.getChild(childIndex), -1);
        else
            pruningProbaSecondChild = calculatePruningProb(node.getChild(childIndex), -1);


        double[] probsEndBranch = new double[numberOfTypes];
        for (int i = 0; i < numberOfTypes; i++) {
            for (int j = 0; j < numberOfTypes; j++) {
                for (int k = 0; k < numberOfTypes; k++) {
                    // TODO there is room for improvement here, at least for getProbabilityAtDivision, because calculations are repeated since symmetric on j and k.
                    // for now getProbabilityAtDivision calculations are very lightweight, only two types, so it's not an issue, but it may become one if really n types
                    probsEndBranch[i] += getProbabilityAtDivisionNode(node, i, j, k) * pruningProbaFirstChild[j] * pruningProbaSecondChild[k];
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

    public double[] getProbaAtLeaves(Node leaf, int leafType) { // here leafType refers to the type of the cell at the end of branch
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

    }


}



