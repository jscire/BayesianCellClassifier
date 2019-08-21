package generalclassifier.core;

import beast.core.*;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import generalclassifier.lineagetree.Cell;
import generalclassifier.lineagetree.CellTree;
import generalclassifier.lineagetree.LineageTree;
import generalclassifier.parametrization.DistributionForMeasurement;
import generalclassifier.parametrization.ExperimentalMeasurements;
import generalclassifier.parametrization.Parametrization;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

import java.util.List;
import java.util.Random;


public class LineageTreeProb extends Distribution {

    public Input<CellTree> lineageTreeInput = new Input<>("tree",
            "Lineage tree.",
            Input.Validate.REQUIRED);

    public Input<Parametrization> parametrizationInput = new Input<>("parametrization",
            "",
            Input.Validate.REQUIRED);

    public Input<IntegerParameter> cellTypeInput = new Input<>("cellType",
            "Each number corresponds to the type of a cell. The index of a cell type is its tracknumber times the index of its tree." +
                    "The cell type defined here refers to the type at the beginning of the tree edge that carries the cell." +
                    "This is important if type transitions are to be allowed along tree edges." +
                    "It is not the case yet.");

    public Input<Integer> treeIdxInput = new Input<>("treeIdx",
            "If provided, this is an index into the rootIsHSC boolean " +
                    "parameter specifying which element corresponds to the root HSC status. Default: 0",
            0);

    public Input<Boolean> rootTypeOnlyInput = new Input<>("rootTypeOnly",
            "If true, only the types of the root cells are inferred," +
                    " the types of other cells are marginalized over. Default: true",
            false);


    CellTree tree;

    boolean sumOverDaughterCellTypes;

    int numberOfCellTypes;

    // TrapezoidIntegrator numericalIntegrator = new TrapezoidIntegrator(1e-6, 1e-6, 1, 60);
    IterativeLegendreGaussIntegrator numericalIntegrator = new IterativeLegendreGaussIntegrator(2, 1e-6, 1e-6);


    @Override
    public void initAndValidate() {

        tree = lineageTreeInput.get();

        numberOfCellTypes = parametrizationInput.get().numberOfCellTypes;

        //TODO add check that cellTypeInput contains only values between 0 and numberOfCellTypes -1

        //TODO implement check that treeIdxInput is long enough to get value

        sumOverDaughterCellTypes = rootTypeOnlyInput.get().booleanValue();
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        try {
            double p = 0;
            double[] resultPruning = calculatePruningProb((Cell) tree.getRoot());
            for (int i = 0; i < numberOfCellTypes; i++) {
                p += (resultPruning[i] * parametrizationInput.get().getTypeFreq(i));
            }
            logP = Math.log(p);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return logP;
    }

    public double[] calculatePruningProb(Cell node) { // here nodeType refers to the type of the node at the beginning of the branch

        double[] pruningProb = new double[numberOfCellTypes];

        int nodeType = getFixedCellType(node);

        if(node.isLeaf()) {
            if(nodeType == -1) { // leaf type is not fixed
                for (int i = 0; i < numberOfCellTypes; i++) {
                    pruningProb[i] = getCellProbability(node, i);
                }
            }
            else { // leaf type is fixed to nodeType
                pruningProb[nodeType] = getCellProbability(node, nodeType);
            }
            return pruningProb;
        }

        Cell child1 = (Cell) node.getChild(0);
        Cell child2 = (Cell) node.getChild(1);

        int child1Type = getFixedCellType(child1);
        int child2Type = getFixedCellType(child2);

        double[] pruningProbFirstChild;
        double[] pruningProbSecondChild;

        pruningProbFirstChild = calculatePruningProb(child1);
        pruningProbSecondChild = calculatePruningProb(child2);

        if(nodeType == -1) { // non-fixed type for parent cell
            if(child1Type == -1 && child2Type == -1) { // non-fixed type for either children
                for (int i = 0; i < numberOfCellTypes; i++) {
                    for (int j = 0; j < numberOfCellTypes; j++) {
                        for (int k = 0; k < numberOfCellTypes; k++) {
                            pruningProb[i] += parametrizationInput.get().getTransitionProbability(i, j, k) *
                                    pruningProbFirstChild[j] * pruningProbSecondChild[k];
                        }
                    }
                }
            }
            else if(child1Type > -1 && child2Type == -1) { // fixed type for just one of the two children
                for (int i = 0; i < numberOfCellTypes; i++) {
                    for (int k = 0; k < numberOfCellTypes; k++) {
                        pruningProb[i] += parametrizationInput.get().getTransitionProbability(i, child1Type, k) *
                                pruningProbFirstChild[child1Type] * pruningProbSecondChild[k];
                    }
                }
            }
            else if(child2Type > -1 && child1Type == -1) { // fixed type for just the other children
                for (int i = 0; i < numberOfCellTypes; i++) {
                    for (int j = 0; j < numberOfCellTypes; j++) {
                        pruningProb[i] += parametrizationInput.get().getTransitionProbability(i, j, child2Type) *
                                pruningProbFirstChild[j] * pruningProbSecondChild[child2Type];
                    }
                }
            }
            else if(child2Type > -1 && child1Type > -1) { // fixed type for both children
                for (int i = 0; i < numberOfCellTypes; i++) {
                        pruningProb[i] += parametrizationInput.get().getTransitionProbability(i, child1Type, child2Type) *
                                pruningProbFirstChild[child1Type] * pruningProbSecondChild[child2Type];
                }
            }
            else {
                throw new IllegalStateException("Invalid types for parentCell, child1 or child2.");
            }

            for (int i = 0; i < numberOfCellTypes; i++) {
                pruningProb[i] *= getCellProbability(node, i);
            }
        }
        else { // fixed type for parent cell
            if(child1Type == -1 && child2Type == -1) { // non-fixed type for either children
                for (int j = 0; j < numberOfCellTypes; j++) {
                    for (int k = 0; k < numberOfCellTypes; k++) {
                        pruningProb[nodeType] += parametrizationInput.get().getTransitionProbability(nodeType, j, k) *
                                pruningProbFirstChild[j] * pruningProbSecondChild[k];
                    }
                }
            }
            else if(child1Type > -1 && child2Type == -1) { // fixed type for just one of the two children
                for (int k = 0; k < numberOfCellTypes; k++) {
                    pruningProb[nodeType] += parametrizationInput.get().getTransitionProbability(nodeType, child1Type, k) *
                            pruningProbFirstChild[child1Type] * pruningProbSecondChild[k];
                }
            }
            else if(child2Type > -1 && child1Type == -1) { // fixed type for just the other children
                for (int j = 0; j < numberOfCellTypes; j++) {
                    pruningProb[nodeType] += parametrizationInput.get().getTransitionProbability(nodeType, j, child2Type) *
                            pruningProbFirstChild[j] * pruningProbSecondChild[child2Type];
                }
            }
            else if(child2Type > -1 && child1Type > -1) { // fixed type for both children
                    pruningProb[nodeType] += parametrizationInput.get().getTransitionProbability(nodeType, child1Type, child2Type) *
                            pruningProbFirstChild[child1Type] * pruningProbSecondChild[child2Type];
            }
            else {
                throw new IllegalStateException("Invalid types for child1 or child2.");
            }
            //TODO remove
            double d = getCellProbability(node, nodeType);
            pruningProb[nodeType] *= getCellProbability(node, nodeType);
        }

        return pruningProb;
    }

    public double getCellProbability(Cell cell, int cellType) {

        if(cell.isLostCell()) return parametrizationInput.get().getLossProbability();

        double p = 1.0;

        for (DistributionForMeasurement d : parametrizationInput.get().getDistributions()) {

            double measuredValue = cell.getValueMeasured(d.getMeasurementTag());

            p *= d.getProbability(measuredValue, cellType, cell.getIsIncompletelyMeasured(), cell.getFate(), cell.isRootCell());
        }

        p *= parametrizationInput.get().getFateProbability(cell.getFate(), cellType);

        return p;
    }

    /**
     * Return the fixed cell type at beginning of tree edge
     * returns -1 if type of cell is not fixed
     * !!! The type of the root must always be fixed !!!! (with the current code, but not necessary to keep it that way)
     * @param cell
     * @return
     */
    private int getFixedCellType(Cell cell){
        try {
            if(!sumOverDaughterCellTypes) {
                return cellTypeInput.get().getValue(treeIdxInput.get() + cell.getTrackNumber() - 1); // track numbers start at 1, so we substract 1
            }

            else if(cell.getTrackNumber() == 1) // cell is root cell
                return cellTypeInput.get().getValue(treeIdxInput.get());
            else return -1;
        }
        catch (Exception e) {
            throw new IllegalStateException("Attempting to get the type of a cell but failed.");
        }

    }

    public String getTreeID() {
        if (treeIdxInput.get() != null)
            return treeIdxInput.get() + "";
        else return "";
    }

    public String getCorrespondanceLabelsToNodeArray() {

        String result = "";
        for (Node node : tree.getNodesAsArray()) {
            int newNodeNumber = Integer.parseInt(node.getID().substring(0, 1));
            result += newNodeNumber;
        }
        return result;
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
        try {
            if (cellTypeInput.get().somethingIsDirty()) {
                if (sumOverDaughterCellTypes) {
                    if (cellTypeInput.get().isDirty(treeIdxInput.get()))
//                        System.out.println(treeIdxInput.get());
                        return true;
                    //TODO can we add "else return false"? we would need to know for sure that cellTypeInput is the only dirty element
                } else {
                    for (int cellLabel :  lineageTreeInput.get().getLabelsOfAllCellsInTree()) {
                        if (cellTypeInput.get().isDirty(treeIdxInput.get() + cellLabel - 1))
                            return true;
                    }
                }
            }

            for (BEASTInterface beastObject : listActiveBEASTObjects()) {
                if (beastObject == cellTypeInput.get()) continue;

                if (beastObject instanceof StateNode && ((StateNode) beastObject).somethingIsDirty()) {
                    return true;
                }

                if (beastObject instanceof CalculationNode && ((CalculationNode) beastObject).isDirtyCalculation()) {
                    return true;
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return false;

    }

    public static void main(String[] args) {

// first test commented out to test stuff further below, but it works as is
//        LineageTree tree =  new LineageTree();
//
//        ExperimentalMeasurements measuresSca1 = new ExperimentalMeasurements();
//        measuresSca1.initByName("measurementTag", "Sca1",
//                "values", "1:534, 2:543624.534432, 4:00.32, 0007: 012.32174839");
//        tree.setInputValue("measurement", measuresSca1);
//
//        ExperimentalMeasurements measuresLifetime = new ExperimentalMeasurements();
//        measuresLifetime.initByName("measurementTag", "lifetime",
//                "values", "1:4, 3:02, 5:54, 6:47");
//        tree.setInputValue("measurement", measuresLifetime);
//
//        ExperimentalMeasurements measuresTMRM = new ExperimentalMeasurements();
//        measuresTMRM.initByName("measurementTag", "TMRM",
//                "values", "3:543.789, 2:43.278, 1:243.43, 9:123");
//        tree.setInputValue("measurement", measuresTMRM);
//
//        tree.setInputValue("cellsInTree", "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16");
//        tree.setInputValue("lastGenerationOfInterest", 3);
//
//        tree.initAndValidate();
//
//        String tag1 = "lifetime";
//
//        DistributionForMeasurement distr1 = new DistributionForMeasurement();
//
//        distr1.initByName("measurementTag", tag1,
//                "parm1Distribution", new RealParameter("1.0 2.3"),
//                "parm2Distribution", new RealParameter("0.5 0.3"),
//                "zeroFraction", new RealParameter("0.01 0.01"),
//                "distributionType", "lognormal",
//                "estimateType", "max",
//                "isAppliedToRootCells", true);
//
//        distr1.initAndValidate();
//
//        String tag2 = "Sca1";
//
//        DistributionForMeasurement distr2 = new DistributionForMeasurement();
//
//        distr2.initByName("measurementTag", tag2,
//                "parm1Distribution", new RealParameter("1.0 2.3"),
//                "parm2Distribution", new RealParameter("0.5 0.3"),
//                "zeroFraction", new RealParameter("0.01 0.01"),
//                "distributionType", "lognormal",
//                "estimateType", "max",
//                "isAppliedToRootCells", true);
//
//        distr2.initAndValidate();
//
//        Parametrization parametrization = new Parametrization();
//
//        parametrization.initByName("distribution", distr1,
//                "distribution", distr2,
//                "transitionUponDivisionProbs", new RealParameter("0.2 0.25 0.55"),
//                "transitionUponDivisionProbs", new RealParameter("0.3 0.1 0.6"));
//
//        LineageTreeProb treeProb  =new LineageTreeProb();
//
//        treeProb.setInputValue("tree", tree);
//        treeProb.setInputValue("parametrization", parametrization);
//        treeProb.setInputValue("cellType", new IntegerParameter("0 1 0 0 0 1 1 1 1 1 1 1 1 1 1 1 0 1"));
//
//        treeProb.initAndValidate();
//
//        double logP = treeProb.calculateLogP();
//
//        System.out.println(logP);

        LineageTree tree =  new LineageTree();

        ExperimentalMeasurements measuresMot50 = new ExperimentalMeasurements();
        measuresMot50.initByName("measurementTag", "mot_50",
                "values", "1:0.5, 2:0.55, 3:0.45");
        tree.setInputValue("measurement", measuresMot50);

        tree.setInputValue("cellsInTree", "1,2,3,4,5,6,7");
        tree.setInputValue("lastGenerationOfInterest", 2);

        tree.initAndValidate();

        String tag1 = "lifetime";

        DistributionForMeasurement distr1 = new DistributionForMeasurement();

        distr1.initByName("measurementTag", "mot_50",
                "parm1Distribution", new RealParameter("3.68 131"),
                "parm2Distribution", new RealParameter("3.6 10.8"),
                "distributionType", "beta",
                "estimateType", "mean",
                "isAppliedToRootCells", true);

        distr1.initAndValidate();

        Parametrization parametrization = new Parametrization();

        parametrization.initByName("distribution", distr1,
                "transitionUponDivisionProbs", new RealParameter("0.2 0.25 0.55"),
                "transitionUponDivisionProbs", new RealParameter("0 0 1"));

        LineageTreeProb treeProb  =new LineageTreeProb();

        treeProb.setInputValue("tree", tree);
        treeProb.setInputValue("parametrization", parametrization);
        treeProb.setInputValue("cellType", new IntegerParameter("0"));
        treeProb.setInputValue("rootTypeOnly", "true");

        treeProb.initAndValidate();

        double logP = treeProb.calculateLogP();

        System.out.println(logP);

    }


}



