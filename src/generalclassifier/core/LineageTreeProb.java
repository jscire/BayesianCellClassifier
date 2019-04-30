package generalclassifier.core;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import generalclassifier.lineagetree.Cell;
import generalclassifier.lineagetree.LineageTree;
import generalclassifier.parametrization.DistributionForMeasurement;
import generalclassifier.parametrization.Parametrization;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

import java.io.IOException;
import java.util.*;


public class LineageTreeProb extends Distribution {

    public Input<LineageTree> lineageTreeInput = new Input<>("tree",
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
                    "parameter specifying which element corresponds to the root " +
                    "HSC status.");

    public Input<Boolean> rootTypeOnlyInput = new Input<>("rootTypeOnly", "If true, only the types of the root cells are inferred, the types of other cells are marginalized over. Default: true", true);


    LineageTree tree;

    boolean sumOverDaughterCellTypes;

    int numberOfCellTypes;

    // TrapezoidIntegrator numericalIntegrator = new TrapezoidIntegrator(1e-6, 1e-6, 1, 60);
    IterativeLegendreGaussIntegrator numericalIntegrator = new IterativeLegendreGaussIntegrator(2, 1e-6, 1e-6);


    @Override
    public void initAndValidate() {

        tree = lineageTreeInput.get();

        numberOfCellTypes = parametrizationInput.get().numberOfCellTypes;

        //TODO remove when tested with more than 2 types.
        if(numberOfCellTypes != 2)
            throw new IllegalStateException("Number of types is not two. Such a configuration is not tested yet.");


        //TODO add check that cellTypeInput contains only values between 0 and numberOfCellTypes -1

        //TODO implement check that treeIdxInput is long enough to get value

        sumOverDaughterCellTypes = rootTypeOnlyInput.get().booleanValue();
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        //for now prior is that frequency is equal for all types
        //TODO allow for more flexibility
        try {
            double p = 0;
            double[] resultPruning = calculatePruningProb((Cell) tree.getRoot());
            for (int i = 0; i < numberOfCellTypes; i++) {
                p += (resultPruning[i] * 1.0/ numberOfCellTypes);
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
            else {
                throw new IllegalStateException("Invalid types for parentCell, child1 or child2.");
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
            else {
                throw new IllegalStateException("Invalid types for child1 or child2.");
            }
        }

        return pruningProb;
    }

    public double getCellProbability(Cell cell, int cellType) {

        if(cell.isLostCell()) return parametrizationInput.get().getLossProbability();

        double p = 1.0;

        for (DistributionForMeasurement d : parametrizationInput.get().getDistributions()) {

            if(cell.isRootCell() && !d.isAppliedOnRootCells()) continue; // cell is root cell and measure is not taken into account for root cells

            double measuredValue = cell.getValueMeasured(d.getMeasurementTag());

            if(Double.isNaN(measuredValue)) continue; // no data on this cell for this particular measure

            p *= d.getProbability(measuredValue, cellType, cell.getIsIncompletelyMeasured(), cell.getFate());
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
            if(!sumOverDaughterCellTypes)
                return cellTypeInput.get().getValue(treeIdxInput.get() * cell.getTrackNumber());
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
        return true;
    } //TODO implement real check, do not recalculate if nothing changed.

    public static void main(String[] args) throws IOException {
        //TODO update example
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
        fateProbabilities.add(new RealParameter("1.0 0. 0."));

        List<RealParameter> scaleWeibull = new ArrayList<>();
        scaleWeibull.add(new RealParameter("10 10"));
        scaleWeibull.add(new RealParameter("10 10"));
        scaleWeibull.add(new RealParameter("10 10"));

        List<RealParameter> shapeWeibull = new ArrayList<>();
        shapeWeibull.add(new RealParameter("1 1"));
        shapeWeibull.add(new RealParameter("1 1"));
        shapeWeibull.add(new RealParameter("1 1"));

        List<RealParameter> transitionUponDivisionProbs = new ArrayList<>();
        transitionUponDivisionProbs.add(new RealParameter("0.3 0.1 0.1 0.15 0.15 0.2"));
        transitionUponDivisionProbs.add(new RealParameter("0 0 0 1.0 0 0 "));
        transitionUponDivisionProbs.add(new RealParameter("0 0 0 0 0 1.0 "));


        RealParameter meanLogNormalAreaGrowthRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter sdLogNormalAreaGrowthRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter zeroFractionAreaGrowthRateInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanNormalEccentricityInput = new RealParameter("0.5 .5 .5");
        RealParameter sdNormalEccentricityInput = new RealParameter("0.3 0.3 .3");

        RealParameter meanLogNormalInstantSpeedInput = new RealParameter("1. 1. 1.");
        RealParameter sdLogNormalInstantSpeedInput = new RealParameter("1.0 1.0 1.");
        RealParameter zeroFractionInstantSpeedInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanLogNormalTMRMProductionRateInput = new RealParameter("3.0 3.0 3.0");
        RealParameter sdLogNormalTMRMProductionRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter zeroFractionTMRMProductionRateInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanLogNormalTMRMMaxRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter sdLogNormalTMRMMaxRateInput = new RealParameter("1.0 1.0 1.0");
        RealParameter zeroFractionTMRMMaxRateInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanLogNormalROSProductionRateInput = new RealParameter("1.0 1.0 1.0");
        RealParameter sdLogNormalROSProductionRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter zeroFractionROSProductionRateInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanLogNormalCD71APCProductionRateInput = new RealParameter("5.0 5.0 5.0");
        RealParameter sdLogNormalCD71APCProductionRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter zeroFractionCD71APCProductionRateInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanLogNormalCD71PEProductionRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter sdLogNormalCD71PEProductionRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter zeroFractionCD71PEProductionRateInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanLogNormalcMycGFPMaxRateInput = new RealParameter("1.0 1.0 1.0");
        RealParameter sdLogNormalcMycGFPMaxRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter zeroFractioncMycGFPMaxRateInput = new RealParameter("0.1 0.1 0.1");


//        RealParameter meanNormalPerimeterGrowthRateInput = new RealParameter("0. 1.");
//        RealParameter sdNormalPerimeterGrowthRateInput = new RealParameter("5.0 10.0");

        probTree.setInputValue("transitionUponDivisionProbs", transitionUponDivisionProbs);
        probTree.setInputValue("fateProbabilities", fateProbabilities);
        probTree.setInputValue("scaleWeibull",scaleWeibull);
        probTree.setInputValue("shapeWeibull",shapeWeibull);
        probTree.setInputValue("lossProb", new RealParameter("0."));

//        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.setInputValue("rootType", new IntegerParameter("2"));

        probTree.setInputValue("meanAreaGrowthRate", meanLogNormalAreaGrowthRateInput);
        probTree.setInputValue("sdAreaGrowthRate", sdLogNormalAreaGrowthRateInput);
        probTree.setInputValue("zeroFractionAreaGrowthRate", zeroFractionAreaGrowthRateInput);

        probTree.setInputValue("meanEccentricity", meanNormalEccentricityInput);
        probTree.setInputValue("sdEccentricity", sdNormalEccentricityInput);

        probTree.setInputValue("meanInstantSpeed", meanLogNormalInstantSpeedInput);
        probTree.setInputValue("sdInstantSpeed", sdLogNormalInstantSpeedInput);
        probTree.setInputValue("zeroFractionInstantSpeed", zeroFractionInstantSpeedInput);

        probTree.setInputValue("meanTMRMProductionRate", meanLogNormalTMRMProductionRateInput);
        probTree.setInputValue("sdTMRMProductionRate", sdLogNormalTMRMProductionRateInput);
        probTree.setInputValue("zeroFractionTMRMProductionRate", zeroFractionROSProductionRateInput);

        probTree.setInputValue("meanTMRMMaxRate", meanLogNormalTMRMMaxRateInput);
        probTree.setInputValue("sdTMRMMaxRate", sdLogNormalTMRMMaxRateInput);
        probTree.setInputValue("zeroFractionTMRMMaxRate", zeroFractionTMRMMaxRateInput);

        probTree.setInputValue("meanROSProductionRate", meanLogNormalROSProductionRateInput);
        probTree.setInputValue("sdROSProductionRate", sdLogNormalROSProductionRateInput);
        probTree.setInputValue("zeroFractionROSProductionRate", zeroFractionROSProductionRateInput);

        probTree.setInputValue("meanCD71APCProductionRate", meanLogNormalCD71APCProductionRateInput);
        probTree.setInputValue("sdCD71APCProductionRate", sdLogNormalCD71APCProductionRateInput);
        probTree.setInputValue("zeroFractionCD71APCProductionRate", zeroFractionCD71APCProductionRateInput);

        probTree.setInputValue("meanCD71PEProductionRate", meanLogNormalCD71PEProductionRateInput);
        probTree.setInputValue("sdCD71PEProductionRate", sdLogNormalCD71PEProductionRateInput);
        probTree.setInputValue("zeroFractionCD71PEProductionRate", zeroFractionCD71PEProductionRateInput);

        probTree.setInputValue("meancMycGFPMaxRate", meanLogNormalcMycGFPMaxRateInput);
        probTree.setInputValue("sdcMycGFPMaxRate", sdLogNormalcMycGFPMaxRateInput);
        probTree.setInputValue("zeroFractioncMycGFPMaxRate", zeroFractioncMycGFPMaxRateInput);

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



