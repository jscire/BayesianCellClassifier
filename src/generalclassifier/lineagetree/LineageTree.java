package generalclassifier.lineagetree;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import generalclassifier.parametrization.ExperimentalMeasurements;

import java.io.IOException;
import java.util.*;

public class LineageTree extends Tree {

    public Input<List<ExperimentalMeasurements>> measuresOnCellsInput = new Input<>("measurement",
            "List of experimental measurements made on cells in the tree.");

    public Input<String> cellsInTreeInput = new Input<>("cellsInTree",
            "Comma separated list of labels of cells in tree. " +
            "Cells must be numbered following the convention: root cell is '1'," +
                    " the daughters of a cell 'n' are '2n' and '2n+1'. ", Input.Validate.REQUIRED);

    //TODO ask about convention for numbering generations
    public Input<Integer> lastGenerationOfInterestInput = new Input<>("lastGenerationOfInterest",
            "Root cell is considered to be at generation 1." +
                    "If not specified, all cells in 'cellsInTree' will be used in the analysis." +
                    "The point of having more cells in cellsInTree than what the lastGenerationOfInterest imposes" +
                    " is to provide information what cells do at the end of their lifetime: do they divide or not?.",
            Input.Validate.REQUIRED);

    SortedSet<Integer> labelsOfAllCellsInTree;

    HashSet<String> uniqueMeasurementTags;

    int numberOfCellsInTree;

    @Override
    public void initAndValidate() {


        String tag;

        //check for duplication in measurement tags
        for(ExperimentalMeasurements measure : measuresOnCellsInput.get()) {

            tag = measure.getMeasurementTag();

            if(uniqueMeasurementTags.contains(tag))
                throw new IllegalArgumentException("Duplicated tag: " + tag + " in list of measurements.");
            else
                uniqueMeasurementTags.add(tag);
        }

        //TODO add check fir correct formatting og the cellsInTree string
        String[] labelsOfCellsInTreeAsStrings = cellsInTreeInput.get().replaceAll("\\s","").split(",");

        labelsOfAllCellsInTree = new TreeSet<Integer>();

        for(String s : labelsOfCellsInTreeAsStrings) {
            labelsOfAllCellsInTree.add(Integer.parseInt(s));
        }

        Map<Integer, Cell> cellsOfInterest =  buildAllCellsOfInterest();

        Cell rootCell = buildTreeAndGetRoot(1, cellsOfInterest);
        rootCell.labelNodesInTree();

        setRoot(rootCell);
        initArrays();

        super.initAndValidate();
    }

    // naive way of getting a node by its label number
    public Node getNodeByLabel(int nodeLabel){
        for (Node node : getNodesAsArray()) {
            if (node.getNr() == nodeLabel) return node;
        }
        return null;
    }

    public Map<Integer, Cell> buildAllCellsOfInterest() {

        int maxCellNumber = (int) Math.pow(2, lastGenerationOfInterestInput.get()) -1;

        Map<Integer, Cell> cells  = new HashMap<>();

        for(Integer cellLabel : labelsOfAllCellsInTree) {
            if (cellLabel > maxCellNumber || cells.containsKey(cellLabel))
                continue;

            Cell newCell = new Cell(cellLabel);
            newCell.addExperimentalMeasures(measuresOnCellsInput.get());

            cells.put(cellLabel, newCell);
        }

        return cells;
    }

    public Cell buildTreeAndGetRoot(int rootKey, Map<Integer, Cell> cellsOfInterest) {
        if(!cellsOfInterest.containsKey(rootKey))
            throw new IllegalArgumentException("There is no cell with label " + rootKey + " in the set of cells of interest.");

        Cell rootCell = cellsOfInterest.get(rootKey);

        if(labelsOfAllCellsInTree.contains(2*rootKey) || labelsOfAllCellsInTree.contains(2*rootKey + 1))
            rootCell.setAsDivider();

        if(cellsOfInterest.containsKey(2*rootKey)) {
            Cell child1 = buildTreeAndGetRoot(2*rootKey, cellsOfInterest);
            child1.setParent(rootCell);
            rootCell.addChild(child1);
        }
        if(cellsOfInterest.containsKey(2*rootKey + 1)) {
            Cell child2 = buildTreeAndGetRoot(2*rootKey + 1, cellsOfInterest);
            child2.setParent(rootCell);
            rootCell.addChild(child2);
        }

        return rootCell;
    }

    //TODO implement proper toString method to at least print the cell tracknumbers in the metadata.
    @Override
    public String toString() {
        return super.toString();
    }


    //LineageTrees are not treated as StateNodes, so we override the store/restore with empty methods.
    @Override
    protected void store() {}
    @Override
    public void restore() {}

    public static void main(String[] parms) {
        //TODO update example to adapt to new code.

//        String fileName = "../Data/Examples/toyFile.csv";
        String fileName = "../Data/Examples/testFile_shortLife.csv";
        LineageTree tree =  new LineageTree();
        tree.setInputValue("measuresCSVFile", fileName);

        tree.initAndValidate();

        System.out.println(tree.toString());
    }

}
