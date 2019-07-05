package generalclassifier.lineagetree;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import generalclassifier.parametrization.ExperimentalMeasurements;

import java.io.IOException;
import java.util.*;

public class LineageTree extends CellTree {

    public Input<List<ExperimentalMeasurements>> measuresOnCellsInput = new Input<>("measurement",
            "List of experimental measurements made on cells in the tree.",
            new ArrayList<>());

    public Input<String> cellsInTreeInput = new Input<>("cellsInTree",
            "Comma separated list of labels of cells in tree. " +
            "Cells must be numbered following the convention: root cell is '1'," +
                    " the daughters of a cell 'n' are '2n' and '2n+1'. ");

    //TODO ongoing work, continue
    public Input<Boolean> isLastgenerationCompleteInput = new Input<>("isLastGenerationComplete",
            "Have the cells in the last generation of interest been tracked for their entire lifespan?" +
                    "", Boolean.TRUE );

    SortedSet<Integer> labelsOfAllCellsInTree;


    @Override
    public void initAndValidate() {

        String tag;
        uniqueMeasurementTags = new TreeSet<>();

        //check for duplication in measurement tags
        for(ExperimentalMeasurements measure : measuresOnCellsInput.get()) {

            tag = measure.getMeasurementTag();

            if(uniqueMeasurementTags.contains(tag))
                throw new IllegalArgumentException("Duplicated tag: " + tag + " in list of measurements.");
            else
                uniqueMeasurementTags.add(tag);
        }

        if(cellsInTreeInput.get() != null) {
            //TODO add check for correct formatting of the cellsInTree string (all cell numbers must be > 0 integers)
            String[] labelsOfCellsInTreeAsStrings = cellsInTreeInput.get().replaceAll("\\s","").split(",");

            labelsOfAllCellsInTree = new TreeSet<>();

            for(String s : labelsOfCellsInTreeAsStrings) {
                labelsOfAllCellsInTree.add(Integer.parseInt(s));
            }
        }
        else {
            labelsOfAllCellsInTree = new TreeSet<>();
            for (int i = 1; i < Math.pow(2, lastGenerationOfInterestInput.get()); i++) {
                labelsOfAllCellsInTree.add(i);
            }
        }

        Map<Integer, Cell> cellsOfInterest =  buildAllCellsOfInterest();

        labelsOfCellsOfInterest = new TreeSet<Integer>();

        for(Integer i : cellsOfInterest.keySet()) {
            labelsOfCellsOfInterest.add(i);
        }

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
            rootCell.setFate(Cell.Fate.D);

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

    public String toCSV() {
        return super.toString();
    }

    public SortedSet<Integer> getLabelsOfCellsOfInterest(){
        return labelsOfCellsOfInterest;
    }


    //LineageTrees are not treated as StateNodes, so we override the store/restore with empty methods.
    @Override
    protected void store() {}
    @Override
    public void restore() {}

    public static void main(String[] parms) {

        LineageTree tree =  new LineageTree();

        ExperimentalMeasurements measuresSca1 = new ExperimentalMeasurements();
        measuresSca1.initByName("measurementTag", "Sca1",
                "values", "1:534, 2:543624.534432, 4:00.32, 0007: 012.32174839");
        tree.setInputValue("measurement", measuresSca1);

        ExperimentalMeasurements measuresLifetime = new ExperimentalMeasurements();
        measuresLifetime.initByName("measurementTag", "lifetime",
                "values", "1:4, 3:02, 5:54, 6:47");
        tree.setInputValue("measurement", measuresLifetime);

        ExperimentalMeasurements measuresTMRM = new ExperimentalMeasurements();
        measuresTMRM.initByName("measurementTag", "TMRM",
                "values", "3:543.789, 2:43.278, 1:243.43, 9:123");
        tree.setInputValue("measurement", measuresTMRM);

        tree.setInputValue("cellsInTree", "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16");
        tree.setInputValue("lastGenerationOfInterest", 3);

        tree.initAndValidate();

        System.out.println(tree.toString());
    }

}
