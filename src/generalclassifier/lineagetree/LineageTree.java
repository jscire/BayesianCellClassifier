package generalclassifier.lineagetree;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class LineageTree extends Tree {

    public Input<String> experimentalMeasuresFileInput = new Input<>("measures",
            "CSV file containing all the measures performed on this tree.");

    public Input<String> frameRateInput = new Input<>("frameRate",
            "Specify only if time is not explicitly measured in input csv file." +
                    "Accepted values:  'day', 'hour', 'minute', 'second'");

    static final int rootCellNumberInInput = 1;

    private Cell[] cells;


    @Override
    public void initAndValidate() {

        if(experimentalMeasuresFileInput.get() == null) throw new IllegalArgumentException("No input file to initialize the tree.");

        Map<Integer, Cell> cells  = new HashMap<>();
        try {
            LineageTreeParser parser = new LineageTreeParser(experimentalMeasuresFileInput.get());
            cells = parser.parseRawCells();
        } catch (IOException e) {
            e.printStackTrace();
        }

        Cell rootCell = Cell.buildTreeAndGetRoot(cells, rootCellNumberInInput);
        rootCell.labelNodesInTree();


        setRoot(rootCell);
        initArrays();

        // apply an offset to the nodeheights so that they conform to the standard format
        double maxHeight = 0;
        for(Node node : getNodesAsArray()) {
            if(node.getHeight() > maxHeight) maxHeight = node.getHeight();
        }

        for(Node node : getNodesAsArray()) {
            node.setHeight(maxHeight - node.getHeight());
        }

        super.initAndValidate();

        //TODO add fates
    }

    // naive way of getting a node by its label number
    public Node getNodeByLabel(int nodeLabel){
        for (Node node : getNodesAsArray()) {
            if (node.getNr() == nodeLabel) return node;
        }
        return null;
    }


    //TODO that's a dirty attempt to go around the fact that nodes in Tree should be numbered the standard way for BEAST, not the way they are in the input data
    //TODO does it even make sense to do things that way?
    /**
     * StateNode implementation *
     */
    @Override
    protected void store() {}
    @Override
    public void restore() {}

    public static void main(String[] parms) {

        String fileName = "../Data/Examples/toyFile.csv";
        LineageTree tree =  new LineageTree();
        tree.setInputValue("measures", fileName);

        tree.initAndValidate();

        System.out.println(tree.toString());
    }

}
