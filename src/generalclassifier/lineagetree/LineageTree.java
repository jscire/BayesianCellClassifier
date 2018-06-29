package generalclassifier.lineagetree;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class LineageTree extends Tree {

    public Input<String> newickInput = new Input<>("newick",
            "Newick string describing lineage tree.");

    public Input<String> experimentalMeasuresFileInput = new Input<>("measures",
            "CSV file containing all the measures performed on this tree.");

    public Input<String> frameRateInput = new Input<>("frameRate",
            "Specify only if time is not explicitly measured in input csv file." +
                    "Accepted values:  'day', 'hour', 'minute', 'second'");


    public enum Fate {
        // D: divides, A: apoptoses, L: lost cell, N: does nothing (not sure we'll keep that), U: unobserved
        D, A, L, N, U
    }

    static final int rootCellNumberInInput = 1;

    private Cell[] cells;

    private Fate[] nodeFates;

    private double edgeLengths[];

    //TODO change that
    static final int maxTreeSize = 50; // maximal number of nodes a LineageTree can have



    @Override
    public void initAndValidate() {
        // keep the previous method with newick string inputs
        if(newickInput.get() != null) {
            initWithNewickString();
            return;
        }

        if(experimentalMeasuresFileInput.get() == null) throw new IllegalArgumentException("No input file nor newick string to initialize the tree.");

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

        //TODO add root branch
        //TODO add fates
    }

    public void initWithNewickString() {
        String newickString = newickInput.get();

        assignFromWithoutID(new TreeParser(newickString, false, false, true, 0));

        //TODO remove if below works
//        nodeFates = new Fate[getNodeCount()]; // original version
//        edgeLengths = new double[getNodeCount()];

        nodeFates = new Fate[maxTreeSize];
        edgeLengths = new double[maxTreeSize];

        for (Node node : getNodesAsArray()) {

            // TODO check that ok to always do this renumbering (helpful when looking at transitions on branches)
            int newNodeNumber = Integer.parseInt(node.getID().substring(0, 1)) -1; // the "-1" is there to account for the fact that cell indices start at 1 in the datasets
            node.setNr(newNodeNumber);

            String fateStr = node.getID().substring(1);
            Fate fate;
            switch(fateStr) {
                case "D":
                    fate = Fate.D;
                    break;
                case "A":
                    fate = Fate.A;
                    break;
                case "L":
                    fate = Fate.L;
                    break;
                case "N":
                    fate = Fate.N;
                    break;
                case "U":
                    throw new IllegalStateException("Unobserved fate 'U' is not fully implemented yet.");
                    //fate = Fate.U;
                    //break;
                default:
                    throw new IllegalArgumentException("Unknown cell fate '" + fateStr + "'");
            }

            nodeFates[node.getNr()] = fate;

            if (!node.isRoot())
                edgeLengths[node.getNr()] = node.getParent().getHeight()-node.getHeight();
        }

        String rootEdgeLengthStr = newickString.substring(1+newickString.lastIndexOf(":"));
        if (rootEdgeLengthStr.endsWith(";"))
            rootEdgeLengthStr = rootEdgeLengthStr.substring(0, rootEdgeLengthStr.length()-1);

        edgeLengths[getRoot().getNr()] = Double.parseDouble(rootEdgeLengthStr);
    }

    public Fate getNodeFate(int nodeNr) {
        return nodeFates[nodeNr];
    }

    public double getEdgeLength(int nodeNr) {
        return edgeLengths[nodeNr];
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
