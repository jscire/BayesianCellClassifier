package generalclassifier.lineagetree;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

public class LineageTree extends Tree {

    public Input<String> newickInput = new Input<>("newick",
            "Newick string describing lineage tree.",
            Input.Validate.REQUIRED);


    public enum Fate {
        // D: divides, A: apoptoses, L: lost cell, N: does nothing (not sure we'll keep that), U: unobserved fate
        D, A, L, N, U
    }

    private Fate[] nodeFates;

    private double edgeLengths[];

    public boolean allowTransitionOnEdge = false;

    //TODO change that before going for bigger trees
    static final int maxTreeSize = 7; // maximal number of nodes a LineageTree can have

    @Override
    public void initAndValidate() {

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

}
