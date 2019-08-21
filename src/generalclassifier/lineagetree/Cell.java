package generalclassifier.lineagetree;

import beast.evolution.tree.Node;
import generalclassifier.parametrization.ExperimentalMeasurements;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.SortedSet;

public class Cell extends Node {

    //TODO implement way to pass fate of cells as being A or L, and take the fates into account in calculations
    // so far only fates D and U are implemented

    public enum Fate {
        // D: divides, A: apoptoses, L: lost cell, U: unobserved
        D, A, U, L
    }

    Fate fate;
    protected final int trackNumber;
    boolean isIncompletelyMeasured;
    int childrenNum;
    NumberFormat numberFormat = new DecimalFormat("####0.000");

    HashMap<String, Double> experimentalMeasures = new HashMap<>();


    public Cell(int trackNumber) {
        this.trackNumber = trackNumber;
        this.childrenNum = 0;
        this.setFate(Fate.U);
        isIncompletelyMeasured = true;
    }

    public Fate getFate(){
        return this.fate;
    }

    public void setFate(Fate f) {
        switch (f) {
            case D:
                /*
                 * Set the cell as being a dividing cell.
                 * If not root, set isIncompletelyMeasured to false
                 * to mark that the cell has been measured over its entire lifetime
                 */
                if(!this.isRootCell()) this.setCompletelyMeasured();
                break;
            default:
        }
        this.fate = f;
    }

    public void setCompletelyMeasured(){
        this.isIncompletelyMeasured = false;
    }

    public boolean getIsCompletelyMeasured(){
        return !isIncompletelyMeasured;
    }

    public boolean isLostCell(){
        return fate == Fate.L;
    }

    public int getTrackNumber(){
        return this.trackNumber;
    }

    public int getCellGeneration() {
        return (int) Math.floor(Math.log(this.trackNumber)/Math.log(2)) + 1;
    }

    public boolean isRootCell(){
        return this.trackNumber == 1;
    }

    public static int findParent(int cellTrackNumber) {
        if(cellTrackNumber < 2)
            return 0; // no parent
        if(cellTrackNumber % 2 == 0)
            return cellTrackNumber/2;
        else
            return (cellTrackNumber -1)/2;
    }

    public void addExperimentalMeasures(List<ExperimentalMeasurements> measuresOnTree){
        for(ExperimentalMeasurements measure : measuresOnTree) {

            if(experimentalMeasures.containsKey(measure.getMeasurementTag())) continue;

            experimentalMeasures.put(measure.getMeasurementTag(), measure.getMeasuredValues().get(trackNumber));
        }
    }

    public double getValueMeasured(String measurementTag){
        if(experimentalMeasures.containsKey(measurementTag) && experimentalMeasures.get(measurementTag) != null)
            return experimentalMeasures.get(measurementTag);
        else
            return Double.NaN;
    }

    public boolean getIsIncompletelyMeasured(){
        return isIncompletelyMeasured;
    }

    /**
     * Renumber all the nodes in the tree following the conventions dictated by BEAST2.
     * Leaves must be numbered 0 to n. Internal nodes are n+1 to 2n-2. the root is numbered 2n -1.
     * Warning: this method is written with application to root cells only in mind.
     */
    public void labelNodesInTree() {

        if(!this.isRootCell())
            throw new IllegalArgumentException("labelNodesInTree method can only be applied to a root cell.");

        List<Cell> leaves = (List<Cell>) (List) this.getAllLeafNodes();
        int count = 0;
        for (Cell cell : leaves) {
            cell.setNr(count);
            count ++;
        }

        List<Cell> internalNodes = this.getAllInternalNodes();
        for (Cell cell : internalNodes) {
            cell.setNr(count);
            count ++;
        }
        // last, label the root
        this.setNr(count);
    }

    /**
     * get all internal nodes under this node, does not include the node it is applied to.
     */
    public List<Cell> getAllInternalNodes() {
        final List<Cell> internalNodes = new ArrayList<>();
        if(!this.isLeaf()) {
            for(Cell child : (List<Cell>) (List) this.getChildren())
                child.getAllInternalNodes(internalNodes);
        }
        return internalNodes;
    }

    // recursive
    public void getAllInternalNodes(final List<Cell> internalNodes) {
        if(!this.isLeaf()) {
            internalNodes.add(this);
        }

        for (Cell child : (List<Cell>) (List) getChildren())
            child.getAllInternalNodes(internalNodes);
    }

    public String toCSVRecord(SortedSet<String> sortedTags){
        String res = "";
        res += trackNumber + "," + fate.toString();
        for(String tag : sortedTags) {
            res += "," + numberFormat.format(experimentalMeasures.get(tag));
        }
        return res;
    }

}
