package generalclassifier.lineagetree;

import beast.evolution.tree.Node;
import org.apache.commons.csv.CSVRecord;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

//TODO remove everything that is not used
public class Cell extends Node {

    public enum Fate {
        // D: divides, A: apoptoses, L: lost cell, N: does nothing (not sure we'll keep that), U: unobserved
        D, A, L, N, U
    }

    Fate fate;
    int trackNumber;
    double edgeLength;

   // List<double[]> timePoints = new ArrayList<> ();

    List<MeasureType> measureTypes = new ArrayList<>();

    protected List<ExperimentalMeasure> measures = new ArrayList<> ();

    List<Double> timePoints = new ArrayList<>();

    public Cell(int trackNumber) {
        this.trackNumber = trackNumber;
    }

    public Cell(int trackNumber, List<MeasureType> measureTypes) {
        this.trackNumber = trackNumber;
        this.measureTypes = measureTypes;
        for (MeasureType measureType : measureTypes) {
            this.measures.add(new ExperimentalMeasure(measureType));
        }
    }

    public double getEdgeLength(){
        return this.edgeLength;
    }

    public void setEdgeLength() {
        //TODO remove
//        this.edgeLength = timePoints.get(timePoints.size()-1) -  timePoints.get(0);
        //TODO maybe replace +1 with +0.5?
        this.edgeLength = timePoints.get(timePoints.size()-1) -  timePoints.get(0) + 1; // add 1 to account for the fact that a cell can only last 1 timeInterval (only seen once, but doesnot have a lifetime of zero)
    }

    public Fate getFate(){
        return this.fate;
    }

    //TODO not enough to allow for more fates than just dividers (need to read actually it from the csv file)
    public void setFate(Fate f) {
        this.fate = f;
    }

    public boolean hasMeasure(MeasureType measureType) {
        return measureTypes.contains(measureType); // todo check that actually returns true when measure is in the list
    }

    /**
     * Only apply to a measure that was measured for this cell
     * and a measure that has been summarized (ExperimentalMeasure.summarize() was called)
     * @param measureType
     * @return summary value
     */
    public double getSummaryValue(MeasureType measureType) {
        double measuredValue = Double.NaN;

        if(this.isRoot() && !measureType.isAccurateMeasureForRootCell()) // return Double.NaN if root cell and inaccurate measure for roots
            return measuredValue;

        for(ExperimentalMeasure measure : measures) {
            if(measure.measureType == measureType) {
                measuredValue = measure.summaryValue;
                break;
            }
        }

        return measuredValue;
    }



    public List<MeasureType> getMeasureTypes(){
        return this.measureTypes;
    }

    public void addDataPoint(CSVRecord record) {
        for (int i = 0; i < measures.size(); i++) {
            measures.get(i).addDataPointFromRecord(record);
        }
    }

    /**
     * Note: this method relies on the fact that daughter cells are numbered 2n and 2n+1, with n the label of the mother cell.
     * @param cells is the output of LineageTreeParser.parseRawCells()
     * @param rootKey
     * @return
     */
    public static Cell buildTreeAndGetRoot (Map<Integer, Cell> cells, int rootKey) {

        if(cells.size() == 0) throw new IllegalArgumentException("Empty set of cells.");

        Cell rootCell = cells.get(rootKey);

        // set everything relative to node height and branch lengths
        rootCell.setTimePointsAndCleanUpMeasures();
        rootCell.setHeight(rootCell.timePoints.get(rootCell.timePoints.size() - 1));
        rootCell.setEdgeLength();

        // summarize each measure following the appropriate calculation method
        rootCell.summarizeMeasures();

        //TODO rework on how the fate is attributed to the cell
        // only dividers are allowed for now
        rootCell.setFate(Fate.D);

        if(cells.containsKey(2*rootKey)) {
            Cell child1 = buildTreeAndGetRoot(cells, 2*rootKey);
            child1.setParent(rootCell);
            rootCell.addChild(child1);
        }
        if(cells.containsKey(2*rootKey + 1)) {
            Cell child2 = buildTreeAndGetRoot(cells, 2*rootKey + 1);
            child2.setParent(rootCell);
            rootCell.addChild(child2);
        }

        return rootCell;
    }

    /**
     * Renumber all the nodes in the tree following the conventions dictated by Beast2.
     * Leaves must be numbered 0 to n. Internal nodes are n+1 to 2n-2. the root is numbered 2n -1.
     * Warning: this method is written with application to root cells only in mind.
     */
    public void labelNodesInTree() {


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

    /**
     * Try to set timePoints to be the proper time measures,
     * If they don't exist,
     */
    //TODO for now
    //TODO add height to the cell
    //TODO propagate the fact that unscaled time, with a certain frame rate, if no proper time measure.
    public void setTimePointsAndCleanUpMeasures(){
        int unscaledTimeInd = -1;
        int properTimeInd = -1;

        for (int i = 0; i < measures.size(); i++) {
            if(measures.get(i).measureType == MeasureType.ElapsedTime)
                properTimeInd = i;
            else if(measures.get(i).measureType == MeasureType.Lifetime)
                unscaledTimeInd = i;
        }

        if(properTimeInd > -1) { // the proper time measure is present
            // be careful with the order in which elements are removed to avoid problems with renumbering
            if(properTimeInd > unscaledTimeInd) {
//                timePoints = measures.remove(properTimeInd).dataPoints;
//                measures.remove(unscaledTimeInd);

                //For now, the proper time points are never taken into account
                //This is done to simplify things when dealing with cells where only one time point is taken
                //This makes us assume that all time points are separated by the same amount of time
                measures.remove(properTimeInd);
                timePoints = measures.remove(unscaledTimeInd).dataPoints;
            }
            else {
                measures.remove(unscaledTimeInd);
                timePoints = measures.remove(properTimeInd).dataPoints;
            }
        }
        else {
            if(unscaledTimeInd <0) throw new IllegalStateException("No time point column found. This is not compatible with the current implementation.");
            timePoints = measures.remove(unscaledTimeInd).dataPoints;
        }
    }

    public void summarizeMeasures() {
        int inputForSpeedCalculation = 0;

        ExperimentalMeasure xPos = new ExperimentalMeasure();
        ExperimentalMeasure yPos = new ExperimentalMeasure();

        if(timePoints.size() == 0)
            throw new IllegalStateException("Summary is impossible. Time points have not been defined, call setTimePointsAndCleanUpMeasures method first.");

        for (ExperimentalMeasure measure : measures) {
            if(measure.calculationMethod != ExperimentalMeasure.CalculationMethod.undefined &&
                    measure.calculationMethod != ExperimentalMeasure.CalculationMethod.averageInstantSpeed) {
                measure.summarize(timePoints);
            }
            // perform checks to see if we'll have to
            if(measure.measureType == MeasureType.XPosition) {
                xPos = measure;
                inputForSpeedCalculation ++;
            }
            else if(measure.measureType == MeasureType.YPosition){
                yPos = measure;
                inputForSpeedCalculation ++;
            }
        }


        // assuming no duplicates here (only one measure each for xPos and yPos)
        if (inputForSpeedCalculation == 2) {
            for(ExperimentalMeasure measure: measures) {
                if(measure.measureType == MeasureType.InstantSpeed) {
                    this.setInstantSpeedAndSummarize(xPos, yPos, measure);
                    break;
                }
            }
        }
    }

    /**
     * Note: If only one timePoint, no speed is calculated.
     */
    void setInstantSpeedAndSummarize(ExperimentalMeasure xPos, ExperimentalMeasure yPos, ExperimentalMeasure instantSpeed) {
        if((timePoints.size() != xPos.dataPoints.size()) || (yPos.dataPoints.size() != timePoints.size()))
            throw new IllegalArgumentException("All inputs should be of the same size.");

        double dist;
        double time;

        double sumOfSpeeds = 0;

        if(timePoints.size() < 2)  {
            instantSpeed.dataPoints.add(Double.NaN);
            return;
        }

        for (int i = 1; i < timePoints.size(); i++) {

            dist = Math.sqrt(Math.pow(xPos.dataPoints.get(i) - xPos.dataPoints.get(i-1),2)
                    + Math.pow(yPos.dataPoints.get(i) - yPos.dataPoints.get(i-1),2));
            time = timePoints.get(i) - timePoints.get(i-1);

            sumOfSpeeds += dist/time;
            instantSpeed.dataPoints.add(dist/time);
        }

        instantSpeed.summaryValue = sumOfSpeeds/(timePoints.size() -1);
    }


}
