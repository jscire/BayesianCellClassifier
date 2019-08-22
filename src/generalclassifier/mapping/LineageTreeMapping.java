package generalclassifier.mapping;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Loggable;
import generalclassifier.core.LineageTreeProb;
import generalclassifier.lineagetree.CellTree;
import generalclassifier.parametrization.Parametrization;
import generalclassifier.utils.Pair;

import java.io.PrintStream;
import java.util.*;

public class LineageTreeMapping extends CalculationNode implements Loggable {

    public Input <LineageTreeProb> lineageTreeProbInput = new Input<>("lineageTreeProb", "",
            Input.Validate.REQUIRED);

    public Input<CellTree> lineageTreeInput = new Input<>("tree",
            "Lineage tree.",
            Input.Validate.REQUIRED);

    public Input<Parametrization> parametrizationInput = new Input<>("parametrization", "",
            Input.Validate.REQUIRED);


    int numberOfCellTypes;

    Map<Integer, Double[]> storedPruningProb;

    Map<Integer, Integer> mappedCellTypes;

    Random random;

    @Override
    public void initAndValidate() {

        mappedCellTypes = new HashMap<>();
        random = new Random();
        numberOfCellTypes = parametrizationInput.get().numberOfCellTypes;
        paintTree();
    }

    public void paintTree(){

        mappedCellTypes.clear(); // reinitialize painted types on tree

        storedPruningProb = lineageTreeProbInput.get().updateStoredPruningProb();

        int rootTrackNumber = 1;

        if(!lineageTreeInput.get().getLabelsOfAllCellsInTree().contains(rootTrackNumber))
            throw new IllegalStateException("Types on tree cannot be painted as it does not contain a root cell," +
                    " at least not with track number " + rootTrackNumber + ".");

        // draw root type and store it.
        int drawnRootType = drawRootType(rootTrackNumber);
        mappedCellTypes.put(rootTrackNumber, drawnRootType);

        // draw types in rest of the tree.
        paintDaughterCells(rootTrackNumber, drawnRootType);

        if(!mappedCellTypes.keySet().containsAll(lineageTreeInput.get().getLabelsOfAllCellsInTree())) {
            throw new IllegalStateException("A cell in the tree was not painted with a cell type." +
                    " This may mean that there are cells which are not 'connected'" +
                    " to others in the declared set of cells in the tree.");
        }
    }

    void paintDaughterCells(int motherTrackNumber, int motherType){

        int child1TrackNumber = 2 * motherTrackNumber;
        int child2TrackNumber = 2 * motherTrackNumber + 1;

        // check that motherCells has two daughters.
        // if it does, draw types of two daughters.
        if(lineageTreeInput.get().getLabelsOfAllCellsInTree().contains(child1TrackNumber) &&
                lineageTreeInput.get().getLabelsOfAllCellsInTree().contains(child2TrackNumber)) {

            Pair drawnTypes = drawSisterTypes(child1TrackNumber, child2TrackNumber, motherType);

            // store drawn types into mappedCellTypes
            mappedCellTypes.put(child1TrackNumber, drawnTypes.getFirstInt());
            mappedCellTypes.put(child2TrackNumber, drawnTypes.getSecondInt());

            // go to next generation
            paintDaughterCells(child1TrackNumber, drawnTypes.getFirstInt());
            paintDaughterCells(child2TrackNumber, drawnTypes.getSecondInt());
        }

        return;
    }

    int drawRootType(int rootTrackNumber){

        // if root cell type is already fixed, do not draw, just return it.
        if(lineageTreeProbInput.get().getFixedCellType(rootTrackNumber) > -1)
            return lineageTreeProbInput.get().getFixedCellType(rootTrackNumber);

        double[] typeProbs = new double[numberOfCellTypes];
        double sumTypeProbs = 0;

        for (int i = 0; i < numberOfCellTypes; i++) {
            typeProbs[i] = parametrizationInput.get().getTypeFreq(i) * storedPruningProb.get(rootTrackNumber)[i];
            sumTypeProbs += typeProbs[i];
        }

        // draw random value between 0 and sumTypeProbs
        double randValue = random.nextDouble() * sumTypeProbs;
        // find drawn cell type
        for (int i = 0; i < numberOfCellTypes; i++) {
            if(randValue <= typeProbs[i])
                return i;
            else
                randValue -= typeProbs[i];
        }

        return -1;
    }

    Pair drawSisterTypes(int child1TrackNumber, int child2TrackNumber, int motherType){

        int child1Type = lineageTreeProbInput.get().getFixedCellType(child1TrackNumber);
        int child2Type = lineageTreeProbInput.get().getFixedCellType(child2TrackNumber);

        if(child1Type > -1 && child2Type > -1) // daughter types are already fixed, no need to draw them.
            return new Pair(child1Type, child2Type);

        LinkedHashMap<Pair, Double> intermediateProbs = new LinkedHashMap<>();
        double sumTypeProbs = 0;

        if(child1Type == -1 && child2Type == -1) { // neither daughter types are fixed, draw both of them.

            for (int j = 0; j < numberOfCellTypes; j++) {
                for (int k = 0; k < numberOfCellTypes; k++) {

                    Double typeProb = storedPruningProb.get(child1TrackNumber)[j] *
                            storedPruningProb.get(child2TrackNumber)[k] *
                            parametrizationInput.get().getTransitionProbability(motherType, j, k);

                    intermediateProbs.put(new Pair(j,k), typeProb);

                    sumTypeProbs += typeProb;
                }
            }

        } else if(child1Type > -1 && child2Type == -1) { // fixed type for just one of the two children

            for (int k = 0; k < numberOfCellTypes; k++) {
                Double typeProb = storedPruningProb.get(child1TrackNumber)[child1Type] *
                        storedPruningProb.get(child2TrackNumber)[k] *
                        parametrizationInput.get().getTransitionProbability(motherType, child1Type, k);

                intermediateProbs.put(new Pair(child1Type,k), typeProb);

                sumTypeProbs += typeProb;
            }

        } else if(child1Type == -1 && child2Type > -1) { // fixed type for just one of the two children

            for (int j = 0; j < numberOfCellTypes; j++) {

                Double typeProb = storedPruningProb.get(child1TrackNumber)[j] *
                        storedPruningProb.get(child2TrackNumber)[child2Type] *
                        parametrizationInput.get().getTransitionProbability(motherType, j, child2Type);

                intermediateProbs.put(new Pair(j, child2Type), typeProb);

                sumTypeProbs += typeProb;
            }
        } else {
            throw new IllegalStateException("Child1's and/or Child2's type are/is fixed to invalid value(s).");
        }

        // draw random value between 0 and sumTypeProbs
        double randValue = random.nextDouble() * sumTypeProbs;
        // find drawn cell type
        for (Pair p : intermediateProbs.keySet()) {
            if(randValue <= intermediateProbs.get(p))
                return p;
            else
                randValue -= intermediateProbs.get(p);
        }

        return null;
    }

    /**
     * Loggable interface implementation follows.
     */

    @Override
    public void init(final PrintStream out) {
        int treeIdx = lineageTreeProbInput.get().treeIdxInput.get();

        for (Integer label : lineageTreeInput.get().getLabelsOfAllCellsInTree()) {
            out.print("cellType" + treeIdx + "_" + label + "\t");
        }
    }

    @Override
    public void log(long sampleNr, PrintStream out) {
        paintTree();
        for (Integer label : lineageTreeInput.get().getLabelsOfAllCellsInTree()) {
            out.print(mappedCellTypes.get(label) + "\t");
        }
    }

    @Override
    public void close(final PrintStream out) {
        // nothing to do
    }

}
