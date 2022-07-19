package BayesianCellClassifier.lineagetree;

import BayesianCellClassifier.parametrization.DistributionForMeasurement;
import BayesianCellClassifier.parametrization.Parametrization;
import BayesianCellClassifier.utils.Utils;

import java.util.HashMap;
import java.util.List;
import java.util.SortedSet;

public class SimulatedCell extends Cell {

    private int cellType;

    public SimulatedCell(int trackNumber) {
        super(trackNumber);
        this.cellType = 0;
    }

    public SimulatedCell(int trackNumber, int cellType) {
        super(trackNumber);
        this.cellType = cellType;
    }

    public int getCellType() {
        return cellType;
    }

    public void setCellType(int cellType){
        this.cellType = cellType;
    }

    /**
     * Simulate experimental measures under distributions defined in paramterizationInput
     * Note that we do not test whether a cell is a root cell before simulating.
     * So, for root cells we may end up simulating some useless values, (for measures not applied to root cells).
     */
    public void simulateExperimentalMeasures(Parametrization parametrizationInput, int cellType){
        //wipe previous experimentalMeasures attached to this cell.
        experimentalMeasures = new HashMap<>();

        List<DistributionForMeasurement> distributionsForMeasures = parametrizationInput.getDistributions();

        for(DistributionForMeasurement distribution : distributionsForMeasures) {

            if(distribution.getHasZeroFraction() && Math.random() < distribution.getZeroFraction(cellType)) {
                // is simulated value zero?
                experimentalMeasures.put(distribution.getMeasurementTag(), 0.0);
            }
            else {
                double drawnValue = Utils.getRandomValueFromDistribution(
                        distribution.getDistributionType(),
                        distribution.getParm1(cellType),
                        distribution.getParm2(cellType));

                experimentalMeasures.put(distribution.getMeasurementTag(), drawnValue);
            }

        }
    }

    public void simulateCell(Parametrization parametrization, int maxGenerationNumber) {

        this.simulateExperimentalMeasures(parametrization, cellType);

        //Draw whether cell was lost or not
        double rand = Math.random();
        if(rand < parametrization.getLossProbability()) {
            this.setFate(Cell.Fate.L);
            return;
        }

        //Draw cell fate (either divides or dies)
        double[] fateProbs = parametrization.getFateProbability(cellType);
        rand = Math.random();
        int idx = 0;
        while(idx < (fateProbs.length - 1)) {
            if(rand < fateProbs[idx])
                break;
            else {
                rand -= fateProbs[idx];
                idx ++;
            }
        }

        if(idx == 0)
            this.setFate(Cell.Fate.D);
        else if (idx == 1) {
            this.setFate(Cell.Fate.A);
            return;
        }
        else
            throw new IllegalArgumentException("Fate of idx " + idx + " is not supported in simulation yet.");

        // if here, cell divides

        // check if we are at the max generation
        if(Math.floor(Math.log(2*this.trackNumber)/Math.log(2.0)) >= maxGenerationNumber) return;

        // draw daughters types
        int[] daughterTypes  = drawDaughterCellTypes(parametrization, cellType);
        int idxTypeChild1 = Math.random()  < 0.5 ? 0 : 1; // shuffle which of the two daughters has which of the two types
        int idxTypeChild2 = idxTypeChild1 ^ 1; // idxTypeChild is 0 if idxTypeChild1 is 1 and vice-versa

        // simulate first child
        SimulatedCell child1 = new SimulatedCell(2 * this.trackNumber, daughterTypes[idxTypeChild1]);
        child1.simulateCell(parametrization, maxGenerationNumber);
        //simulate second child
        SimulatedCell child2 = new SimulatedCell(2 * this.trackNumber + 1, daughterTypes[idxTypeChild2]);
        child2.simulateCell(parametrization, maxGenerationNumber);

        child1.setParent(this);
        child2.setParent(this);

        this.addChild(child1);
        this.addChild(child2);
    }


    private int[] drawDaughterCellTypes(Parametrization parametrization, int motherCellType){

        double rand = Math.random();
        double[] transitionProbs = parametrization.getTransitionProbabilitiesForMotherType(motherCellType);

        int idx = 0;
        int typeChild1 = 0;
        int typeChild2 = 0;

        while(idx < (transitionProbs.length - 1)){

            if(rand < transitionProbs[idx])
                break;
            else{
                rand -= transitionProbs[idx];
                idx ++;

                if(typeChild2 < (parametrization.getNumberOfCellTypes() -1))
                    typeChild2 ++;
                else {
                    typeChild1 ++;
                    typeChild2=typeChild1;
                }
            }
        }
        return new int[]{typeChild1, typeChild2};
    }


    @Override
    public String toCSVRecord(SortedSet<String> sortedTags){

        String res = "";
        res += trackNumber + "," + fate.toString() + ","  + cellType;
        for(String tag : sortedTags) {
            res += "," + numberFormat.format(experimentalMeasures.get(tag));
        }
        return res;
    }

}
