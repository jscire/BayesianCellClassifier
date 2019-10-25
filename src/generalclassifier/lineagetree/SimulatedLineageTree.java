package generalclassifier.lineagetree;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import generalclassifier.parametrization.DistributionForMeasurement;
import generalclassifier.parametrization.Parametrization;

import java.util.SortedSet;
import java.util.TreeSet;

public class SimulatedLineageTree extends CellTree {

    public Input<Parametrization> parametrizationInput = new Input<>("parametrization",
            "",
            Input.Validate.REQUIRED);

    public Input<Integer> lastGenerationInput = new Input<>("lastGeneration",
            "",
            Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() {
        if(lastGenerationInput.get() < 1)
            throw new IllegalArgumentException("lastGenerationOfInterest must be at least 1.");

        uniqueMeasurementTags = parametrizationInput.get().getMeasurementTags();

        this.simulate();
        initArrays();

        labelsOfAllCellsInTree = buildSetOfCellLabels((SimulatedCell) this.getRoot());

        super.initAndValidate();
    }

    private void simulate(){

        int rootType = drawRootCellType();
        SimulatedCell rootCell = new SimulatedCell(1, rootType);
        rootCell.simulateCell(parametrizationInput.get(), lastGenerationInput.get());
        rootCell.labelNodesInTree();
        setRoot(rootCell);
    }

    private int drawRootCellType(){
        double rand = Math.random();
        int cellType = 0;
        while(cellType < (parametrizationInput.get().getNumberOfCellTypes() -1)){
            if(rand < parametrizationInput.get().getTypeFreq(cellType))
                break;
            else {
                rand -= parametrizationInput.get().getTypeFreq(cellType);
                cellType ++;
            }
        }

        return cellType;
    }

    @Override
    public String getHeaderOfCSV(){
        String res="TrackNumber,Fate,Type";
        for(String tag : uniqueMeasurementTags) {
            res += "," + tag;
        }
        return res + "";
    }

    public SortedSet<Integer> buildSetOfCellLabels(Cell cell){
        SortedSet<Integer> resultSubTree = new TreeSet<>();
        resultSubTree.add(cell.getTrackNumber());
        if(cell.isLeaf())
            return resultSubTree;
        else {
            Cell child1 = (Cell) cell.getChild(0);
            Cell child2 = (Cell) cell.getChild(1);
            resultSubTree.addAll(buildSetOfCellLabels(child1));
            resultSubTree.addAll(buildSetOfCellLabels(child2));
            return resultSubTree;
        }
    }

    public SortedSet<Integer> getLabelsOfAllCellsInTree(){
        return labelsOfAllCellsInTree;
    };

    public static void main(String[] parms) {

        SimulatedLineageTree tree =  new SimulatedLineageTree();

        Parametrization parametrization = new Parametrization();

        DistributionForMeasurement distributionSca1 = new DistributionForMeasurement();
        distributionSca1.initByName("parm1Distribution", new RealParameter("2 -0.1"),
                "parm2Distribution", new RealParameter("0.5 0.4"),
                "measurementTag", "Sca1",
                "distributionType", "normal");

        DistributionForMeasurement distributionCD71 = new DistributionForMeasurement();
        distributionCD71.initByName("parm1Distribution", new RealParameter("-0.5 0.4"),
                "parm2Distribution", new RealParameter("1.5 0.9"),
                "measurementTag", "CD71",
                "distributionType", "normal");

        DistributionForMeasurement distributionArea = new DistributionForMeasurement();
        distributionArea.initByName("parm1Distribution", new RealParameter("-1.3 0.9"),
                "parm2Distribution", new RealParameter("0.3 0.28"),
                "measurementTag", "Area",
                "distributionType", "normal",
                "isAppliedToRootCells", "false");


        DistributionForMeasurement distributionLifetime = new DistributionForMeasurement();
        distributionLifetime.initByName("parm1Distribution", new RealParameter("29 25"),
                "parm2Distribution", new RealParameter("3.4 5.1"),
                "measurementTag", "Lifetime",
                "distributionType", "normal",
                "estimateType", "min");

        distributionSca1.initAndValidate();
        distributionCD71.initAndValidate();
        distributionArea.initAndValidate();
        distributionLifetime.initAndValidate();

        parametrization.setInputValue("distribution", distributionSca1);
        parametrization.setInputValue("distribution", distributionCD71);
        parametrization.setInputValue("distribution", distributionArea);
        parametrization.setInputValue("distribution", distributionLifetime);

        parametrization.initByName(
//                "transitionUponDivisionProbs", new RealParameter("0.6 0.15 0.25"),
                "transitionUponDivisionProbs", new RealParameter("0.0 1.0 0.0"),
                "transitionUponDivisionProbs", new RealParameter("0.0 1.0 0.0"));

        tree.setInputValue("parametrization", parametrization);

        tree.initAndValidate();

        System.out.println("Done");
//
//        tree.initAndValidate();
//
//        System.out.println(tree.toString());
    }
}
