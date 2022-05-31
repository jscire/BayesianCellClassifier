package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import generalclassifier.core.LineageTreeProb;
import generalclassifier.lineagetree.LineageTree;
import generalclassifier.mapping.LineageTreeMapping;
import generalclassifier.parametrization.DistributionForMeasurement;
import generalclassifier.parametrization.ExperimentalMeasurements;
import generalclassifier.parametrization.Parametrization;
import junit.framework.TestCase;
import org.junit.Test;

import java.util.LinkedList;
import java.util.List;


public class LikelihoodTests extends TestCase{


    @Test
    /**
     * Reference values come from R calculations ('by hand')
     */
    public void testBasicLikelihood_1() throws Exception {
        ////////////////////////////////// Tree with 1 cell, 1 normal measure
        LineageTree tree = new LineageTree();

        LineageTreeProb treeProb = new LineageTreeProb();

        DistributionForMeasurement distr_measure1 = new DistributionForMeasurement();

        distr_measure1.initByName("measurementTag", "measure_1",
                "parm1Distribution", new RealParameter("0.5 -0.5"),
                "parm2Distribution", new RealParameter("1.0 1.0"),
                "distributionType", "normal",
                "estimateType", "mean",
                "isAppliedToRootCells", true);
        distr_measure1.initAndValidate();

        Parametrization parametrization = new Parametrization();

        parametrization.initByName("distribution", distr_measure1,
                "transitionUponDivisionProbs", new RealParameter("0.5 0.4 0.1"),
                "transitionUponDivisionProbs", new RealParameter("0 0 1"));

        ExperimentalMeasurements measure_1 = new ExperimentalMeasurements();
        measure_1.initByName("measurementTag", "measure_1", "values", "1:1.0");
        tree.setInputValue("measurement", measure_1);

        tree.setInputValue("cellsInTree", "1");
        tree.setInputValue("cellsAreFullyTracked", "true");

        tree.initAndValidate();

        treeProb.setInputValue("tree", tree);
        treeProb.setInputValue("parametrization", parametrization);
        treeProb.setInputValue("cellType", new IntegerParameter("0"));
        treeProb.setInputValue("rootTypeOnly", "true");

        treeProb.initAndValidate();

        double logP = treeProb.calculateLogP();

        assertEquals(logP, -1.737086, 1e-4);
    }

    @Test
    /**
     * Reference values come from R calculations ('by hand')
     */
    public void testBasicLikelihood_2() throws Exception {

        /// Tree with 1 cell, 3 normal measures including lifetime and one measure excluded from root
        LineageTree tree = new LineageTree();

        LineageTreeProb treeProb = new LineageTreeProb();

        ExperimentalMeasurements lifetime = new ExperimentalMeasurements();
        lifetime.initByName("measurementTag", "lifetime", "values", "1:-0.5");
        ExperimentalMeasurements measure_1 = new ExperimentalMeasurements();
        measure_1.initByName("measurementTag", "measure_1", "values", "1:0.1");
        ExperimentalMeasurements measure_2 = new ExperimentalMeasurements();
        measure_2.initByName("measurementTag", "measure_2", "values", "1:1.3");

        List<ExperimentalMeasurements> experimentalMeasurements = new LinkedList<>();
        experimentalMeasurements.add(lifetime);
        experimentalMeasurements.add(measure_1);
        experimentalMeasurements.add(measure_2);

        DistributionForMeasurement distr_lifetime = new DistributionForMeasurement();

        distr_lifetime.initByName("measurementTag", "lifetime",
                "parm1Distribution", new RealParameter("-0.5 0.5"),
                "parm2Distribution", new RealParameter("1.0 0.9"),
                "distributionType", "normal",
                "estimateType", "max",
                "isAppliedToRootCells", true);
        distr_lifetime.initAndValidate();

        DistributionForMeasurement distr_measure1 = new DistributionForMeasurement();

        distr_measure1.initByName("measurementTag", "measure_1",
                "parm1Distribution", new RealParameter("-0.4 -0.5"),
                "parm2Distribution", new RealParameter("0.4 0.3"),
                "distributionType", "normal",
                "estimateType", "mean",
                "isAppliedToRootCells", true);
        distr_measure1.initAndValidate();

        DistributionForMeasurement distr_measure2 = new DistributionForMeasurement();

        distr_measure2.initByName("measurementTag", "measure_2",
                "parm1Distribution", new RealParameter("-0.9 -0.9"),
                "parm2Distribution", new RealParameter("0.6 0.5"),
                "distributionType", "normal",
                "estimateType", "mean",
                "isAppliedToRootCells", false);
        distr_measure2.initAndValidate();

        List<DistributionForMeasurement> distributions = new LinkedList<>();
        distributions.add(distr_lifetime);
        distributions.add(distr_measure1);
        distributions.add(distr_measure2);

        Parametrization parametrization  = new Parametrization();

        parametrization.initByName("distribution", distributions,
                "transitionUponDivisionProbs", new RealParameter("0.5 0.4 0.1"),
                "transitionUponDivisionProbs", new RealParameter("0 0 1"));

        tree.setInputValue("measurement", experimentalMeasurements);

        tree.setInputValue("cellsInTree", "1");
        tree.initAndValidate();

        treeProb.setInputValue("tree", tree);
        treeProb.setInputValue("parametrization", parametrization);
        treeProb.setInputValue("cellType", new IntegerParameter("1"));
        treeProb.setInputValue("rootTypeOnly", "true");

        treeProb.initAndValidate();

        double logP = treeProb.calculateLogP();
        assertEquals( -2.551129, logP, 1e-4);

        treeProb.setInputValue("cellType", new IntegerParameter("0"));
        treeProb.setInputValue("rootTypeOnly", "true");

        treeProb.initAndValidate();

        logP = treeProb.calculateLogP();
        assertEquals(-2.170192, logP, 1e-4);
    }


    @Test
    /**
     * Reference values come from R calculations ('by hand')
     */
    public void testBasicLikelihood_3() throws Exception {

        /// Tree with 3 cells, 4 normal measures including lifetime
        /// one measure excluded from root
        /// one measure only for roots

        LineageTree tree = new LineageTree();

        LineageTreeProb treeProb = new LineageTreeProb();

        ExperimentalMeasurements lifetime = new ExperimentalMeasurements();
        lifetime.initByName("measurementTag", "lifetime", "values", "1:-2.3,2:1.1,3:0.1");
        ExperimentalMeasurements measure_1 = new ExperimentalMeasurements();
        measure_1.initByName("measurementTag", "measure_1", "values", "1:0.1,2:0.5,3:-0.3");
        ExperimentalMeasurements measure_2 = new ExperimentalMeasurements();
        measure_2.initByName("measurementTag", "measure_2", "values", "1:1.3,2:2.1,3:-0.1");
        ExperimentalMeasurements measure_3 = new ExperimentalMeasurements();
        measure_3.initByName("measurementTag", "measure_3", "values", "1:-0.1,2:0.5,3:-1.1");
        ExperimentalMeasurements measure_4 = new ExperimentalMeasurements();
        measure_4.initByName("measurementTag", "measure_4", "values", "1:-0.6,2:0.2,3:0.3");

        List<ExperimentalMeasurements> experimentalMeasurements = new LinkedList<>();
        experimentalMeasurements.add(lifetime);
        experimentalMeasurements.add(measure_1);
        experimentalMeasurements.add(measure_2);
        experimentalMeasurements.add(measure_3);
        experimentalMeasurements.add(measure_4);

        DistributionForMeasurement distr_lifetime = new DistributionForMeasurement();

        distr_lifetime.initByName("measurementTag", "lifetime",
                "parm1Distribution", new RealParameter("0.1 0.4"),
                "parm2Distribution", new RealParameter("1.0 0.9"),
                "distributionType", "normal",
                "estimateType", "max",
                "isAppliedToRootCells", true);
        distr_lifetime.initAndValidate();

        DistributionForMeasurement distr_measure1 = new DistributionForMeasurement();

        distr_measure1.initByName("measurementTag", "measure_1",
                "parm1Distribution", new RealParameter("-0.4 -0.1"),
                "parm2Distribution", new RealParameter("0.4 0.3"),
                "distributionType", "normal",
                "estimateType", "mean",
                "isAppliedToRootCells", true);
        distr_measure1.initAndValidate();

        DistributionForMeasurement distr_measure2 = new DistributionForMeasurement();

        distr_measure2.initByName("measurementTag", "measure_2",
                "parm1Distribution", new RealParameter("-0.9 -0.3"),
                "parm2Distribution", new RealParameter("0.6 0.5"),
                "distributionType", "normal",
                "estimateType", "mean",
                "isAppliedToRootCells", false);
        distr_measure2.initAndValidate();

        DistributionForMeasurement distr_measure3 = new DistributionForMeasurement();

        distr_measure3.initByName("measurementTag", "measure_3",
                "parm1Distribution", new RealParameter("-0.5 -0.2"),
                "parm2Distribution", new RealParameter("0.5 0.8"),
                "distributionType", "normal",
                "estimateType", "mean",
                "isAppliedToRootCells", true);
        distr_measure3.initAndValidate();

        DistributionForMeasurement distr_measure4 = new DistributionForMeasurement();

        distr_measure4.initByName("measurementTag", "measure_4",
                "parm1Distribution", new RealParameter("1.1 -0.5"),
                "parm2Distribution", new RealParameter("2.0 1.5"),
                "distributionType", "normal",
                "estimateType", "mean",
                "isAppliedToRootCells", true,
                "isAppliedToRootCellsOnly", true);
        distr_measure4.initAndValidate();

        List<DistributionForMeasurement> distributions = new LinkedList<>();
        distributions.add(distr_lifetime);
        distributions.add(distr_measure1);
        distributions.add(distr_measure2);
        distributions.add(distr_measure3);
        distributions.add(distr_measure4);

        Parametrization parametrization  = new Parametrization();

        parametrization.initByName("distribution", distributions,
                "transitionUponDivisionProbs", new RealParameter("0.5 0.4 0.1"),
                "transitionUponDivisionProbs", new RealParameter("0 0 1"));

        tree.setInputValue("measurement", experimentalMeasurements);

        tree.setInputValue("cellsInTree", "1,2,3");
        tree.setInputValue("cellsAreFullyTracked", "true");

        tree.initAndValidate();

        //// celltypes 0,0,0
        treeProb.setInputValue("tree", tree);
        treeProb.setInputValue("parametrization", parametrization);
        treeProb.setInputValue("cellType", new IntegerParameter("0 0 0"));
        treeProb.setInputValue("rootTypeOnly", "false");

        treeProb.initAndValidate();

        double logP = treeProb.calculateLogP();
        assertEquals(logP, -26.97992, 1e-4);


        //// celltypes 1,1,1
        treeProb.setInputValue("tree", tree);
        treeProb.setInputValue("parametrization", parametrization);
        treeProb.setInputValue("cellType", new IntegerParameter("1 1 1"));
        treeProb.setInputValue("rootTypeOnly", "false");

        treeProb.initAndValidate();

        logP = treeProb.calculateLogP();
        assertEquals(logP, -20.75805, 1e-4);

        //// celltypes 0,1,0
        treeProb.setInputValue("tree", tree);
        treeProb.setInputValue("parametrization", parametrization);
        treeProb.setInputValue("cellType", new IntegerParameter("0 1 0"));
        treeProb.setInputValue("rootTypeOnly", "false");

        treeProb.initAndValidate();

        logP = treeProb.calculateLogP();
        assertEquals(logP, -24.46488, 1e-4);

        //// celltypes 1,0,1
        parametrization  = new Parametrization();

        parametrization.initByName("distribution", distributions,
                "transitionUponDivisionProbs", new RealParameter("0.5 0.4 0.1"),
                "transitionUponDivisionProbs", new RealParameter("0.3 0.2 0.5"));

        treeProb.setInputValue("tree", tree);
        treeProb.setInputValue("parametrization", parametrization);
        treeProb.setInputValue("cellType", new IntegerParameter("1 0 1"));
        treeProb.setInputValue("rootTypeOnly", "false");

        treeProb.initAndValidate();

        logP = treeProb.calculateLogP();
        assertEquals(logP, -26.49197, 1e-4);

        //// cells are not fully tracked
        parametrization  = new Parametrization();

        parametrization.initByName("distribution", distributions,
                "transitionUponDivisionProbs", new RealParameter("0.5 0.4 0.1"),
                "transitionUponDivisionProbs", new RealParameter("0.3 0.2 0.5"));

        tree.setInputValue("cellsAreFullyTracked", "false");
        tree.initAndValidate();

        treeProb.setInputValue("tree", tree);
        treeProb.setInputValue("parametrization", parametrization);
        treeProb.setInputValue("cellType", new IntegerParameter("1 0 1"));
        treeProb.setInputValue("rootTypeOnly", "false");

        treeProb.initAndValidate();

        logP = treeProb.calculateLogP();
        assertEquals(logP, -26.50606, 1e-4);

        //// root type only, root type 1
        parametrization  = new Parametrization();

        parametrization.initByName("distribution", distributions,
                "transitionUponDivisionProbs", new RealParameter("0.5 0.4 0.1"),
                "transitionUponDivisionProbs", new RealParameter("0.3 0.2 0.5"));

        tree.setInputValue("cellsAreFullyTracked", "true");
        tree.initAndValidate();

        treeProb.setInputValue("tree", tree);
        treeProb.setInputValue("parametrization", parametrization);
        treeProb.setInputValue("cellType", new IntegerParameter("1 0 1"));
        treeProb.setInputValue("rootTypeOnly", "true");

        treeProb.initAndValidate();
        logP = treeProb.calculateLogP();
        assertEquals(logP, -21.3472, 1e-4);

        //// root type only, root type 0
        parametrization  = new Parametrization();

        parametrization.initByName("distribution", distributions,
                "transitionUponDivisionProbs", new RealParameter("0.5 0.4 0.1"),
                "transitionUponDivisionProbs", new RealParameter("0.3 0.2 0.5"));

        tree.setInputValue("cellsAreFullyTracked", "true");
        tree.initAndValidate();

        treeProb.setInputValue("tree", tree);
        treeProb.setInputValue("parametrization", parametrization);
        treeProb.setInputValue("cellType", new IntegerParameter("0 0 1"));
        treeProb.setInputValue("rootTypeOnly", "true");

        treeProb.initAndValidate();
        logP = treeProb.calculateLogP();
        assertEquals(logP, -23.67036, 1e-4);
    }

    @Test
    /**
     * Reference values from package itself, to check internal consistency
     */
    public void testBasicLikelihood_4() throws Exception {

        /// Tree with 7 cells, 3 normal measures, one is lifetime (incomplete for root and tips)

        LineageTree tree = new LineageTree();

        LineageTreeProb treeProb = new LineageTreeProb();

        ExperimentalMeasurements lifetime = new ExperimentalMeasurements();
        lifetime.initByName("measurementTag", "lifetime", "values", "1:-2.3,2:1.1,3:0.1,4:-0.5,5:-0.8,6:-0.2,7:-0.3");
        ExperimentalMeasurements measure_1 = new ExperimentalMeasurements();
        measure_1.initByName("measurementTag", "measure_1", "values", "1:0.1,2:0.5,3:-0.3,4:0.89,5:-0.1,6:-1.3,7:1.5");
        ExperimentalMeasurements measure_2 = new ExperimentalMeasurements();
        measure_2.initByName("measurementTag", "measure_2", "values", "1:1.3,2:2.1,3:-0.1,4:-0.45,5:0.59,6:-2.5,7:0.4");
        ExperimentalMeasurements measure_3 = new ExperimentalMeasurements();
        measure_3.initByName("measurementTag", "measure_3", "values", "1:-0.1,2:0.5,3:-1.1,4:-0.2,5:0.8,6:1.2,7:0.33");

        List<ExperimentalMeasurements> experimentalMeasurements = new LinkedList<>();
        experimentalMeasurements.add(lifetime);
        experimentalMeasurements.add(measure_1);
        experimentalMeasurements.add(measure_2);
        experimentalMeasurements.add(measure_3);

        DistributionForMeasurement distr_lifetime = new DistributionForMeasurement();

        distr_lifetime.initByName("measurementTag", "lifetime",
                "parm1Distribution", new RealParameter("0.1 0.4"),
                "parm2Distribution", new RealParameter("1.0 0.9"),
                "distributionType", "normal",
                "estimateType", "max",
                "isAppliedToRootCells", true);
        distr_lifetime.initAndValidate();

        DistributionForMeasurement distr_measure1 = new DistributionForMeasurement();

        distr_measure1.initByName("measurementTag", "measure_1",
                "parm1Distribution", new RealParameter("-0.4 -0.1"),
                "parm2Distribution", new RealParameter("0.4 0.3"),
                "distributionType", "normal",
                "estimateType", "mean",
                "isAppliedToRootCells", true);
        distr_measure1.initAndValidate();

        DistributionForMeasurement distr_measure2 = new DistributionForMeasurement();

        distr_measure2.initByName("measurementTag", "measure_2",
                "parm1Distribution", new RealParameter("-0.9 0.3"),
                "parm2Distribution", new RealParameter("0.6 0.5"),
                "distributionType", "normal",
                "estimateType", "mean",
                "isAppliedToRootCells", false);
        distr_measure2.initAndValidate();

        DistributionForMeasurement distr_measure3 = new DistributionForMeasurement();

        distr_measure3.initByName("measurementTag", "measure_3",
                "parm1Distribution", new RealParameter("-0.5 0.7"),
                "parm2Distribution", new RealParameter("0.5 0.8"),
                "distributionType", "normal",
                "estimateType", "mean",
                "isAppliedToRootCells", true);
        distr_measure3.initAndValidate();


        List<DistributionForMeasurement> distributions = new LinkedList<>();
        distributions.add(distr_lifetime);
        distributions.add(distr_measure1);
        distributions.add(distr_measure2);
        distributions.add(distr_measure3);

        Parametrization parametrization = new Parametrization();

        parametrization.initByName("distribution", distributions,
                "haveGenerationSpecificTransitionProbs", true,
                "transitionUponDivisionProbs", new RealParameter("0.8 0.1 0.1"),
                "transitionUponDivisionProbs", new RealParameter("0 0 1"),
                "transitionUponDivisionProbs", new RealParameter("0.2 0.1 0.7"),
                "transitionUponDivisionProbs", new RealParameter("0 0 1"));

        tree.setInputValue("measurement", experimentalMeasurements);

        tree.setInputValue("cellsInTree", "1,2,3,4,5,6,7");
        tree.setInputValue("cellsAreFullyTracked", "false");

        tree.initAndValidate();

        //// celltypes 0,0,0,0,1,1,1
        treeProb.setInputValue("tree", tree);
        treeProb.setInputValue("parametrization", parametrization);
        treeProb.setInputValue("cellType", new IntegerParameter("0 0 0 0 1 1 1"));
        treeProb.setInputValue("rootTypeOnly", "false");

        treeProb.initAndValidate();

        double logP = treeProb.calculateLogP();
        assertEquals(logP, -75.45081537472011, 1e-4);


        //mapping
        treeProb.setInputValue("cellType", new IntegerParameter("-1 -1 -1 -1 -1 -1 -1"));
        treeProb.initAndValidate();

        LineageTreeMapping mapping = new LineageTreeMapping();

        mapping.setInputValue("lineageTreeProb", treeProb);
        mapping.setInputValue("tree", tree);
        mapping.setInputValue("parametrization", parametrization);

        mapping.initAndValidate();
        for (int i = 0; i < 10; i++) {
            System.out.println(mapping.printMapping(true));
        }
    }
}
