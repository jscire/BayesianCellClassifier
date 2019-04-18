package test;

import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import com.sun.org.apache.regexp.internal.RE;
import generalclassifier.core.LineageTreeProb;
import generalclassifier.lineagetree.LineageTree;
import junit.framework.TestCase;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;


public class LikelihoodTests extends TestCase{

    //TODO obsolete tests remake unit tests
    //TODO rework on all the tests are we don't take newick inputs anymore
    //TODO the reference values are modified compared to the SCclassify tests. Check if it makes sense.
    // The prob is multiplied by .5 at the root to reflect the expected purity of the pop and by .5 at each (non lost) leaf to reflect the naked expectancy we have on each cell being a certain type (not sure about that part)
    // I left the + k log(0.5) so that it's easily changeable later.

    @Test
    /**
     * Reference values come from calculations "by hand" in R
     */
    public void testBasicLikelihood() throws Exception{
        LineageTree tree = new LineageTree();
        LineageTreeProb probTree = new LineageTreeProb();

        List<RealParameter> fateProbabilities = new ArrayList<>();
        fateProbabilities.add(new RealParameter("0.4 0.4 0.2"));
        fateProbabilities.add(new RealParameter("0.3 0.3 0.4"));

        List<RealParameter> scaleWeibull = new ArrayList<>();
        scaleWeibull.add(new RealParameter("2 2"));
        scaleWeibull.add(new RealParameter("3 3"));

        List<RealParameter> shapeWeibull = new ArrayList<>();
        shapeWeibull.add(new RealParameter("1 1"));
        shapeWeibull.add(new RealParameter("1 1"));

        List<RealParameter> probsTransitionUponDivision = new ArrayList<>();
        probsTransitionUponDivision.add(new RealParameter("0.5 0.3 0.2"));
        probsTransitionUponDivision.add(new RealParameter("0 0 1.0"));

        probTree.setInputValue("probsOfTransitionUponDivision", probsTransitionUponDivision);
        probTree.setInputValue("fateProbabilities", fateProbabilities);
        probTree.setInputValue("scaleWeibull",scaleWeibull);
        probTree.setInputValue("shapeWeibull",shapeWeibull);
        probTree.setInputValue("lossProb", new RealParameter("0.1"));
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));

        double logP;

        // Basic tree: 1A:34;
        tree.setInputValue("newick", "1A:34;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -18.7148 + 2*Math.log(0.5), 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -13.74128 + 2*Math.log(0.5), 1e-4);


        // Basic tree: 1N:26;
        tree.setInputValue("newick", "1N:26;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -1.714789 + 2*Math.log(0.5), 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -1.021393 + 2*Math.log(0.5), 1e-4);


        // Basic tree: 1L:22;
        tree.setInputValue("newick", "1L:22;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.302585 + Math.log(0.5), 1e-4);

        // Tree with 2 leaves: (2L:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2L:4,3N:3)1D:2;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -5.611862 + 2*Math.log(0.5), 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -5.959421 + 2*Math.log(0.5), 1e-4);

        // Tree with 2 leaves: (2A:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2A:4,3N:3)1D:2;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -7.034257 + 3*Math.log(0.5), 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -7.398115 + 3* Math.log(0.5), 1e-4);

        // Tree with 4 leaves: ((4A:7,5A:2)2D:1,(6D:6,7A:9)3D:3)1D:4;
        tree.setInputValue("newick", "((4A:7,5A:2)2D:1,(6D:6,7A:9)3D:3)1D:4;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -26.78519 + 5* Math.log(0.5), 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -27.52229 + 5* Math.log(0.5), 1e-4);

        // Tree with 4 leaves: ((4L:7,5A:2)2D:1,(6N:6,7L:9)3D:3)1D:4;
        tree.setInputValue("newick", "((4L:7,5A:2)2D:1,(6N:6,7L:9)3D:3)1D:4;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -17.32192 + 3*Math.log(0.5), 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -18.40712 + 3*Math.log(0.5), 1e-4);
    }


    @Test
    /**
     * Reference values come from calculations "by hand" in R
     */
    public void testLikelihoodWithIntegrationOverTransitionTime() throws Exception{

        LineageTree tree = new LineageTree();
        LineageTreeProb probTree = new LineageTreeProb();

        List<RealParameter> fateProbabilities = new ArrayList<>();
        fateProbabilities.add(new RealParameter("0.6 0.2 0.1 0.1"));
        fateProbabilities.add(new RealParameter("0.3 0.3 0.4"));

        List<RealParameter> scaleWeibull = new ArrayList<>();
        scaleWeibull.add(new RealParameter("2 2.2 2.5"));
        scaleWeibull.add(new RealParameter("3 3.1"));

        List<RealParameter> shapeWeibull = new ArrayList<>();
        shapeWeibull.add(new RealParameter("1 1.5 1.3"));
        shapeWeibull.add(new RealParameter("1 1.2"));

        List<RealParameter> probsTransitionUponDivision = new ArrayList<>();
        probsTransitionUponDivision.add(new RealParameter("0.5 0.3 0.2"));
        probsTransitionUponDivision.add(new RealParameter("0 0 1.0"));

        probTree.setInputValue("probsOfTransitionUponDivision", probsTransitionUponDivision);
        probTree.setInputValue("fateProbabilities", fateProbabilities);
        probTree.setInputValue("scaleWeibull",scaleWeibull);
        probTree.setInputValue("shapeWeibull",shapeWeibull);
        probTree.setInputValue("lossProb", new RealParameter("0.1"));
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.setInputValue("allowedTransitions", new BooleanParameter("1 0"));

        double logP;

        // Basic tree: 1A:2;
        tree.setInputValue("newick", "1A:2;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.975158 + 2*Math.log(0.5), 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.937085 + 2*Math.log(0.5), 1e-4);

        // Basic tree: 1D:4.5

        tree.setInputValue("newick", "1D:4.5;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP,  -3.53223 + 2*Math.log(0.5), 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -3.907946 + 2*Math.log(0.5), 1e-4);


        // Basic tree: 1N:3;
        tree.setInputValue("newick", "1N:3;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -1.241888 + 2*Math.log(0.5), 1e-4);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -0.5752535 + 2*Math.log(0.5), 1e-4);



        // Basic tree: 1L:22;
        tree.setInputValue("newick", "1L:22;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.302585 + Math.log(0.5), 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.302585 + Math.log(0.5), 1e-4);

        // Small tree: (2D:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2D:4,3N:3)1D:2;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -6.401941 + 3*Math.log(0.5), 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -7.391145 + 3*Math.log(0.5), 1e-4);


        // Small tree: (2A:5,3L:6)1D:3
        tree.setInputValue("newick", "(2A:5,3L:6)1D:3;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -9.274919 + 2 * Math.log(0.5), 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -9.648058 + 2 * Math.log(0.5), 1e-4);

        // Bigger tree: ((4D:7,5N:2)2D:1,(6D:6,7A:9)3D:3)1D:4
        tree.setInputValue("newick", "((4D:7,5N:2)2D:1,(6D:6,7A:9)3D:3)1D:4;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -23.66497 + 5*Math.log(0.5), 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -25.11163 + 5*Math.log(0.5), 1e-4);

        // Bigger tree: ((4L:7,5A:2)2D:1,(6L:6,7N:9)3D:3)1D:4
        tree.setInputValue("newick", "((4L:7,5A:2)2D:1,(6L:6,7N:9)3D:3)1D:4;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -16.63671 + 3 * Math.log(0.5), 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -18.39805 + 3 * Math.log(0.5), 1e-4);
    }
}
