package test;

import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import junit.framework.TestCase;
import org.junit.Test;
import generalclassifier.core.LineageTree;
import generalclassifier.core.LineageTreeProb;
import generalclassifier.core.LineageTreeProbForLeafType;

public class LikelihoodTests_InferenceLeafType extends TestCase{

    @Test
    /**
     * Reference values come from calculations "by hand" in R
     */
    public void testLikelihoodWithLeavesOfInterest() throws Exception{
        LineageTree tree = new LineageTree();
        tree.setInputValue("isInferenceOfLeafType", new BooleanParameter("true"));
        LineageTreeProbForLeafType probTree = new LineageTreeProbForLeafType();
        probTree.setInputValue("fateProbsHSC", new RealParameter("0.4 0.35 0.25"));
        probTree.setInputValue("fateProbsMPP", new RealParameter("0.3 0.25 0.45"));
        probTree.setInputValue("scaleHSC", new RealParameter("2 2"));
        probTree.setInputValue("scaleMPP", new RealParameter("3 3"));
        probTree.setInputValue("shapeHSC", new RealParameter("1 1"));
        probTree.setInputValue("shapeMPP", new RealParameter("1 1"));
        probTree.setInputValue("lossProb", new RealParameter("0.1"));
        probTree.setInputValue("transitionProbs", new RealParameter("0.2 0.3 0.5"));

        double logP;


        // Basic tree: 1N:26;
        tree.setInputValue("newick", "1N:26;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("leafIsHSC", new BooleanParameter("true"));
        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("0"));
        probTree.setInputValue("leavesOfInterest", new IntegerParameter("1"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -1.491655, logP, 1e-4);

        probTree.setInputValue("leafIsHSC", new BooleanParameter("false"));
        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("0"));
        probTree.setInputValue("leavesOfInterest", new IntegerParameter("1"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -0.9038682, logP, 1e-4);


        // Tree with 2 leaves: (2L:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2L:4,3N:3)1D:2;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("leafIsHSC", new BooleanParameter("true"));
        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("0"));
        probTree.setInputValue("leavesOfInterest", new IntegerParameter("3"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -6.732182, logP, 1e-4);

        probTree.setInputValue("leafIsHSC", new BooleanParameter("false"));
        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("0"));
        probTree.setInputValue("leavesOfInterest", new IntegerParameter("3"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -5.74076, logP, 1e-4);

        // Tree with 2 leaves: (2N:4,3A:3)1D:2;
        tree.setInputValue("newick", "(2N:4,3A:3)1D:2;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("leafIsHSC", new BooleanParameter("true"));
        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("0"));
        probTree.setInputValue("leavesOfInterest", new IntegerParameter("2"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -7.861945, logP, 1e-4);

        probTree.setInputValue("leafIsHSC", new BooleanParameter("false"));
        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("0"));
        probTree.setInputValue("leavesOfInterest", new IntegerParameter("2"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -6.962136, logP, 1e-4);

        // Tree with 4 leaves: ((4A:7,5N:2)2D:1,(6D:6,7A:9)3D:3)1D:4;
        tree.setInputValue("newick", "((4A:7,5N:2)2D:1,(6D:6,7A:9)3D:3)1D:4;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("leafIsHSC", new BooleanParameter("false true false"));
        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("1")); // the offset of 1 is here to test that offsetting works
        probTree.setInputValue("leavesOfInterest", new IntegerParameter("5 6"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -26.89796, logP, 1e-4);

        probTree.setInputValue("leafIsHSC", new BooleanParameter("true false true true"));
        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("2")); // the offset of 2 is here to test that offsetting works
        probTree.setInputValue("leavesOfInterest", new IntegerParameter("5 6"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -27.65154, logP, 1e-4);

        // Tree with 4 leaves: ((4L:7,5N:2)2D:1,(6D:6,7L:9)3D:3)1D:4;
        tree.setInputValue("newick", "((4L:7,5N:2)2D:1,(6D:6,7L:9)3D:3)1D:4;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("leafIsHSC", new BooleanParameter("false true false"));
        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("1")); // the offset of 1 is here to test that offsetting works
        probTree.setInputValue("leavesOfInterest", new IntegerParameter("5 6"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -20.57085, logP, 1e-4);

        probTree.setInputValue("leafIsHSC", new BooleanParameter("true false true true"));
        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("2")); // the offset of 2 is here to test that offsetting works
        probTree.setInputValue("leavesOfInterest", new IntegerParameter("5 6"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -21.09545, logP, 1e-4);
    }

    public void testLikelihoodWithFixedRootType() throws Exception{
        LineageTree tree = new LineageTree();
        tree.setInputValue("isInferenceOfLeafType", new BooleanParameter("true"));
        LineageTreeProbForLeafType probTree = new LineageTreeProbForLeafType();
        probTree.setInputValue("fateProbsHSC", new RealParameter("0.4 0.4 0.2"));
        probTree.setInputValue("fateProbsMPP", new RealParameter("0.3 0.3 0.4"));
        probTree.setInputValue("scaleHSC", new RealParameter("2 2"));
        probTree.setInputValue("scaleMPP", new RealParameter("3 3"));
        probTree.setInputValue("shapeHSC", new RealParameter("1 1"));
        probTree.setInputValue("shapeMPP", new RealParameter("1 1"));
        probTree.setInputValue("lossProb", new RealParameter("0.1"));
        probTree.setInputValue("transitionProbs", new RealParameter("0.2 0.3 0.5"));

        double logP;

        // Basic tree: 1L:22;
        tree.setInputValue("newick", "1L:22;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(-2.302585, logP, 1e-4);

        // Basic tree: 1A:34;
        tree.setInputValue("newick", "1A:34;");
        tree.initAndValidate();
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(-18.7148, logP, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(-13.74128, logP, 1e-4);


        // Basic tree: 1N:26;
        tree.setInputValue("newick", "1N:26;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(-1.714798, logP, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(-1.021651, logP, 1e-4);


        // Tree with 2 leaves: (2L:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2L:4,3N:3)1D:2;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -6.144395, logP, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -6.398849, logP, 1e-4);

        // Tree with 2 leaves: (2A:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2A:4,3N:3)1D:2;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -7.566823, logP, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -7.837542, logP, 1e-4);

        // Tree with 4 leaves: ((4A:7,5A:2)2D:1,(6D:6,7A:9)3D:3)1D:4;
        tree.setInputValue("newick", "((4A:7,5A:2)2D:1,(6D:6,7A:9)3D:3)1D:4;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -26.78519, logP, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -27.52229, logP, 1e-4);

        // Tree with 4 leaves: ((4L:7,5A:2)2D:1,(6N:6,7L:9)3D:3)1D:4;
        tree.setInputValue("newick", "((4L:7,5A:2)2D:1,(6N:6,7L:9)3D:3)1D:4;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -17.50582, logP, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals( -18.59194, logP, 1e-4);
    }

    public void testLikelihoodWithNoLeafOfInterest() throws Exception{

        LineageTree tree = new LineageTree();
        tree.setInputValue("isInferenceOfLeafType", new BooleanParameter("true"));
        LineageTreeProbForLeafType probTree = new LineageTreeProbForLeafType();
        probTree.setInputValue("fateProbsHSC", new RealParameter("0.4 0.4 0.2"));
        probTree.setInputValue("fateProbsMPP", new RealParameter("0.3 0.3 0.4"));
        probTree.setInputValue("scaleHSC", new RealParameter("2 2"));
        probTree.setInputValue("scaleMPP", new RealParameter("3 3"));
        probTree.setInputValue("shapeHSC", new RealParameter("1 1"));
        probTree.setInputValue("shapeMPP", new RealParameter("1 1"));
        probTree.setInputValue("lossProb", new RealParameter("0.1"));
        probTree.setInputValue("transitionProbs", new RealParameter("0.2 0.3 0.5"));

        double logP;

        // Basic tree: 1L:22;
        tree.setInputValue("newick", "1L:22;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(-2.302585, logP, 1e-4);

        // Basic tree: 1A:34;
        tree.setInputValue("newick", "1A:34;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(-13.73439, logP, 1e-4);

        // Basic tree: 1N:26;
        tree.setInputValue("newick", "1N:26;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(-0.6161858, logP, 1e-4);

        // Tree with 2 leaves: (2L:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2L:4,3N:3)1D:2;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(-5.570403, logP, 1e-4);


        // Tree with 2 leaves: (2A:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2A:4,3N:3)1D:2;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(-6.999902, logP, 1e-4);


        // Tree with 4 leaves: ((4A:7,5A:2)2D:1,(6D:6,7A:9)3D:3)1D:4;
        tree.setInputValue("newick", "((4A:7,5A:2)2D:1,(6D:6,7A:9)3D:3)1D:4;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(-26.39416, logP, 1e-4);


        // Tree with 4 leaves: ((4L:7,5A:2)2D:1,(6N:6,7L:9)3D:3)1D:4;
        tree.setInputValue("newick", "((4L:7,5A:2)2D:1,(6N:6,7L:9)3D:3)1D:4;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(-17.2150, logP, 1e-4);
    }
}
