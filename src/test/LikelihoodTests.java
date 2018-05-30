package test;

import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import junit.framework.TestCase;
import org.junit.Test;
import scclassify.LineageTree;
import scclassify.LineageTreeProb;

public class LikelihoodTests extends TestCase{

    @Test
    /**
     * Reference values come from calculations "by hand" in R
     */
    public void testBasicLikelihood() throws Exception{
        LineageTree tree = new LineageTree();
        LineageTreeProb probTree = new LineageTreeProb();
        probTree.setInputValue("fateProbsHSC", new RealParameter("0.4 0.4 0.2"));
        probTree.setInputValue("fateProbsMPP", new RealParameter("0.3 0.3 0.4"));
        probTree.setInputValue("scaleHSC", new RealParameter("2 2"));
        probTree.setInputValue("scaleMPP", new RealParameter("3 3"));
        probTree.setInputValue("shapeHSC", new RealParameter("1 1"));
        probTree.setInputValue("shapeMPP", new RealParameter("1 1"));
        probTree.setInputValue("lossProb", new RealParameter("0.1"));
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.setInputValue("transitionProbs", new RealParameter("0.2 0.3 0.5"));

        double logP;

        // Basic tree: 1A:34;
        tree.setInputValue("newick", "1A:34;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -18.7148, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -13.74128, 1e-4);


        // Basic tree: 1N:26;
        tree.setInputValue("newick", "1N:26;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -1.714789, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -1.021393, 1e-4);


        // Basic tree: 1L:22;
        tree.setInputValue("newick", "1L:22;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.302585, 1e-4);

        // Tree with 2 leaves: (2L:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2L:4,3N:3)1D:2;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -5.611862, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -5.959421, 1e-4);

        // Tree with 2 leaves: (2A:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2A:4,3N:3)1D:2;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -7.034257, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -7.398115, 1e-4);

        // Tree with 4 leaves: ((4A:7,5A:2)2D:1,(6D:6,7A:9)3D:3)1D:4;
        tree.setInputValue("newick", "((4A:7,5A:2)2D:1,(6D:6,7A:9)3D:3)1D:4;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -26.78519, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -27.52229, 1e-4);

        // Tree with 4 leaves: ((4L:7,5A:2)2D:1,(6N:6,7L:9)3D:3)1D:4;
        tree.setInputValue("newick", "((4L:7,5A:2)2D:1,(6N:6,7L:9)3D:3)1D:4;");
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -17.32192, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -18.40712, 1e-4);
    }

    @Test
    /**
     * Reference values come from calculations "by hand" in R
     */
    public void testLikelihoodWithStateTransitionOnBranchesOnly() throws Exception{
        LineageTree tree = new LineageTree();
        LineageTreeProb probTree = new LineageTreeProb();
        probTree.setInputValue("fateProbsHSC", new RealParameter("0.6 0.2 0.1 0.1"));
        probTree.setInputValue("fateProbsMPP", new RealParameter("0.3 0.3 0.4"));
        probTree.setInputValue("scaleHSC", new RealParameter("2 2.2 2.5"));
        probTree.setInputValue("scaleMPP", new RealParameter("3 3.1"));
        probTree.setInputValue("shapeHSC", new RealParameter("1 1.5 1.3"));
        probTree.setInputValue("shapeMPP", new RealParameter("1 1.2"));
        probTree.setInputValue("lossProb", new RealParameter("0.1"));
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));

        double logP;

        // Bigger tree: ((4D:7,5N:2)2D:1,(6D:6,7A:9)3D:3)1D:4
        tree.setInputValue("newick", "((4D:7,5N:2)2D:1,(6D:6,7A:9)3D:3)1D:4;");
        tree.setInputValue("transitionPosition", new RealParameter("0.3 0.5 0.7 0.3 0.32 0.8 0.5"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -25.94084, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -25.11163, 1e-4);

        // Bigger tree: ((4L:7,5A:2)2D:1,(6L:6,7N:9)3D:3)1D:4
        tree.setInputValue("newick", "((4L:7,5A:2)2D:1,(6L:6,7N:9)3D:3)1D:4;");
        tree.setInputValue("transitionPosition", new RealParameter("0.3 0.5 0.7 0.3 0.32 0.8 0.5"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -17.68962, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -18.39805, 1e-4);


        // Small tree: (2D:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2D:4,3N:3)1D:2;");
        tree.setInputValue("transitionPosition", new RealParameter("0.4 0.55 0.71"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP,  -6.809311, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -7.391145, 1e-4);


        // Small tree: (2A:5,3L:6)1D:3
        tree.setInputValue("newick", "(2A:5,3L:6)1D:3;");
        tree.setInputValue("transitionPosition", new RealParameter("0.89 0.6 0.71"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -10.0361, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -9.648058, 1e-4);


        // Basic tree: 1A:2;
        tree.setInputValue("newick", "1A:2;");
        tree.setInputValue("transitionPosition", new RealParameter("0.4"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.974327, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.937085, 1e-4);


        // Basic tree: 1N:26;
        tree.setInputValue("newick", "1N:26;");
        tree.setInputValue("transitionPosition", new RealParameter("0.35"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.406493, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -1.02152, 1e-4);

        // Basic tree: 1D:12

        tree.setInputValue("newick", "1D:12;");
        tree.setInputValue("transitionPosition", new RealParameter("0.45"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -7.247093, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -6.407946, 1e-4);

        // Basic tree: 1L:22;
        tree.setInputValue("newick", "1L:22;");
        tree.setInputValue("transitionPosition", new RealParameter("0.2"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.302585, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.302585, 1e-4);




    }

    @Test
    /**
     * Reference values come from calculations "by hand" in R
     */
    public void testLikelihoodWithStateTransitionOnBranchesAndUponDivision() throws Exception{

        LineageTree tree = new LineageTree();
        LineageTreeProb probTree = new LineageTreeProb();
        probTree.setInputValue("fateProbsHSC", new RealParameter("0.6 0.2 0.1 0.1"));
        probTree.setInputValue("fateProbsMPP", new RealParameter("0.3 0.3 0.4"));
        probTree.setInputValue("scaleHSC", new RealParameter("2 2.2 2.5"));
        probTree.setInputValue("scaleMPP", new RealParameter("3 3.1"));
        probTree.setInputValue("shapeHSC", new RealParameter("1 1.5 1.3"));
        probTree.setInputValue("shapeMPP", new RealParameter("1 1.2"));
        probTree.setInputValue("lossProb", new RealParameter("0.1"));
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.setInputValue("transitionProbs", new RealParameter("0.2 0.3 0.5"));

        double logP;

        // Basic tree: 1A:2;
        tree.setInputValue("newick", "1A:2;");
        tree.setInputValue("transitionPosition", new RealParameter("0.4"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.974327, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.937085, 1e-4);


        // Basic tree: 1N:3;
        tree.setInputValue("newick", "1N:3;");
        tree.setInputValue("transitionPosition", new RealParameter("0.35"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -1.23475, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -0.5752535, 1e-4);

        // Basic tree: 1D:4.5

        tree.setInputValue("newick", "1D:4.5;");
        tree.setInputValue("transitionPosition", new RealParameter("0.45"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP,  -3.528195, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -3.907946, 1e-4);

        // Basic tree: 1L:22;
        tree.setInputValue("newick", "1L:22;");
        tree.setInputValue("transitionPosition", new RealParameter("0.2"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.302585, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.302585, 1e-4);

        // Small tree: (2D:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2D:4,3N:3)1D:2;");
        tree.setInputValue("transitionPosition", new RealParameter("0.4 0.55 0.71"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -6.398898, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -7.391145, 1e-4);


        // Small tree: (2A:5,3L:6)1D:3
        tree.setInputValue("newick", "(2A:5,3L:6)1D:3;");
        tree.setInputValue("transitionPosition", new RealParameter("0.89 0.6 0.71"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -9.268804, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -9.648058, 1e-4);

        // Bigger tree: ((4D:7,5N:2)2D:1,(6D:6,7A:9)3D:3)1D:4
        tree.setInputValue("newick", "((4D:7,5N:2)2D:1,(6D:6,7A:9)3D:3)1D:4;");
        tree.setInputValue("transitionPosition", new RealParameter("0.3 0.5 0.7 0.3 0.32 0.8 0.5"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -23.65285, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -25.11163, 1e-4);

        // Bigger tree: ((4L:7,5A:2)2D:1,(6L:6,7N:9)3D:3)1D:4
        tree.setInputValue("newick", "((4L:7,5A:2)2D:1,(6L:6,7N:9)3D:3)1D:4;");
        tree.setInputValue("transitionPosition", new RealParameter("0.3 0.5 0.7 0.3 0.32 0.8 0.5"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -16.6362, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -18.39805, 1e-4);


    }

    @Test
    /**
     * Reference values come from calculations "by hand" in R
     */
    public void testLikelihoodWithIntegrationOverTransitionTime() throws Exception{

        LineageTree tree = new LineageTree();
        LineageTreeProb probTree = new LineageTreeProb();
        probTree.setInputValue("fateProbsHSC", new RealParameter("0.6 0.2 0.1 0.1"));
        probTree.setInputValue("fateProbsMPP", new RealParameter("0.3 0.3 0.4"));
        probTree.setInputValue("scaleHSC", new RealParameter("2 2.2 2.5"));
        probTree.setInputValue("scaleMPP", new RealParameter("3 3.1"));
        probTree.setInputValue("shapeHSC", new RealParameter("1 1.5 1.3"));
        probTree.setInputValue("shapeMPP", new RealParameter("1 1.2"));
        probTree.setInputValue("lossProb", new RealParameter("0.1"));
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.setInputValue("transitionProbs", new RealParameter("0.2 0.3 0.5"));

        double logP;

        // Basic tree: 1A:2;
        tree.setInputValue("newick", "1A:2;");
        tree.setInputValue("transitionPosition", new RealParameter("0.4"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.975158, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.937085, 1e-4);

        // Basic tree: 1D:4.5

        tree.setInputValue("newick", "1D:4.5;");
        tree.setInputValue("transitionPosition", new RealParameter("0.45"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP,  -3.53223, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -3.907946, 1e-4);


        // Basic tree: 1N:3;
        tree.setInputValue("newick", "1N:3;");
        tree.setInputValue("transitionPosition", new RealParameter("0.35"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -1.241888, 1e-4);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -0.5752535, 1e-4);



        // Basic tree: 1L:22;
        tree.setInputValue("newick", "1L:22;");
        tree.setInputValue("transitionPosition", new RealParameter("0.2"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.302585, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -2.302585, 1e-4);

        // Small tree: (2D:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2D:4,3N:3)1D:2;");
        tree.setInputValue("transitionPosition", new RealParameter("0.4 0.55 0.71"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -6.401941, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -7.391145, 1e-4);


        // Small tree: (2A:5,3L:6)1D:3
        tree.setInputValue("newick", "(2A:5,3L:6)1D:3;");
        tree.setInputValue("transitionPosition", new RealParameter("0.89 0.6 0.71"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -9.274919, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -9.648058, 1e-4);

        // Bigger tree: ((4D:7,5N:2)2D:1,(6D:6,7A:9)3D:3)1D:4
        tree.setInputValue("newick", "((4D:7,5N:2)2D:1,(6D:6,7A:9)3D:3)1D:4;");
        tree.setInputValue("transitionPosition", new RealParameter("0.3 0.5 0.7 0.3 0.32 0.8 0.5"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -23.66497, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -25.11163, 1e-4);

        // Bigger tree: ((4L:7,5A:2)2D:1,(6L:6,7N:9)3D:3)1D:4
        tree.setInputValue("newick", "((4L:7,5A:2)2D:1,(6L:6,7N:9)3D:3)1D:4;");
        tree.setInputValue("transitionPosition", new RealParameter("0.3 0.5 0.7 0.3 0.32 0.8 0.5"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -16.63671, 1e-4);

        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        assertEquals(logP, -18.39805, 1e-4);


    }
}
