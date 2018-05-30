package generalclassifier.core;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import generalclassifier.utils.MiscUtils;

import java.util.*;

//TODO have this class extend LineageTreeProb to avoid duplicated code
public class LineageTreeProbForLeafType extends Distribution {

    public Input<LineageTree> lineageTreeInput = new Input<>("tree",
            "Lineage tree.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> fateProbsHSCInput = new Input<>("fateProbsHSC",
            "Vector containing probabilities of HSC being divider, apoptoser or non-divider.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> fateProbsMPPInput = new Input<>("fateProbsMPP",
            "Vector containing probabilities of MPP cell being divider, apoptoser or non-divider.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> scaleHSCInput = new Input<>("scaleHSC",
            "Vector containing Weibull scale parameters for HSC division and apoptosis times.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> scaleMPPInput = new Input<>("scaleMPP",
            "Vector containing Weibull scale parameters for MPP cell division and apoptosis times.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> shapeHSCInput = new Input<>("shapeHSC",
            "Vector containing Weibull shape parameters for HSC division and apoptosis times.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> shapeMPPInput = new Input<>("shapeMPP",
            "Vector containing Weibull shape parameters for MPP cell division and apoptosis times.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> lossProbInput = new Input<>("lossProb",
            "Probability of losing any given cell",
            Input.Validate.REQUIRED);

    public Input<BooleanParameter> leafIsHSCInput = new Input<>("leafIsHSC",
            "True if leaf is an HSC cell, false otherwise.");

    public Input<IntegerParameter> leafIsHSCidxInput = new Input<>("leafIsHSCidx",
            "If provided, this is an index into the leafIsHSC boolean " +
                    "parameter specifying which element corresponds to the first leaf of interest " +
                    "HSC status.");


    public Input<IntegerParameter> leavesOfInterestInput = new Input<>("leavesOfInterest",
            "If provided, this is the indices of the leaves whose type is being inferred " +
                    "HSC status.");

    public Input<BooleanParameter> rootIsHSCInput = new Input<>("rootIsHSC",
            "True if root is an HSC cell, false otherwise.");


    public Input<RealParameter> transitionProbsInput = new Input<>("transitionProbs",
            "Vector containing transition probabilities (q0, q1, q2).",
            Input.Validate.REQUIRED);

    LineageTree tree;

    List<BitSet> typeConfigurations;
    double[] configLogProbs;

    int[] leavesOfInterest = new int[]{-1}; // by default, there is no leaf of interest

    public LineageTreeProbForLeafType() { }

    @Override
    public void initAndValidate() {

        tree = lineageTreeInput.get();

        if(leafIsHSCidxInput.get() != null) {
            leavesOfInterest = new int[leavesOfInterestInput.get().getDimension()];

            for (int i = 0; i < leavesOfInterest.length; i++) {
                leavesOfInterest[i] = leavesOfInterestInput.get().getValue(i) -1; // the "-1" is there to account for the fact that cell indices start at 1 in the datasets
            }
        }
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        //TODO can this be rationalized a bit better?
        if (rootIsHSCInput.get() != null) { // root type is known
            typeConfigurations = getTypeConfigurations(tree, rootIsHSCInput.get().getValue());
        } else if(tree.getNodeCount()==1 && isLeafOfInterest(tree.getRoot())) { // if only one cell in the tree and is the cell of interest
            typeConfigurations = getTypeConfigurations(tree, isHSC(tree.getRoot()));
        } else {
            typeConfigurations = getTypeConfigurations(tree);
        }

        configLogProbs = new double[typeConfigurations.size()];

        double maxLogProb = Double.NEGATIVE_INFINITY;
        for (int i=0; i<typeConfigurations.size(); i++) {
            configLogProbs[i] = getProbabilityForConfiguration(typeConfigurations.get(i));
            maxLogProb = Math.max(configLogProbs[i], maxLogProb);
        }

        // Compute sum of scaled configuration probabilities
        double scaledP = 0.0;
        for (int i = 0; i< configLogProbs.length; i++)
            scaledP += Math.exp(configLogProbs[i]-maxLogProb);

        // Transform to log of unscaled marginal probability.
        logP = Math.log(scaledP) + maxLogProb;

        return logP;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        throw new RuntimeException("Not implemented.");
    }

    double getProbabilityForConfiguration(BitSet configuration) {
        double logConfProb = 0.0;

        //for (int nodeIdx = 0; nodeIdx<tree.getNodeCount(); nodeIdx++) { // for debugging // TODO remove if below doesn't work
        for (Node node: tree.getNodesAsArray()) {
            int nodeIdx = node.getNr();

            boolean isHSC = configuration.get(nodeIdx);

            // Fate-dependent probabilities
            switch(tree.getNodeFate(nodeIdx)) { // debugging
                case D:
                    if (isHSC) {
                        logConfProb += Math.log(fateProbsHSCInput.get().getValue(0))
                                + getWeibullLogDensity(tree.getEdgeLength(nodeIdx),
                                scaleHSCInput.get().getValue(0),
                                shapeHSCInput.get().getValue(0));

                        if (!node.isLeaf()) {
                            double q0 = transitionProbsInput.get().getValue(0);
                            double q1 = transitionProbsInput.get().getValue(1);
                            double q2 = transitionProbsInput.get().getValue(2);

                            boolean child1HSC = configuration.get(node.getLeft().getNr());
                            boolean child2HSC = configuration.get(node.getRight().getNr());

                            if (child1HSC && child2HSC) {
                                logConfProb += Math.log(q2);
                            } else {
                                if (child1HSC || child2HSC) {
                                    logConfProb += Math.log(q1);
                                } else {
                                    logConfProb += Math.log(q0);
                                }
                            }
                        }
                    } else {
                        logConfProb += Math.log(fateProbsMPPInput.get().getValue(0))
                                + getWeibullLogDensity(tree.getEdgeLength(nodeIdx),
                                scaleMPPInput.get().getValue(0),
                                shapeMPPInput.get().getValue(0));
                    }

                    logConfProb += Math.log(1.0 - lossProbInput.get().getValue());
                    break;

                case A:
                    if (isHSC) {
                        logConfProb += Math.log(fateProbsHSCInput.get().getValue(1))
                                + getWeibullLogDensity(tree.getEdgeLength(nodeIdx),
                                scaleHSCInput.get().getValue(1),
                                shapeHSCInput.get().getValue(1));
                    } else {
                        logConfProb += Math.log(fateProbsMPPInput.get().getValue(1))
                                + getWeibullLogDensity(tree.getEdgeLength(nodeIdx),
                                scaleMPPInput.get().getValue(1),
                                shapeMPPInput.get().getValue(1));
                    }

                    logConfProb += Math.log(1.0 - lossProbInput.get().getValue());
                    break;

                case N:
                    if (isHSC) {
                        logConfProb += Math.log(fateProbsHSCInput.get().getValue(2));
                    } else {
                        logConfProb += Math.log(fateProbsMPPInput.get().getValue(2));
                    }
                    logConfProb += Math.log(1.0 - lossProbInput.get().getValue());
                    break;

                case L:
                    logConfProb += Math.log(lossProbInput.get().getValue());
                    break;

                default:
                    throw new IllegalStateException("Invalid cell fate.");
            }
        }

        return logConfProb;
    }

    public List<BitSet> getTypeConfigurations (LineageTree tree){

        List<BitSet> configurations = getTypeConfigurationsForSubtree(tree.getRoot(), true);

        // in case the tree only has one leaf and it was leaf,
        // we don't add the configuration with root is MPP as the type of the root is indifferent
        if (tree.getNodeCount() == 1 && tree.getNodeFate(tree.getRoot().getNr()) == LineageTree.Fate.L){
            return configurations;
        }

        List<BitSet> configurationsFalse = getTypeConfigurationsForSubtree(tree.getRoot(), false);
        configurations.addAll(configurationsFalse); // there is only one configuration with root is MPP
        return configurations;
    }

    public List<BitSet> getTypeConfigurations (Tree tree, boolean rootIsHSC){ // in case the root is of a know type (MPP or HSC)
        if (rootIsHSC) return getTypeConfigurationsForSubtree(tree.getRoot(), rootIsHSC);

        else{
            List<BitSet> configurations = new ArrayList<BitSet>();
            configurations.add(new BitSet());
            return configurations;}
    }

    public List<BitSet> getTypeConfigurationsForSubtree (Node subtreeRoot, boolean nodeIsHSC) {
        List<BitSet> configurations = new ArrayList<>();

        if (subtreeRoot.isLeaf()) {
            if(!isLeafOfInterest(subtreeRoot) || isHSC(subtreeRoot) == nodeIsHSC) { // if leaf is not of set type (it's not a leaf of interest) or if it's the right type
                BitSet config = new BitSet();
                config.set(subtreeRoot.getNr(), nodeIsHSC);
                configurations.add(config);
            }

        } else {
            List<BitSet> leftConfigsTrue, leftConfigsFalse, rightConfigsTrue, rightConfigsFalse;
            leftConfigsFalse = getTypeConfigurationsForSubtree(subtreeRoot.getLeft(), false);
            rightConfigsFalse = getTypeConfigurationsForSubtree(subtreeRoot.getRight(), false);

            if(!leftConfigsFalse.isEmpty() && !rightConfigsFalse.isEmpty()) { // if configurations are possible with both subtrees' root nodes being MPP, add this configuration
                BitSet newConfig = new BitSet();
                newConfig.set(subtreeRoot.getNr(), nodeIsHSC);
                configurations.add(newConfig);
            }

            if (nodeIsHSC) {

                leftConfigsTrue = getTypeConfigurationsForSubtree(subtreeRoot.getLeft(), true);
                rightConfigsTrue = getTypeConfigurationsForSubtree(subtreeRoot.getRight(), true);

                for (BitSet configLeft : leftConfigsTrue) { // the nested for loop should eliminate all configurations for which there is no possible configuration in the subtree
                    for (BitSet configRight : rightConfigsTrue) {
                        BitSet newConfig = new BitSet();
                        newConfig.or(configLeft);
                        newConfig.or(configRight);
                        newConfig.set(subtreeRoot.getNr(), true);
                        configurations.add(newConfig);
                    }
                }

                for (BitSet configLeft : leftConfigsTrue) {
                    for (BitSet configRight : rightConfigsFalse) {
                        BitSet newConfig = new BitSet();
                        newConfig.or(configLeft);
                        newConfig.or(configRight);
                        newConfig.set(subtreeRoot.getNr(), true);
                        configurations.add(newConfig);
                    }
                }

                for (BitSet configLeft : leftConfigsFalse) {
                    for (BitSet configRight : rightConfigsTrue) {
                        BitSet newConfig = new BitSet();
                        newConfig.or(configLeft);
                        newConfig.or(configRight);
                        newConfig.set(subtreeRoot.getNr(), true);
                        configurations.add(newConfig);
                    }
                }
            }
        }

        return configurations;
    }

    boolean isHSC(Node currentLeaf) {
        int numLeaf = currentLeaf.getNr(); // for debugging

        int idxHSCLeaf = leafIsHSCidxInput.get().getValue() + MiscUtils.match(leavesOfInterest, currentLeaf.getNr());
        return leafIsHSCInput.get().getValue(idxHSCLeaf);
    }

    boolean isLeafOfInterest(Node leaf){
        return MiscUtils.contains(leavesOfInterest, leaf.getNr());
    }

    /**
     * Compute density at x under a Weibull distribution with mean lambda
     * and shape parameter k.
     *
     * @param x parameter with distribution
     * @param lambda scale parameter
     * @param k shape parameter
     * @return density at x
     */
    public double getWeibullLogDensity(double x, double lambda, double k) {
        double loglambda = Math.log(lambda);
        double logk = Math.log(k);
        double logx = Math.log(x);

        return x<0.0
                ? Double.NEGATIVE_INFINITY
                : logk - loglambda + (k-1)*(logx - loglambda) - Math.pow(x/lambda, k);
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    public static void main(String[] args) {


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
//        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));

        probTree.setInputValue("transitionProbs", new RealParameter("0.2 0.3 0.5"));

        double logP;

        // Tree with 3 leaves: (2A:2,3D:1(4N:1.5, 5A:2))1D:1.5;
//        tree.setInputValue("newick", "(2A:2,(4N:1.5, 5A:2)3D:1)1D:1.5;");
//        tree.initAndValidate();
//        probTree.setInputValue("tree", tree);
//        probTree.initAndValidate();
//        logP = probTree.calculateLogP();
//        System.out.println(logP + " versus -11.20425");


        // Basic tree: 1A:34;
//        tree.setInputValue("newick", "1A:34;");
//        tree.initAndValidate();
//        probTree.setInputValue("leafIsHSC", new BooleanParameter("true"));
//        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("0"));
//        probTree.setInputValue("leavesOfInterest", new IntegerParameter("1"));
//        probTree.setInputValue("tree", tree);
//        probTree.initAndValidate();
//        logP = probTree.calculateLogP();
//        System.out.println(logP + " versus -18.7148");

        // Tree with 2 leaves: (2L:4,3N:3)1D:2;
        tree.setInputValue("newick", "(2L:4,3N:3)1D:2;");
        tree.initAndValidate();
        probTree.setInputValue("leafIsHSC", new BooleanParameter("true"));
        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("0"));
        probTree.setInputValue("leavesOfInterest", new IntegerParameter("3"));
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();
        System.out.println(logP);

//        assertEquals(logP, -6.955325, 1e-4);

        probTree.setInputValue("leafIsHSC", new BooleanParameter("false"));
        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("0"));
        probTree.setInputValue("leavesOfInterest", new IntegerParameter("3"));
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();
        System.out.println(logP + " versus -5.858543");

        // Tree with 3 leaves: (2A:2,3D:1(4N:1.5, 5A:2))1D:1.5;
        tree.setInputValue("newick", "(2A:2,(4N:1.5, 5A:2)3D:1)1D:1.5;");
        tree.initAndValidate();
        probTree.setInputValue("leafIsHSC", new BooleanParameter("true"));
        probTree.setInputValue("leafIsHSCidx", new IntegerParameter("0"));
        probTree.setInputValue("leavesOfInterest", new IntegerParameter("4"));
        probTree.setInputValue("tree", tree);
        probTree.initAndValidate();
        logP = probTree.calculateLogP();
        System.out.println(logP + " versus -12.51083");





        //Tree tree = new TreeParser("((A:1,B:1):1,(C:1,D:1):1):1;");
        //Tree tree = new TreeParser("(A:1,B:1):1;");
//        Tree tree = new TreeParser("((A:1,B:1):1,C:2):1;");
//        Tree tree = new TreeParser("((A:1,B:1):1,(C:2,D:3):2):1;");
//        System.out.println(tree);
//
//        LineageTreeProbForLeafType prob = new LineageTreeProbForLeafType();
//        List<BitSet> configs = prob.getTypeConfigurations(tree, true);
//
//        int count = 0;
//        for (BitSet config : configs) {
//            System.out.println(config);
//            count ++;
//        }
//        System.out.println(count);
    }
}
