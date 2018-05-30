package generalclassifier.core;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

import generalclassifier.utils.Utils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;


import java.io.*;
import java.util.*;

public class LineageTreeProb extends Distribution {

    public Input<LineageTree> lineageTreeInput = new Input<>("tree",
            "Lineage tree.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> fateProbsHSCInput = new Input<>("fateProbsHSC",
            "Vector containing probabilities of HSC being divider, apoptoser, non-divider or type-transitioner (if defined).",
            Input.Validate.REQUIRED);

    public Input<RealParameter> fateProbsMPPInput = new Input<>("fateProbsMPP",
            "Vector containing probabilities of MPP cell being divider, apoptoser or non-divider.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> scaleHSCInput = new Input<>("scaleHSC",
            "Vector containing Weibull scale parameters for HSC division, apoptosis and type-transition (if defined) times.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> scaleMPPInput = new Input<>("scaleMPP",
            "Vector containing Weibull scale parameters for MPP cell division and apoptosis times.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> shapeHSCInput = new Input<>("shapeHSC",
            "Vector containing Weibull shape parameters for HSC division, apoptosis and type-transition (if defined) times.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> shapeMPPInput = new Input<>("shapeMPP",
            "Vector containing Weibull shape parameters for MPP cell division and apoptosis times.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> lossProbInput = new Input<>("lossProb",
            "Probability of losing any given cell",
            Input.Validate.REQUIRED);

    public Input<BooleanParameter> rootIsHSCInput = new Input<>("rootIsHSC",
            "True if root is an HSC cell, false otherwise.",
            Input.Validate.REQUIRED);

    public Input<Integer> rootIsHSCidxInput = new Input<>("rootIsHSCidx",
            "If provided, this is an index into the rootIsHSC boolean " +
                    "parameter specifying which element corresponds to the root " +
                    "HSC status.");

    public Input<RealParameter> transitionUponDivisionProbsInput = new Input<>("transitionProbs",
            "Vector containing transition-upon-division probabilities (q0, q1, q2). These correspond to transitions upon division."); // to do rename input to be more explicit when this set-up works. For now I leave 'transitionProbs' to not mess up with what we already have that works.


    LineageTree tree;

    List<BitSet> typeConfigurationsHSC;
    BitSet typeConfigurationMPP;
    boolean rootIsHSC;

    boolean transitionUponDivisionIsAllowed;
    boolean transitionDuringLifetimeIsAllowed;

    /**
     * current and stored log probability of the potential configurations *
     */
    protected double[] configsLogProbsHSC;
    protected double[] storedConfigsLogProbsHSC;

    //TODO remove
//    final static double offsetFromBoundaries = 1E-9; // sets how far from the 0 and 1 boundaries should the numerical integration start and stop (to avoid numerical issues and because events cannot be simultaneous)
//
//
//    double lowerBound = offsetFromBoundaries;
//    double upperBound = 1-offsetFromBoundaries;
//
//    UniformRealDistribution uniformDistr = new UniformRealDistribution(lowerBound, upperBound);
//
//    RombergIntegrator numericalIntegrator = new RombergIntegrator();
      // TrapezoidIntegrator numericalIntegrator = new TrapezoidIntegrator(1e-6, 1e-6, 1, 60);
    IterativeLegendreGaussIntegrator numericalIntegrator = new IterativeLegendreGaussIntegrator(1, 1e-6, 1e-6);
    //  RombergIntegrator numericalIntegrator = new RombergIntegrator();

    int maxEval = 10000;


    public LineageTreeProb() {
    }

    @Override
    public void initAndValidate() {

        tree = lineageTreeInput.get();

        typeConfigurationsHSC = getTypeConfigurationsForSubtree(tree.getRoot(), true);
        configsLogProbsHSC = new double[typeConfigurationsHSC.size()];
        Arrays.fill(configsLogProbsHSC, Double.NaN);
        storedConfigsLogProbsHSC = new double[typeConfigurationsHSC.size()];
        Arrays.fill(storedConfigsLogProbsHSC, Double.NaN);
        typeConfigurationMPP = new BitSet();

        transitionUponDivisionIsAllowed = transitionUponDivisionProbsInput.get() != null;
        transitionDuringLifetimeIsAllowed = (fateProbsHSCInput.get().getDimension() == 4
                || shapeHSCInput.get().getDimension() == 3
                || scaleHSCInput.get().getDimension() == 3);

        if (!transitionUponDivisionIsAllowed && !transitionDuringLifetimeIsAllowed) {
            throw new IllegalStateException("Specify at least one type of transition (during lifetime or at division time).");
        }

        int formatFateProbHSC = fateProbsHSCInput.get().getDimension() == 4 ? 1 : 0;
        int formatShapeHSC = shapeHSCInput.get().getDimension() == 3 ? 1 : 0;
        int formatScaleMPP = scaleHSCInput.get().getDimension() == 3 ? 1 : 0;

        if ((formatFateProbHSC + formatScaleMPP + formatShapeHSC) == 1 ||
                (formatFateProbHSC + formatScaleMPP + formatShapeHSC) == 2) {
            throw new IllegalStateException("Faulty set-up. All 'transitionDuringLifetimeProbs' , 'fateProbsHSCInput' must be specified if one of the three is specified." +
                    "If that's the case, check that fateProbsHSC, shapeHSC and scaleHSC have the right dimensions (should be 4/3/3).");
        }

    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        if (rootIsHSCidxInput.get() != null)
            rootIsHSC = rootIsHSCInput.get().getValue(rootIsHSCidxInput.get());
        else
            rootIsHSC = rootIsHSCInput.get().getValue();


        if (rootIsHSC) {
            double maxLogProb = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < typeConfigurationsHSC.size(); i++) {
                try {
                    configsLogProbsHSC[i] = getProbabilityForConfiguration(typeConfigurationsHSC.get(i));
                } catch (Exception e) {
                    e.printStackTrace();
                }
                maxLogProb = Math.max(configsLogProbsHSC[i], maxLogProb);
            }

            // Compute sum of scaled configuration log probabilities
            double scaledP = 0.0;
            for (int i = 0; i < configsLogProbsHSC.length; i++)
                scaledP += Math.exp(configsLogProbsHSC[i] - maxLogProb);

            // Transform to log of unscaled marginal probability.
            logP = Math.log(scaledP) + maxLogProb;
        } else {
            try {
                logP = getProbabilityForConfiguration(typeConfigurationMPP);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

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

    //TODO deal with exceptions thrown by failed numerical integrations correctly
    public double getProbabilityForConfiguration(BitSet configuration) {
        double logConfProb = 0.0;


        for (Node node : tree.getNodesAsArray()) {
            int nodeIdx = node.getNr();

            boolean isHSC = configuration.get(nodeIdx);

            // Fate-dependent probabilities
            //TODO rearrange the code to simplify the case with transitions on branches
            //TODO allow for transition at one point or one integrated iver with an input setting
            switch (tree.getNodeFate(nodeIdx)) {
                case D:
                    if (isHSC) {

                        double noTransitionLogProb = Math.log(fateProbsHSCInput.get().getValue(0)) // cell stayed an HSC
                                + Utils.getWeibullLogDensity(tree.getEdgeLength(nodeIdx), // time until division
                                scaleHSCInput.get().getValue(0),
                                shapeHSCInput.get().getValue(0));

                        double q0 = 0;
                        double q1 = 0;
                        double q2 = 0;

                        if (transitionUponDivisionIsAllowed) {
                            q0 = transitionUponDivisionProbsInput.get().getValue(0);
                            q1 = transitionUponDivisionProbsInput.get().getValue(1);
                            q2 = transitionUponDivisionProbsInput.get().getValue(2);
                        }

                        double transitionLogProb = Double.NEGATIVE_INFINITY;

                        if (transitionDuringLifetimeIsAllowed) {
                            // version with transition at one point only
//                            transitionLogProb = Math.log(fateProbsHSCInput.get().getValue(3)) // HSC cell transitions
//                                    + Math.log(fateProbsMPPInput.get().getValue(0)) // MPP cell will divide
//                                    + Utils.getWeibullLogDensity(tree.getTimeToTransition(nodeIdx), scaleHSCInput.get().getValue(2), shapeHSCInput.get().getValue(2)) // time until transition
//                                    + Utils.getConditionalWeibullLogDensity(tree.getEdgeLength(nodeIdx), tree.getTimeToTransition(nodeIdx),
//                                    scaleMPPInput.get().getValue(0), shapeMPPInput.get().getValue(0)); // time until division

                            // version with integration over transition point
                            DensityProbOfStateTransition f = new DensityProbOfStateTransition(tree.getEdgeLength(nodeIdx),
                                    shapeHSCInput.get().getValue(2), scaleHSCInput.get().getValue(2),
                                    shapeMPPInput.get().getValue(0), scaleMPPInput.get().getValue(0));

                            try {
                                double integrationRes = numericalIntegrator.integrate(maxEval, f, 0, 1);
                                transitionLogProb = Math.log(integrationRes * fateProbsHSCInput.get().getValue(3) * fateProbsMPPInput.get().getValue(0));
                            } catch (Exception e) {
                                transitionLogProb = Double.NEGATIVE_INFINITY;
                            }

                        }

                        if (!node.isLeaf()) {
                            boolean child1HSC = configuration.get(node.getLeft().getNr());
                            boolean child2HSC = configuration.get(node.getRight().getNr());

                            if (child1HSC ^ child2HSC) {
                                if (transitionUponDivisionIsAllowed)
                                    logConfProb += Math.log(q1) + noTransitionLogProb;
                                else
                                    logConfProb = Double.NEGATIVE_INFINITY; // if different types and no transition upon division, probability of this config is zero

                            } else if (!child1HSC && !child2HSC) { // both children are MPP so there was a transition

                                if (!transitionUponDivisionIsAllowed && transitionDuringLifetimeIsAllowed)
                                    logConfProb += transitionLogProb;
                                else if (transitionUponDivisionIsAllowed && !transitionDuringLifetimeIsAllowed)
                                    logConfProb += Math.log(q0) + noTransitionLogProb;
                                else { // both types of transition are allowed
                                    logConfProb += Math.log(Math.exp(noTransitionLogProb + Math.log(q0)) + Math.exp(transitionLogProb));
                                }

                            } else { // both children are HSC so there was no transition
                                logConfProb += noTransitionLogProb;

                                if (transitionUponDivisionIsAllowed)
                                    logConfProb += Math.log(q2);
                            }

                        } else { // cell is leaf
                            if (transitionDuringLifetimeIsAllowed)
                                logConfProb += Math.log(Math.exp(noTransitionLogProb) + Math.exp(transitionLogProb));
                            else
                                logConfProb += noTransitionLogProb;
                        }

                    } else { // MPP cell
                        logConfProb += Math.log(fateProbsMPPInput.get().getValue(0))
                                + Utils.getWeibullLogDensity(tree.getEdgeLength(nodeIdx),
                                scaleMPPInput.get().getValue(0),
                                shapeMPPInput.get().getValue(0));
                    }

                    logConfProb += Math.log(1.0 - lossProbInput.get().getValue());
                    break;

                case A:
                    if (isHSC) {

                        double noTransitionLogProb = Math.log(fateProbsHSCInput.get().getValue(1)) // cell stayed an HSC
                                + Utils.getWeibullLogDensity(tree.getEdgeLength(nodeIdx),
                                scaleHSCInput.get().getValue(1),
                                shapeHSCInput.get().getValue(1));

                        double transitionProb;

                        if (transitionDuringLifetimeIsAllowed) {
                            //version with transition at one point
//                            transitionProb = Math.log(fateProbsHSCInput.get().getValue(3)) // HSC cell transitions
//                                    + Math.log(fateProbsMPPInput.get().getValue(1)) // MPP cell will apoptose
//                                    + Utils.getWeibullLogDensity(tree.getTimeToTransition(nodeIdx), scaleHSCInput.get().getValue(2), shapeHSCInput.get().getValue(2)) // time until transition
//                                    + Utils.getConditionalWeibullLogDensity(tree.getEdgeLength(nodeIdx), tree.getTimeToTransition(nodeIdx),
//                                    scaleMPPInput.get().getValue(1), shapeMPPInput.get().getValue(1)); // time until division

                            //logConfProb += Math.log(Math.exp(noTransitionLogProb) + Math.exp(transitionProb));

                            // version with integration over transition point
                            DensityProbOfStateTransition f = new DensityProbOfStateTransition(tree.getEdgeLength(nodeIdx),
                                    shapeHSCInput.get().getValue(2), scaleHSCInput.get().getValue(2),
                                    shapeMPPInput.get().getValue(1), scaleMPPInput.get().getValue(1));

                            try {
                                double integrationRes = numericalIntegrator.integrate(maxEval, f, 0, 1);
                                transitionProb = integrationRes * fateProbsHSCInput.get().getValue(3) * fateProbsMPPInput.get().getValue(1);
                            } catch (Exception e) {
                                transitionProb = 0;
                            }

                            logConfProb += Math.log(Math.exp(noTransitionLogProb) + transitionProb);
                        } else
                            logConfProb += noTransitionLogProb;

                    } else {
                        logConfProb += Math.log(fateProbsMPPInput.get().getValue(1))
                                + Utils.getWeibullLogDensity(tree.getEdgeLength(nodeIdx),
                                scaleMPPInput.get().getValue(1),
                                shapeMPPInput.get().getValue(1));
                    }

                    logConfProb += Math.log(1.0 - lossProbInput.get().getValue());
                    break;

                case N:
                    if (isHSC) {

                        double noTransitionProb = fateProbsHSCInput.get().getValue(2) // cell really does nothing
                                + fateProbsHSCInput.get().getValue(0) * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx) / scaleHSCInput.get().getValue(0), shapeHSCInput.get().getValue(0))) // cell would divide after the end of the branch
                                + fateProbsHSCInput.get().getValue(1) * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx) / scaleHSCInput.get().getValue(1), shapeHSCInput.get().getValue(1))); // cell would die after the end of the branch

                        double transitionProb;

                        if (transitionDuringLifetimeIsAllowed) {
                            // transition at one point
//                            transitionProb = fateProbsHSCInput.get().getValue(3) * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx)/scaleHSCInput.get().getValue(2), shapeHSCInput.get().getValue(2))) // HSC cell transitions after end of branch
//                                    + fateProbsHSCInput.get().getValue(3) * Utils.getWeibullDensity(tree.getTimeToTransition(nodeIdx), scaleHSCInput.get().getValue(2), shapeHSCInput.get().getValue(2)) // HSC cell transitions before end of branch
//                                    * (fateProbsMPPInput.get().getValue(2) //  MPP cell does nothing
//                                    + fateProbsMPPInput.get().getValue(0) * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx)/scaleMPPInput.get().getValue(0), shapeMPPInput.get().getValue(0))
//                                    + Math.pow(tree.getTimeToTransition(nodeIdx)/scaleMPPInput.get().getValue(0), shapeMPPInput.get().getValue(0))) // MPP cell divides after end of branch
//                                    + fateProbsMPPInput.get().getValue(1) * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx)/scaleMPPInput.get().getValue(1), shapeMPPInput.get().getValue(1))
//                                    + Math.pow(tree.getTimeToTransition(nodeIdx)/scaleMPPInput.get().getValue(1), shapeMPPInput.get().getValue(1)))); // MPP cell dies after end of branch

                            // integration over transition point
                            DensityProbTransitionBeforeEndAndDivisionDeathAfter f = new DensityProbTransitionBeforeEndAndDivisionDeathAfter(
                                    tree.getEdgeLength(nodeIdx),
                                    shapeHSCInput.get().getValue(2), scaleHSCInput.get().getValue(2),
                                    shapeMPPInput.get().getValue(0), scaleMPPInput.get().getValue(0),
                                    shapeMPPInput.get().getValue(1), scaleMPPInput.get().getValue(1),
                                    fateProbsMPPInput.get().getValue(0), fateProbsMPPInput.get().getValue(1), fateProbsMPPInput.get().getValue(2));


                            try {
                                double integrationRes = numericalIntegrator.integrate(maxEval, f, 0, 1);
                                transitionProb = fateProbsHSCInput.get().getValue(3)
                                        * (integrationRes // HSC transitions before end of branch, MPP cell does not do anything before end of branch
                                        + Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx) / scaleHSCInput.get().getValue(2), shapeHSCInput.get().getValue(2)))); // HSC cell transitions after end of branch
                            } catch (Exception e) {
                                transitionProb = 0;
                            }

                            logConfProb += Math.log(noTransitionProb + transitionProb);
                        } else
                            logConfProb += Math.log(noTransitionProb);

                    } else { // cell is MPP

                        //TODO remove (for debugging)
//                        double res = Math.log(fateProbsMPPInput.get().getValue(2) // cell really does nothing
//                                + fateProbsMPPInput.get().getValue(0) * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx)/scaleMPPInput.get().getValue(0), shapeMPPInput.get().getValue(0))) // cell would divide after the end of the branch
//                                + fateProbsMPPInput.get().getValue(1) * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx)/scaleMPPInput.get().getValue(1), shapeMPPInput.get().getValue(1))));

                        logConfProb += Math.log(fateProbsMPPInput.get().getValue(2) // cell really does nothing
                                + fateProbsMPPInput.get().getValue(0) * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx) / scaleMPPInput.get().getValue(0), shapeMPPInput.get().getValue(0))) // cell would divide after the end of the branch
                                + fateProbsMPPInput.get().getValue(1) * Math.exp(-Math.pow(tree.getEdgeLength(nodeIdx) / scaleMPPInput.get().getValue(1), shapeMPPInput.get().getValue(1)))); // cell would die after the end of the branch;
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

    public static List<BitSet> getTypeConfigurationsForSubtree(Node subtreeRoot, boolean rootIsHSC) {
        List<BitSet> configurations = new ArrayList<>();

        if (subtreeRoot.isLeaf()) {

            BitSet config = new BitSet();
            config.set(subtreeRoot.getNr(), rootIsHSC);
            configurations.add(config);

        } else {
            List<BitSet> leftConfigsTrue, leftConfigsFalse, rightConfigsTrue, rightConfigsFalse;
            leftConfigsFalse = getTypeConfigurationsForSubtree(subtreeRoot.getLeft(), false);
            rightConfigsFalse = getTypeConfigurationsForSubtree(subtreeRoot.getRight(), false);

            BitSet newConfigF = new BitSet();
            newConfigF.set(subtreeRoot.getNr(), rootIsHSC);
            configurations.add(newConfigF);


            if (rootIsHSC) {

                leftConfigsTrue = getTypeConfigurationsForSubtree(subtreeRoot.getLeft(), true);
                rightConfigsTrue = getTypeConfigurationsForSubtree(subtreeRoot.getRight(), true);

                for (BitSet configLeft : leftConfigsTrue) {
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

//    /**
//     * Compute log density at x under a Weibull distribution with scale lambda
//     * and shape parameter k.
//     *
//     * @param x parameter with distribution
//     * @param lambda scale parameter
//     * @param k shape parameter
//     * @return log density at x
//     */
//    public double getWeibullLogDensity(double x, double lambda, double k) {
//        double loglambda = Math.log(lambda);
//        double logk = Math.log(k);
//        double logx = Math.log(x);
//
//        return x<0.0
//                ? Double.NEGATIVE_INFINITY
//                : logk - loglambda + (k-1)*(logx - loglambda) - Math.pow(x/lambda, k);
//    }
//
//    /**
//     * Compute density at x under a Weibull distribution with scale lambda
//     * and shape parameter k.
//     *
//     * @param x parameter with distribution
//     * @param lambda scale parameter
//     * @param k shape parameter
//     * @return density at x
//     */
//    public double getWeibullDensity(double x, double lambda, double k) {
//
//        return x<0.0
//                ? 0.0
//                : k/lambda * Math.pow(x/lambda, k-1) * Math.exp(-Math.pow(x/lambda, k));
//    }
//
//    /**
//     * Compute density at x under a Weibull distribution with scale lambda
//     * and shape parameter k, conditioned on the lifetime x being higher than y.
//     *
//     * @param x parameter with distribution
//     * @param lambda scale parameter
//     * @param k shape parameter
//     * @return density at x
//     */
//    public double getConditionalWeibullLogDensity(double x, double y, double lambda, double k) {
//        double logOneMinusCDF = -Math.pow(y/lambda, k);
//
//        return (x<y || y<0)
//                ? Double.NEGATIVE_INFINITY
//               : getWeibullLogDensity(x, lambda, k) - logOneMinusCDF;
//    }
//
//
//
//    public double getWeibullLogCDF(double x, double lambda, double k) {
//
//        return (x<0)
//                ? Double.NEGATIVE_INFINITY
//                : Math.log(1 - Math.exp(-Math.exp(k * Math.log(x/lambda))));
//    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    public String getTreeID() {
        if (rootIsHSCidxInput.get() != null)
            return rootIsHSCidxInput.get() + "";
        else return "";
    }

    public boolean isHSCtree() {
        return rootIsHSC;
    }

    public List<BitSet> getTypeConfigurationsHSC() {
        return typeConfigurationsHSC;
    }

    public String getCorrespondanceLabelsToNodeArray() {

        String result = "";
        // for debugging
        int count = 0;

        for (Node node : tree.getNodesAsArray()) {
            int nodeNumber = node.getNr();
            if (nodeNumber != count) { // TODO remove ,for debugging
                System.out.println("Problem");
            }
            int newNodeNumber = Integer.parseInt(node.getID().substring(0, 1));
            count++;
            result += newNodeNumber;
        }
        return result;
    }

    /**
     * get configurations probabilities result from last known calculation, useful for logging
     *
     * @return log probability
     */
    public double[] getCurrentConfigsLogProbsHSC() {
        return configsLogProbsHSC;
    }


    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
        System.arraycopy(configsLogProbsHSC, 0, storedConfigsLogProbsHSC, 0, configsLogProbsHSC.length);
        super.store();
    }

    @Override
    public void restore() {
        System.arraycopy(storedConfigsLogProbsHSC, 0, configsLogProbsHSC, 0, storedConfigsLogProbsHSC.length);
        super.restore();
    }

    private class DensityProbOfStateTransition implements UnivariateFunction {

        double cellLifetime;

        double shapeLifetimeOriginState;
        double shapeLifetimeEndState;

        double scaleLifetimeOriginState;
        double scaleLifetimeEndState;

        DensityProbOfStateTransition(double cellLifetime, double shapeWeibullHSC, double scaleWeibullHSC, double shapeWeibullMPP, double scaleWeibullMPP) {
            this.cellLifetime = cellLifetime;

            this.shapeLifetimeOriginState = shapeWeibullHSC;
            this.scaleLifetimeOriginState = scaleWeibullHSC;

            this.shapeLifetimeEndState = shapeWeibullMPP;
            this.scaleLifetimeEndState = scaleWeibullMPP;
        }

        /**
         * @param x the fraction of the branch length for which the cell-state transition happens
         * @return the probability density value of the cell transitioning at proportion 'x' of the branch length
         */
        @Override
        public double value(double x) {

            double timeUntilTransition = x * cellLifetime;
            double beforeTransition = Utils.getWeibullDensity(timeUntilTransition, scaleLifetimeOriginState, shapeLifetimeOriginState);
            double afterTransition = Utils.getConditionalWeibullDensity(cellLifetime, timeUntilTransition, scaleLifetimeEndState, shapeLifetimeEndState);

            return beforeTransition * afterTransition;
        }
    }

    private class DensityProbBeforeTransition implements UnivariateFunction {

        double cellLifetime;

        double shapeLifetimeOriginState;

        double scaleLifetimeOriginState;

        DensityProbBeforeTransition(double cellLifetime, double shapeWeibull, double scaleWeibull) {
            this.cellLifetime = cellLifetime;

            this.shapeLifetimeOriginState = shapeWeibull;
            this.scaleLifetimeOriginState = scaleWeibull;
        }

        /**
         * @param x the fraction of the branch length for which the cell-state transition happens
         * @return the probability density value of the cell transitioning at proportion 'x' of the branch length
         */
        @Override
        public double value(double x) {

            double timeUntilTransition = x * cellLifetime;
            double beforeTransition = Utils.getWeibullDensity(timeUntilTransition, scaleLifetimeOriginState, shapeLifetimeOriginState);

            return beforeTransition;
        }
    }

    private class DensityProbTransitionBeforeEndAndDivisionDeathAfter implements UnivariateFunction {

        double cellLifetime;

        double shapeLifetimeOriginState;
        double shapeLifetimeEndState1;
        double shapeLifetimeEndState2;

        double scaleLifetimeOriginState;
        double scaleLifetimeEndState1;
        double scaleLifetimeEndState2;

        double probaEndState1;
        double probaEndState2;
        double probaEndState3;


        DensityProbTransitionBeforeEndAndDivisionDeathAfter(double cellLifetime,
                                                            double shapeWeibullTrans, double scaleWeibullTrans,
                                                            double shapeWeibullDiv, double scaleWeibullDiv,
                                                            double shapeWeibullDeath, double scaleWeibullDeath,
                                                            double probaDiv, double probaDeath, double probaNothing) {
            this.cellLifetime = cellLifetime;

            this.shapeLifetimeOriginState = shapeWeibullTrans;
            this.scaleLifetimeOriginState = scaleWeibullTrans;

            this.shapeLifetimeEndState1 = shapeWeibullDiv;
            this.scaleLifetimeEndState1 = scaleWeibullDiv;

            this.shapeLifetimeEndState2 = shapeWeibullDeath;
            this.scaleLifetimeEndState2 = scaleWeibullDeath;

            this.probaEndState1 = probaDiv;
            this.probaEndState2 = probaDeath;
            this.probaEndState3 = probaNothing;
        }

        /**
         * @param x the fraction of the branch length for which the cell-state transition happens
         * @return the probability density value of the cell transitioning at proportion 'x' of the branch length
         */
        @Override
        public double value(double x) {

            double timeUntilTransition = x * cellLifetime;
            double beforeTransition = Utils.getWeibullDensity(timeUntilTransition, scaleLifetimeOriginState, shapeLifetimeOriginState);
            double afterTransition1 = probaEndState1 * Utils.getConditionalWeibullProbaAboveT(cellLifetime, timeUntilTransition, scaleLifetimeEndState1, shapeLifetimeEndState1);
            double afterTransition2 = probaEndState2 * Utils.getConditionalWeibullProbaAboveT(cellLifetime, timeUntilTransition, scaleLifetimeEndState2, shapeLifetimeEndState2);
            double afterTransition3 = probaEndState3; // HSC transitions and MPP does nothing

            return beforeTransition * (afterTransition1 + afterTransition2 + afterTransition3);
        }
    }


    private static class DensityProbTransitionBeforeEndAndDivisionDeathAfterStatic implements UnivariateFunction {

        double cellLifetime;

        double shapeLifetimeOriginState;
        double shapeLifetimeEndState1;
        double shapeLifetimeEndState2;

        double scaleLifetimeOriginState;
        double scaleLifetimeEndState1;
        double scaleLifetimeEndState2;

        double probaEndState1;
        double probaEndState2;
        double probaEndState3;


        DensityProbTransitionBeforeEndAndDivisionDeathAfterStatic(double cellLifetime,
                                                                  double shapeWeibullTrans, double scaleWeibullTrans,
                                                                  double shapeWeibullDiv, double scaleWeibullDiv,
                                                                  double shapeWeibullDeath, double scaleWeibullDeath,
                                                                  double probaDiv, double probaDeath, double probaNothing) {
            this.cellLifetime = cellLifetime;

            this.shapeLifetimeOriginState = shapeWeibullTrans;
            this.scaleLifetimeOriginState = scaleWeibullTrans;

            this.shapeLifetimeEndState1 = shapeWeibullDiv;
            this.scaleLifetimeEndState1 = scaleWeibullDiv;

            this.shapeLifetimeEndState2 = shapeWeibullDeath;
            this.scaleLifetimeEndState2 = scaleWeibullDeath;

            this.probaEndState1 = probaDiv;
            this.probaEndState2 = probaDeath;
            this.probaEndState3 = probaNothing;
        }

        /**
         * @param x the fraction of the branch length for which the cell-state transition happens
         * @return the probability density value of the cell transitioning at proportion 'x' of the branch length
         */
        @Override
        public double value(double x) {

            double timeUntilTransition = x * cellLifetime;
            double beforeTransition = Utils.getWeibullDensity(timeUntilTransition, scaleLifetimeOriginState, shapeLifetimeOriginState);
            double afterTransition1 = probaEndState1 * Utils.getConditionalWeibullProbaAboveT(cellLifetime, timeUntilTransition, scaleLifetimeEndState1, shapeLifetimeEndState1);
            double afterTransition2 = probaEndState2 * Utils.getConditionalWeibullProbaAboveT(cellLifetime, timeUntilTransition, scaleLifetimeEndState2, shapeLifetimeEndState2);
            double afterTransition3 = probaEndState3; // HSC transitions and MPP does nothing

            return beforeTransition * (afterTransition1 + afterTransition2 + afterTransition3);
        }
    }

    public static void main(String[] args) throws IOException {

        //Tree tree = new TreeParser("((A:1,B:1):1,(C:1,D:1):1):1;");
        //Tree tree = new TreeParser("(A:1,B:1):1;");
//        Tree tree = new TreeParser("((A:1,B:1):1,C:2):1;");
//        System.out.println(tree);
//
//        List<BitSet> configs = getTypeConfigurationsForSubtree(tree.getRoot(), false);
//
//        for (BitSet config : configs) {
//            System.out.println(config);
//        }

        // test integrator

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
        tree.setInputValue("newick", "1N:3;");
        tree.setInputValue("transitionPosition", new RealParameter("0.35"));
        tree.initAndValidate();
        probTree.setInputValue("tree", tree);
        probTree.setInputValue("rootIsHSC", new BooleanParameter("true"));
        probTree.initAndValidate();
        logP = probTree.calculateLogP();

        int nodeIdx = 0;

        DensityProbTransitionBeforeEndAndDivisionDeathAfterStatic f = new DensityProbTransitionBeforeEndAndDivisionDeathAfterStatic(3,
                probTree.shapeHSCInput.get().getValue(2), probTree.scaleHSCInput.get().getValue(2),
                probTree.shapeMPPInput.get().getValue(0), probTree.scaleMPPInput.get().getValue(0),
                probTree.shapeMPPInput.get().getValue(1), probTree.scaleMPPInput.get().getValue(1),
                probTree.fateProbsMPPInput.get().getValue(0), probTree.fateProbsMPPInput.get().getValue(1), probTree.fateProbsMPPInput.get().getValue(2));


        IterativeLegendreGaussIntegrator numericalIntegrator = new IterativeLegendreGaussIntegrator(1, 1e-5, 1e-5);
//
        // RombergIntegrator numericalIntegrator = new RombergIntegrator();
//        double integrationRes = numericalIntegrator.integrate(10000 ,f, 0, 1);
//        System.out.println(integrationRes);
//
//        int numberOfSteps = 10001;
//        double[] xValues = new double[numberOfSteps];
//        double[] yValues = new double[numberOfSteps];
//
//
//        double stepSize = 1.0/(numberOfSteps-1);
//        System.out.println(stepSize);
//
//        xValues[0] = 0;
//        for (int i = 1; i < xValues.length; i++) {
//            xValues[i] = xValues[i-1] + stepSize;
////            System.out.println(xValues[i]);
//        }
//
//        StringBuilder builder = new StringBuilder();
//        DecimalFormat df = new DecimalFormat("0.00E0##");
//
//        for (int i = 0; i < xValues.length; i++) {
//            yValues[i] = f.value(xValues[i]);
//            builder.append(df.format(xValues[i]) + "\t" + df.format(yValues[i]) + "\n");
//        }
//
//        try (BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(
//                new FileOutputStream("valuesDensityFunction.txt"), "utf-8"))) {
//            for (int i = 0; i < xValues.length; i++) {
//                yValues[i] = f.value(xValues[i]);
//                writer.write(df.format(xValues[i]) + "\t" + df.format(yValues[i]) + "\n");
//            }
//        }


    }


}



