package generalclassifier.utils;

import java.io.*;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import org.apache.commons.math3.distribution.*;

public class SimulatedLineageTree {

    public static final int maxNumberOfGenerations = 3;
    static final double maxBranchLength = 1E6;
    static DecimalFormat df = new DecimalFormat("#########0");
    static Random randomGenerator = new Random();

    public enum Fate {
        D, A, N, L
    }

    private Fate[] nodeFates;

    private double edgeLengths[];

    public boolean rootIsHSC;

    public static double[] shapeHSC = new double[]{15, 28};
    public static double[] scaleHSC = new double[]{5, 11};
    public static double[] shapeMPP = new double[]{10, 22};
    public static double[] scaleMPP = new double[]{3, 7};

    public static double[] fateProbabilitiesHSC = new double[]{0.7, 0.1, 0.2};
    public static double[] fateProbabilitiesMPP = new double[]{0.4, 0.35, 0.25};

    public static double lostProbability = 0.1;

    static WeibullDistribution divisionHSC = new WeibullDistribution(shapeHSC[0], scaleHSC[0]);
    static WeibullDistribution divisionMPP = new WeibullDistribution(shapeMPP[0], scaleMPP[0]);
    static WeibullDistribution apoptoseHSC = new WeibullDistribution(shapeHSC[1], scaleHSC[1]);
    static WeibullDistribution apoptoseMPP = new WeibullDistribution(shapeMPP[1], scaleMPP[1]);

    public static double[] transitionProbabilities = new double[]{0.1, 0.2, 0.7};

    NodeLineageTree root;
    int nextNodeID;

    public SimulatedLineageTree(boolean treeIsHSC){
        rootIsHSC = treeIsHSC;
        nextNodeID = 1;

        root = new NodeLineageTree(1, rootIsHSC);
    }

    /**
     * set all the parameters that define the simulation of the tree.
     * The only tricky thing is that the input array must be in the right order.
     * The required order is:
     * fateProbabilitiesHSC 1,2,3
     * fateProbabilitiesMPP 1,2,3
     * scaleHSC 1,2
     * scaleMPP 1,2
     * shapeHSC 1,2
     * shapeMPP 1,2
     * lostProbability
     * transitionProbabilities 1,2,3
     *
     * @param setOfParameters
     */
    public static void setAllParameters(double[] setOfParameters){

        if(setOfParameters.length != 18) {
            throw new RuntimeException("Wrong set of input for setAllParameters.");
        }

        fateProbabilitiesHSC[0] = setOfParameters[0];
        fateProbabilitiesHSC[1] = setOfParameters[1];
        fateProbabilitiesHSC[2] = setOfParameters[2];

        fateProbabilitiesMPP[0] = setOfParameters[3];
        fateProbabilitiesMPP[1] = setOfParameters[4];
        fateProbabilitiesMPP[2] = setOfParameters[5];

        scaleHSC[0] = setOfParameters[6];
        scaleHSC[1] = setOfParameters[7];

        scaleMPP[0] = setOfParameters[8];
        scaleMPP[1] = setOfParameters[9];

        shapeHSC[0] = setOfParameters[10];
        shapeHSC[1] = setOfParameters[11];

        shapeMPP[0] = setOfParameters[12];
        shapeMPP[1] = setOfParameters[13];

        lostProbability = setOfParameters[14];

        transitionProbabilities[0] = setOfParameters[15];
        transitionProbabilities[1] = setOfParameters[16];
        transitionProbabilities[2] = setOfParameters[17];


        divisionHSC = new WeibullDistribution(shapeHSC[0], scaleHSC[0]);
        divisionMPP = new WeibullDistribution(shapeMPP[0], scaleMPP[0]);
        apoptoseHSC = new WeibullDistribution(shapeHSC[1], scaleHSC[1]);
        apoptoseMPP = new WeibullDistribution(shapeMPP[1], scaleMPP[1]);
    }

    public String toNewickString(){
        return root.subTreeToNewick() + ";";
    }

    public String toString(){
        return root.subTreeToString() + ";";
    }

    public String getTypesOfLeavesOfInterest(int treeNumber){
        String formattedResult = "";
        ArrayList<String> leavesOfInterest = root.getTypesOfLeavesOfInterestInSubtree();
        for (String leaf : leavesOfInterest) {
            formattedResult += treeNumber + "\t" + leaf + "\n";
        }
        return formattedResult;
    }

    class NodeLineageTree {

        int nodeID;
        Fate cellFate;
        double edgeLength;
        boolean isHSCCell;
        NodeLineageTree leftChild;
        NodeLineageTree rightChild;


        NodeLineageTree(int generationNumber, boolean isAnHSC){
            nodeID = nextNodeID;
            nextNodeID ++;
            isHSCCell = isAnHSC;
            this.setRandomCellFate();
            this.setRandomEdgeLength();

            if(generationNumber >= maxNumberOfGenerations || cellFate == Fate.L || cellFate == Fate.A || cellFate == Fate.N)
                return;
            else {
                if(!isHSCCell) { // MPP parent cell
                    leftChild = generateSubLineageTree(generationNumber + 1, false);
                    rightChild = generateSubLineageTree(generationNumber + 1, false);
                } else { // HSC parent cell
                    // allow for transitions from HSC to MPP during division
                    double randomDivisionEvents = randomGenerator.nextDouble();
                    if (randomDivisionEvents < transitionProbabilities[0]) {// HSC gives 2 MPPS
                        leftChild = generateSubLineageTree(generationNumber + 1, false);
                        rightChild = generateSubLineageTree(generationNumber + 1, false);

                    } else if(randomDivisionEvents < (transitionProbabilities[0] + transitionProbabilities[1])) { // HSC gives one HSC and one MPP
                        leftChild = generateSubLineageTree(generationNumber + 1, true);
                        rightChild = generateSubLineageTree(generationNumber + 1, false);

                    } else { // HSC gives 2 HSC
                        leftChild = generateSubLineageTree(generationNumber + 1, true);
                        rightChild = generateSubLineageTree(generationNumber + 1, true);
                    }
                }
            }
        }

        NodeLineageTree generateSubLineageTree(int currentGenerationNumber, boolean isHSC){
            NodeLineageTree childCell = new NodeLineageTree(currentGenerationNumber, isHSC);
            return childCell;
        }

        void setRandomCellFate() {
            double isLostCell = randomGenerator.nextDouble();
            if(isLostCell < lostProbability) {
                cellFate = Fate.L;
                return;
            }

            // Otherwise, if cell was not lost, sample its fate following the fate probabilities of its type
            double sampledFate = randomGenerator.nextDouble();
            double[] fateProbabilities = isHSCCell? fateProbabilitiesHSC: fateProbabilitiesMPP;

            if(sampledFate < fateProbabilities[0]) {
                cellFate = Fate.D;
            } else if (sampledFate < (fateProbabilities[0] + fateProbabilities[1])) {
                cellFate = Fate.A;
            } else {
                cellFate = Fate.N;
            }
        }

        void setRandomEdgeLength(){
            switch(cellFate){
                case D:
                    edgeLength = isHSCCell? divisionHSC.sample(): divisionMPP.sample();
                    break;
                case A:
                    edgeLength = isHSCCell? apoptoseHSC.sample(): apoptoseMPP.sample();
                    break;

                case N:
                case L:
                    edgeLength = randomGenerator.nextDouble() * maxBranchLength; // draw a random branch length in this case, as we are not interested in the branch length there
                    break;
            }
        }

        int getSubTreeSize(){
            if(leftChild == null)
                return 1;
            else return 1 + leftChild.getSubTreeSize() + rightChild.getSubTreeSize();
        }


        String subTreeToNewick(){

            if(leftChild == null) { // node is a leaf
                return this.nodeToNewick();
            }
            String newickString = "(" + leftChild.subTreeToNewick() + "," + rightChild.subTreeToNewick() + ")" + this.nodeToNewick();
            return newickString;
        }

        String nodeToNewick(){
            return "" + nodeID + cellFate + ":" + df.format(edgeLength);
        }

        String subTreeToString(){ // return a newick format with the type of every cell included.

            if(leftChild == null) { // node is a leaf
                return this.nodeToString();
            }
            String newickFull = "(" + leftChild.subTreeToString() + "," + rightChild.subTreeToString() + ")" + this.nodeToString();
            return newickFull;
        }

        String nodeToString(){
            String type = isHSCCell? "H": "M";
            return "" + nodeID + cellFate + type + ":" + df.format(edgeLength);
        }

        boolean isLeafOfInterest(){
            if(leftChild != null || cellFate == Fate.A || cellFate == Fate.L) return false; // not a leaf or not of interest
            return true;
        }

        ArrayList<String> getTypesOfLeavesOfInterestInSubtree(){
            ArrayList<String> result = new ArrayList<>();
            if(this.isLeafOfInterest())
                result.add(nodeID + "\t" + isHSCCell);
            if(leftChild == null) return result; // leaf but not of interest
            else{ // not leaf
                result.addAll(leftChild.getTypesOfLeavesOfInterestInSubtree());
                result.addAll(rightChild.getTypesOfLeavesOfInterestInSubtree());
            }
            return result;
        }
    }

    public static void main(String[] args) throws IOException {




        /// GENERATES 200 HSC trees with transitions and writes the newick strings in a file
        // also lof the leaves of interest
        int numberOfHSCtrees = 200;
        SimulatedLineageTree[] HSCtrees = new SimulatedLineageTree[numberOfHSCtrees];
        String[] HSCtreesAsNewick = new String[numberOfHSCtrees];
        String[] leavesOfInterest = new String[numberOfHSCtrees];

        for (int i = 0; i < numberOfHSCtrees; i++) {
            HSCtrees[i] = new SimulatedLineageTree(true);
            HSCtreesAsNewick[i] = HSCtrees[i].toNewickString();
            leavesOfInterest[i] = HSCtrees[i].getTypesOfLeavesOfInterest(i);
        }

        String outputFileHSC = "SimulatedTreesHSC.txt";

        List<String> linesHSC = Arrays.asList(HSCtreesAsNewick);
        Path fileHSC = Paths.get(outputFileHSC);
        Files.write(fileHSC, linesHSC, Charset.forName("UTF-8"));

        String outputFileLeavesOfInterest = "SimulatedTreesHSC_leavesOfInterest.txt";

        List<String> linesLeaves = Arrays.asList(leavesOfInterest);
        Path fileLeavesHSC = Paths.get(outputFileLeavesOfInterest);
        Files.write(fileLeavesHSC, linesLeaves, Charset.forName("UTF-8"));

        /// GENERATES 200 trees of random type with transitions and writes the full 'newick' strings in a file
        int numberOftrees = 200;
        SimulatedLineageTree[] trees = new SimulatedLineageTree[numberOftrees];
        String[] treesAsFullNewick = new String[numberOftrees];

        for (int i = 0; i < numberOftrees; i++) {
            trees[i] = new SimulatedLineageTree(randomGenerator.nextBoolean());
            treesAsFullNewick[i] = trees[i].toString();
        }

        String outputFile = "SimulatedTreesWithTransitions.txt";

        List<String> lines = Arrays.asList(treesAsFullNewick);
        Path fileTrees = Paths.get(outputFile);
        Files.write(fileTrees, lines, Charset.forName("UTF-8"));


        // GENERATES 200 pure MPP trees (all cells are MPP) and writes the newick strings in a file
        int numberOfMPPtrees = 200;
        SimulatedLineageTree[] MPPtrees = new SimulatedLineageTree[numberOfMPPtrees];
        String[] MPPtreesAsNewick = new String[numberOfMPPtrees];

        for (int i = 0; i < numberOfMPPtrees; i++) {
            MPPtrees[i] = new SimulatedLineageTree(false);
            MPPtreesAsNewick[i] = MPPtrees[i].toNewickString();
        }

        String outputFileMPP = "SimulatedTreesMPP.txt";

        List<String> linesMPP = Arrays.asList(MPPtreesAsNewick);
        Path fileMPP = Paths.get(outputFileMPP);
        Files.write(fileMPP, linesMPP, Charset.forName("UTF-8"));


        /// COMMENTED OUT, GENERATES 200 pure HSC trees (all cells are HSC) and writes the newick strings in a file

//        int numberOfHSCtrees = 200;
//        SimulatedLineageTree[] pureHSCtrees = new SimulatedLineageTree[numberOfHSCtrees];
//        String[] pureHSCtreesAsNewick = new String[numberOfHSCtrees];
//
//        for (int i = 0; i < numberOfHSCtrees; i++) {
//            pureHSCtrees[i] = new SimulatedLineageTree(true);
//            pureHSCtreesAsNewick[i] = pureHSCtrees[i].toNewickString();
//        }
//
//        String outputFile = "SimulatedTreesPureHSC.txt";
//
//        List<String> lines = Arrays.asList(pureHSCtreesAsNewick);
//        Path file = Paths.get(outputFile);
//        Files.write(file, lines, Charset.forName("UTF-8"));


        // GENERATES 200 pure MPP trees (all cells are MPP) and writes the newick strings in a file
//        int numberOfMPPtrees = 200;
//        SimulatedLineageTree[] pureMPPtrees = new SimulatedLineageTree[numberOfMPPtrees];
//        String[] pureMPPtreesAsNewick = new String[numberOfMPPtrees];
//
//        for (int i = 0; i < numberOfMPPtrees; i++) {
//            pureMPPtrees[i] = new SimulatedLineageTree(false);
//            pureMPPtreesAsNewick[i] = pureMPPtrees[i].toNewickString();
//        }
//
//        String outputFileMPP = "SimulatedTreesPureMPP.txt";
//
//        List<String> linesMPP = Arrays.asList(pureMPPtreesAsNewick);
//        Path fileMPP = Paths.get(outputFileMPP);
//        Files.write(fileMPP, linesMPP, Charset.forName("UTF-8"));
    }
}


