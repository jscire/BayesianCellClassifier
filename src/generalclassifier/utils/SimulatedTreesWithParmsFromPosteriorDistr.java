package generalclassifier.utils;

import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.Scanner;
import java.io.File;


public class SimulatedTreesWithParmsFromPosteriorDistr {

    static Random randomGenerator = new Random();

    String[] treesAsNewick;
    String[] leavesOfInterest;
    String[] rootIsHSC;
    boolean isMixedDataset;

    SimulatedTreesWithParmsFromPosteriorDistr(boolean isMixedDataset, int numberOfTrees) {
        rootIsHSC = new String[numberOfTrees];
        treesAsNewick = new String[numberOfTrees];
        leavesOfInterest = new String[numberOfTrees];
        this.isMixedDataset = isMixedDataset;
    }


    static public SimulatedTreesWithParmsFromPosteriorDistr simulateSetOfTrees(String inputFile, int paramsLineNumber, boolean isMixedDataset, int numberOfTrees) throws Exception {

        SimulatedTreesWithParmsFromPosteriorDistr simulatedTrees = new SimulatedTreesWithParmsFromPosteriorDistr(isMixedDataset, numberOfTrees);
        SimulatedLineageTree tempTree;

        int lineNumber = 0;
        for (Scanner sc = new Scanner(new File(inputFile)); sc.hasNext(); ) {

            String line = sc.nextLine();
            String[] stringSetOfParms = line.split(" ");

            double[] sampledParms = new double[stringSetOfParms.length];
            for (int i = 0; i < sampledParms.length; i++) {
                sampledParms[i] = Double.parseDouble(stringSetOfParms[i]);
                //   System.out.print(stringSetOfParms[i] + " "); // for debugging
            }

            if(lineNumber == paramsLineNumber){
                SimulatedLineageTree.setAllParameters(sampledParms); // assign the set of sampled parameters to the tree simulator
                break;
            }
            lineNumber ++;
        }

        if(lineNumber < paramsLineNumber)
            throw new RuntimeException("Wrong paramsLineNumber, inputFile is too short.");

        for (int i = 0; i < numberOfTrees; i++) {
            boolean nextTreeIsHSC = isMixedDataset? randomGenerator.nextBoolean(): false; // choose type of root at random if it's a mixed dataset, otherwise only MPPs
            tempTree = new SimulatedLineageTree(nextTreeIsHSC); // simulate the tree

            simulatedTrees.treesAsNewick[i] = tempTree.toNewickString(); // store the tree as a string
            simulatedTrees.rootIsHSC[i] = String.valueOf(nextTreeIsHSC); // store the type of the root independently
            simulatedTrees.leavesOfInterest[i] = tempTree.getTypesOfLeavesOfInterest(i);
        }

        return simulatedTrees;
    }



    static public SimulatedTreesWithParmsFromPosteriorDistr simulateSetOfTreesAndChangeParmsForEachTree(String inputFile, boolean isMixedDataset, int numberOfTrees) throws Exception {

        SimulatedTreesWithParmsFromPosteriorDistr simulatedTrees = new SimulatedTreesWithParmsFromPosteriorDistr(isMixedDataset, numberOfTrees);

        int countTrees = 0;
        SimulatedLineageTree tempTree;

        for(Scanner sc = new Scanner(new File(inputFile)); sc.hasNext(); ) {

            if(countTrees == numberOfTrees) break; // break the loop if more sets of parameters were drawn than trees that are to be simulated

            String line = sc.nextLine();
            String[] stringSetOfParms = line.split(" ");

            double[] sampledParms = new double[stringSetOfParms.length];
            for (int i = 0; i < sampledParms.length; i++) {
                sampledParms[i] = Double.parseDouble(stringSetOfParms[i]);
                //   System.out.print(stringSetOfParms[i] + " "); // for debugging
            }
            SimulatedLineageTree.setAllParameters(sampledParms); // assign the set of sampled parameters to the tree simulator

            boolean nextTreeIsHSC = isMixedDataset? randomGenerator.nextBoolean(): false; // choose type of root at random if it's a mixed dataset, otherwise only MPPs
            tempTree = new SimulatedLineageTree(nextTreeIsHSC); // simulate the tree

            simulatedTrees.treesAsNewick[countTrees] = tempTree.toNewickString(); // store the tree as a string
            simulatedTrees.rootIsHSC[countTrees] = String.valueOf(nextTreeIsHSC); // store the type of the root independently

            countTrees ++;
            // System.out.println(""); // for debugging
        }

        return simulatedTrees;
    }

    void writeLeavesOfInterestToFile(String outputFileTrees) throws Exception {
        List<String> linesLeaves = Arrays.asList(leavesOfInterest);
        Path fileLeavesHSC = Paths.get(outputFileTrees + "_leavesOfInterest_types.txt");
        Files.write(fileLeavesHSC, linesLeaves, Charset.forName("UTF-8"));
    }

    void writeTreeSetToFile(String outputFileTrees) throws Exception {
        List<String> linesTrees = Arrays.asList(this.treesAsNewick);
        Path fileTrees = Paths.get(outputFileTrees + ".trees");
        Files.write(fileTrees, linesTrees, Charset.forName("UTF-8"));

        if (isMixedDataset) {
            List<String> linesRootTypes = Arrays.asList(this.rootIsHSC);
            Path fileRootType = Paths.get(outputFileTrees + "_RootTypes.txt");
            Files.write(fileRootType, linesRootTypes, Charset.forName("UTF-8"));
        }
    }

    public static void main(String[] args) throws Exception {

        int numOfMixedTrees = 266;
        int numOfMPPTrees = 267;

        for (int lineNumber = 0; lineNumber < 30; lineNumber++) {
            SimulatedTreesWithParmsFromPosteriorDistr simulatedProperMixedTrees = simulateSetOfTrees("samplesFromPosteriorDistr.txt", lineNumber, true, numOfMixedTrees);
            simulatedProperMixedTrees.writeTreeSetToFile("Simulated_LeafTypeInference_MixedTrees_" + lineNumber);
            simulatedProperMixedTrees.writeLeavesOfInterestToFile("Simulated_LeafTypeInference_MixedTrees_" + lineNumber);

            SimulatedTreesWithParmsFromPosteriorDistr simulatedProperMPPTrees = simulateSetOfTrees("samplesFromPosteriorDistr.txt",lineNumber, false, numOfMPPTrees);
            simulatedProperMPPTrees.writeTreeSetToFile("Simulated_LeafTypeInference_MPPTrees_" + lineNumber);
        }


        //// simulation to check if we can infer the root types properly
//        for (int lineNumber = 0; lineNumber < 100; lineNumber++) {
//            SimulatedTreesWithParmsFromPosteriorDistr simulatedProperMixedTrees = simulateSetOfTrees("samplesFromPosteriorDistr.txt", lineNumber, true, numOfMixedTrees);
//            simulatedProperMixedTrees.writeTreeSetToFile("SimulatedProperMixedTrees_" + lineNumber);
//
//            SimulatedTreesWithParmsFromPosteriorDistr simulatedProperMPPTrees = simulateSetOfTrees("samplesFromPosteriorDistr.txt",lineNumber, false, numOfMPPTrees);
//            simulatedProperMPPTrees.writeTreeSetToFile("SimulatedProperMPPTrees_" + lineNumber);
//        }

//        SimulatedTreesWithParmsFromPosteriorDistr simulatedMixedTrees = simulateSetOfTreesAndChangeParmsForEachTree("samplesFromPosteriorDistr.txt", true, numOfMixedTrees);
//        simulatedMixedTrees.writeTreeSetToFile("SimulatedTreesMixDiffParms");
//
//        SimulatedTreesWithParmsFromPosteriorDistr simulatedMPPTrees = simulateSetOfTreesAndChangeParmsForEachTree("samplesFromPosteriorDistr.txt", false, numOfMPPTrees);
//        simulatedMPPTrees.writeTreeSetToFile("SimulatedTreesMPPDiffParms");



    }
}
