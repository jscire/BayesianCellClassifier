package BayesianCellClassifier.lineagetree;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.SortedSet;
import java.util.TreeMap;

public abstract class CellTree extends Tree {

    public Input<String> filePathInput = new Input<>("filePath", "File to write tree as a csv file in.");

    SortedSet<String> uniqueMeasurementTags;

    SortedSet<Integer> labelsOfAllCellsInTree;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if(filePathInput.get() != null) {
            writeToCSV(filePathInput.get());
        }
    }

    public void writeToCSV(String filePath){

        // build a sorted hashmap with the cell records, to sort them by track numbers
        // as getNodesAsArray will not necessarily spit the cells out in order.
        TreeMap<Integer, String> sortedCSVLines = new TreeMap<>();
        // add header
        sortedCSVLines.put(0, getHeaderOfCSV());
        // add cells
        for(Node cell : this.getNodesAsArray()){
            sortedCSVLines.put(((Cell) cell).getTrackNumber(), ((Cell) cell).toCSVRecord(uniqueMeasurementTags));
        }

        Charset utf8 = StandardCharsets.UTF_8;
        // write to file
        try {
            Files.write(Paths.get(filePath + ".csv"), sortedCSVLines.values(), utf8,
                    StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    abstract public SortedSet<Integer> getLabelsOfAllCellsInTree();

    public String getHeaderOfCSV(){
        String res="TrackNumber,Fate";
        for(String tag : uniqueMeasurementTags) {
            res += "," + tag;
        }
        return res + "";
    }

}
