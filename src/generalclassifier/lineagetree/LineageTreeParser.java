package generalclassifier.lineagetree;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class LineageTreeParser {

    String inputFileName;
    String[] fluoChannelsCorrespondance;
    int maxNumberOfCells = Integer.MAX_VALUE;
    int maxTimePoint = Integer.MAX_VALUE;
    boolean isMaxTimeRelative = true;
    FrameRate frameRate = FrameRate.Hour;

    public LineageTreeParser(String inputFile) {
        this.inputFileName = inputFile;
    }

    public LineageTreeParser(String inputFile, String frameRate) {
        this.frameRate = FrameRate.fromString(frameRate);
        this.inputFileName = inputFile;
    }

    public void setMaxNumberOfCells(int maxNumberOfCells) {
        this.maxNumberOfCells = maxNumberOfCells;
    }

    public void setMaxTimePoint(int maxTimePoint) {
        this.maxTimePoint = maxTimePoint;
    }

    public void setIsMaxTimeRelative(boolean isRelative) {this.isMaxTimeRelative = isRelative;}

    //TODO deal with NA values
    //TODO implement check that tracknumber of root is 1, also force that children are always 2n, 2n+1
    public Map<Integer, Cell> parseRawCells() throws IOException {
        File csvData = new File(inputFileName);
        CSVParser parser = CSVParser.parse(new FileReader(csvData), CSVFormat.EXCEL.withHeader());

        Map<Integer, Cell> cells = new HashMap<>();

        ArrayList<MeasureType> measuresInCSV = new ArrayList<>();
        // add all measure that are actually in the CSV file to the list of measuresInCSV
        for(MeasureType measureType : MeasureType.values()) {
            if(measureType.isMeasureInCSV(parser))
                measuresInCSV.add(measureType);
        }

        boolean isFirstMeasure = true;

        // iterate over the lines of the csv files until an empty line is reached
        for (CSVRecord csvRecord : parser) {
            if(isEmptyRecord(csvRecord)) break;

            // get number of the cell in current line
            int cellNumber = getCellNumberInRecord(csvRecord);

            // get time of the datapoint in current line
            int time = getTimePointInRecord(csvRecord);

            if(isFirstMeasure) { // if it's the first measure and maxtimepoint was input relatively to the first time point, then change the maxtimepoint accordingly.
                if(isMaxTimeRelative) setMaxTimePoint(Math.max(time + maxTimePoint,maxTimePoint)); // use the max to make sure maxTimePoint does not become negative if already at Integer.MAX_VALUE
                isFirstMeasure = false;
            }

            if (time > maxTimePoint) //skip this data point if time is later than the maxTimePoint limit
                continue;

            if(cellNumber < 1)
                continue;

            int parentNumber = Cell.findParent(cellNumber);

            if (parentNumber > 0 && cells.containsKey(parentNumber))
                cells.get(parentNumber).addChild(cellNumber); // record the observed child of the parent, to know if this parent is a divider or not.

            if (cellNumber > maxNumberOfCells) // if cell is in a generation not taken into account, skip the data points
                continue;

            if(!cells.containsKey(cellNumber)) {
                // add a cell with the new track number
                cells.put(cellNumber, new Cell(cellNumber, measuresInCSV));
            }

            // assign the experimental measures in this line to the corresponding cell
            cells.get(cellNumber).addDataPoint(csvRecord);
        }
        return cells;
    }

    boolean isEmptyRecord(CSVRecord record) {
        return record.get("TrackNumber").isEmpty();
    }

    int getCellNumberInRecord(CSVRecord record) {
        return Integer.parseInt(record.get("TrackNumber"));
    }

    int getTimePointInRecord(CSVRecord record) {
        return Integer.parseInt(record.get("TimePoint"));
    }

    boolean isCellNumberInRecord(CSVRecord record, int cellNumber) {
        if(Integer.parseInt(record.get("TrackNumber")) == cellNumber) {
            return true;
        }
        return false;
    }

    public void setFluoChannelsCorrespondance(String[] fluoChannelsCorrespondance) {
        this.fluoChannelsCorrespondance = fluoChannelsCorrespondance;
    }

    public String[] getFluoChannelsCorrespondance() {
        return this.fluoChannelsCorrespondance;
    }

    enum FrameRate {
        Day("day"),
        Hour("hour"),
        Minute("minute"),
        Second("second");

        final String text;

        FrameRate(String s) {
            this.text = s;
        }

        public String getText() {
            return this.text;
        }

        // find the frameRate object that corresponds to the input string
        public static FrameRate fromString(String text) {
            for (FrameRate f : FrameRate.values()) {
                if (f.text.equalsIgnoreCase(text)) {
                    return f;
                }
            }
            throw new IllegalArgumentException("The frame rate input does not follow the mandatory format." +
                    "It should be one of these values: 'day', 'hour', 'minute', 'second'");
        }
    }

    public static void main(String[] parms) {
        LineageTreeParser parser =  new LineageTreeParser("../Data/Examples/testFile_shortLife.csv");

        try {
            Map<Integer, Cell> cells = parser.parseRawCells();
            System.out.println("Job done.");
        } catch (Exception e){
            System.out.println(e.getMessage());
        }
    }
}
