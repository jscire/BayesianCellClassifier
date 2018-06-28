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
    FrameRate frameRate = FrameRate.Hour;

    public LineageTreeParser(String inputFile) {
        this.inputFileName = inputFile;
    }

    public LineageTreeParser(String inputFile, String frameRate) {
        this.frameRate = FrameRate.fromString(frameRate);
        this.inputFileName = inputFile;
    }

    public Map<Integer, Cell> parse() throws IOException {
        File csvData = new File(inputFileName);
//        CSVFormat format = CSVFormat.newFormat(',').withHeader();
        CSVParser parser = CSVParser.parse(new FileReader(csvData), CSVFormat.EXCEL.withHeader());

        int maxCellNumber = 0;
        Map<Integer, Cell> cells = new HashMap<>();


        ArrayList<ExperimentalMeasure> measuresInCSV = new ArrayList<>();
        // add all measure that are actually in the CSV file to the list of measuresInCSV
        for(MeasureType measureType : MeasureType.values()) {
            if(measureType.isMeasureInCSV(parser))
                measuresInCSV.add(new ExperimentalMeasure(measureType));
        }


        for (CSVRecord csvRecord : parser) {

            if(isEmptyRecord(csvRecord)) break;

            // get number of cell in current line
            int cellNumber = getCellNumberInRecord(csvRecord);

            if(!cells.containsKey(cellNumber)) {
                if(maxCellNumber < cellNumber) maxCellNumber = cellNumber;
                // add a cell with the new track number
                cells.put(cellNumber, new Cell(cellNumber, measuresInCSV));
            }

            // assign the experimental measures in this line to the corresponding cell
            cells.get(cellNumber).addTimePoint(csvRecord);
        }

        return cells;
    }

    boolean isEmptyRecord(CSVRecord record) {
        return record.get("TrackNumber").isEmpty();
    }

    int getCellNumberInRecord(CSVRecord record) {
        return Integer.parseInt(record.get("TrackNumber"));
    }


    boolean isCellNumberInRecord(CSVRecord record, int cellNumber) {
        if(Integer.parseInt(record.get("TrackNumber")) == cellNumber) {
            return true;
        }
        return false;
    }

    // not used
    // TODO remove if useless
    public Cell parseCell(int CellNumber) throws IOException {
        File csvData = new File(inputFileName);
        CSVParser parser = CSVParser.parse(new FileReader(csvData), CSVFormat.EXCEL);

        Cell cell = new Cell(CellNumber);

        for (CSVRecord csvRecord : parser) {

            if(isCellNumberInRecord(csvRecord, cell.trackNumber))
                // if current line (record) contains this cellNumber, add a timePoint to the stored observations
                cell.addTimePoint(csvRecord);

            else if(isCellNumberInRecord(csvRecord, 2 * cell.trackNumber)
                    || isCellNumberInRecord(csvRecord, 2 * cell.trackNumber + 1))
                // if one of the daughters of the cell is in the record, stop searching lower, the mother cell does not exist anymore
                break;
        }

        return cell;
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
        LineageTreeParser parser =  new LineageTreeParser("../Data/Examples/toyFile.csv");

        try {
            Map<Integer, Cell> cells = parser.parse();
            System.out.println("Job done.");
        } catch (Exception e){
            System.out.println(e.getMessage());
        }
    }
}
