package generalclassifier.lineagetree;

import org.apache.commons.csv.CSVParser;

public enum MeasureType {

    ElapsedTime(new String[]{"time_sec", "time_hours"}),
    Lifetime(new String[]{"TimePoint"}),
    Area(new String[]{"Area", "AreaMorphologyCh00"}, 0, Double.POSITIVE_INFINITY, false),
    Eccentricity(new String[]{"Eccentricity", "EccentricityMorphologyCh00"}, 0, 1, false),
    Perimeter(new String[]{"Perimeter", "PerimeterMorphologyCh00"}, 0 , Double.POSITIVE_INFINITY, false),
    XPosition(new String[]{"CentroidX","XMorphologyCh00"}),
    YPosition(new String[]{"CentroidY","YMorphologyCh00"}),
    InstantSpeed(new String[]{"Speed"}, 0, Double.POSITIVE_INFINITY, false),

    // fluorescent markers
    Ig2afc(new String[]{"Ig2afc", "Ig2afc_signal"}, 0, Double.POSITIVE_INFINITY),
    Sca1meanValue(new String[]{"Sca1", "Sca1_signal"}, 0, Double.POSITIVE_INFINITY),
    Sca1meanRate(new String[]{"Sca1", "Sca1_signal"}, 0, Double.POSITIVE_INFINITY),
    CD41(new String[]{"CD41", "CD41_signal"}, 0, Double.POSITIVE_INFINITY),
    FcgRIII(new String[]{"FcgRIII", "FcgRIII_signal"}, 0, Double.POSITIVE_INFINITY),
    cMycGFP(new String[]{"cMycGFP", "cMycGFP_signal"}, 0, Double.POSITIVE_INFINITY),
    ROS(new String[]{"ROS", "ROS_signal"}, 0, Double.POSITIVE_INFINITY),
    TMRMmean(new String[]{"TMRM", "TMRM_signal"}, 0, Double.POSITIVE_INFINITY),
    TMRMmax(new String[]{"TMRM", "TMRM_signal"}, 0, Double.POSITIVE_INFINITY),
    CD71APC(new String[]{"CD71APC", "CD71APC_signal"}, 0, Double.POSITIVE_INFINITY),
    CD71PE(new String[]{"CD71PE", "CD71PE_signal"}, 0, Double.POSITIVE_INFINITY);

    private String[] namesInInputFile;

    private boolean isAccurateForRootCell = true;

//    The hard lower and upper bounds refer to the acceptable values for the experimental measures that
    private double hardLowerBound = Double.NEGATIVE_INFINITY;
    private double hardUpperBound = Double.POSITIVE_INFINITY;

    MeasureType(String[] namesInCSV) {
        this.namesInInputFile = namesInCSV;
    }

    MeasureType(String[] namesInCSV, double lower, double upper) {
        this.namesInInputFile = namesInCSV;
        this.hardLowerBound = lower;
        this.hardUpperBound = upper;
    }

    MeasureType(String[] namesInCSV, double lower, double upper, boolean isAccurateForRoot) {
        this.namesInInputFile = namesInCSV;
        this.hardLowerBound = lower;
        this.hardUpperBound = upper;
        this.isAccurateForRootCell = isAccurateForRoot;
    }

    MeasureType() {this.namesInInputFile = new String[]{""}; }

    public String[] getNames() {
        return namesInInputFile;
    }

    public boolean isAccurateMeasureForRootCell(){
        return this.isAccurateForRootCell;
    }

    public double getHardLowerBound(){
        return this.hardLowerBound;
    }

    public double getHardUpperBound(){
        return this.hardUpperBound;
    }

    /**
     * Check if experimental measure is present in the input csv file.
     * @param parser
     * @return
     */
    boolean isMeasureInCSV(CSVParser parser){
        if(this == InstantSpeed)
            return (XPosition.isMeasureInCSV(parser) && YPosition.isMeasureInCSV(parser));

        for (String name : this.namesInInputFile) {
            if(parser.getHeaderMap().containsKey(name))
                return true;
        }
        return false;
    }
}
