package generalclassifier.lineagetree;

import org.apache.commons.csv.CSVParser;

public enum MeasureType {

    ElapsedTime(new String[]{"time_sec", "time_hours"}),
    Lifetime(new String[]{"TimePoint"}),
    Area(new String[]{"Area", "AreaMorphologyCh00"}, 0, Double.POSITIVE_INFINITY),
    Eccentricity(new String[]{"Eccentricity", "EccentricityMorphologyCh00"}, 0, 1),
    Perimeter(new String[]{"Perimeter", "PerimeterMorphologyCh00"}, 0 , Double.POSITIVE_INFINITY),
    XPosition(new String[]{"CentroidX","XMorphologyCh00"}),
    YPosition(new String[]{"CentroidY","YMorphologyCh00"}),
    InstantSpeed(new String[]{"Speed"}, 0, Double.POSITIVE_INFINITY),

    // fluorescent markers
    Ig2afc(new String[]{"Ig2afc", "Ig2afc_signal"}),
    Sca1(new String[]{"Sca1", "Sca1_signal"}),
    CD41(new String[]{"CD41", "CD41_signal"}),
    FcgRIII(new String[]{"FcgRIII", "FcgRIII_signal"}),
    cMycGFP(new String[]{"cMycGFP", "cMycGFP_signal"}),
    ROS(new String[]{"ROS", "ROS_signal"}),
    TMRMmean(new String[]{"TMRM", "TMRM_signal"}),
    TMRMmax(new String[]{"TMRM", "TMRM_signal"}),
    CD71APC(new String[]{"CD71APC", "CD71APC_signal"}),
    CD71PE(new String[]{"CD71PE", "CD71PE_signal"});

    private String[] namesInInputFile;
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

    MeasureType() {this.namesInInputFile = new String[]{""}; }

    public String[] getNames() {
        return namesInInputFile;
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
        for (String name : this.namesInInputFile) {
            if(parser.getHeaderMap().containsKey(name))
                return true;
        }
        return false;
    }
}
