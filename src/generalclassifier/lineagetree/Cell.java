package generalclassifier.lineagetree;

import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.math3.analysis.function.Exp;

import java.util.ArrayList;
import java.util.List;

public class Cell {

    public enum Fate {
        // D: divides, A: apoptoses, L: lost cell, N: does nothing (not sure we'll keep that), U: unobserved
        D, A, L, N, U
    }

    Fate fate;
    int trackNumber;

    List<double[]> timePoints = new ArrayList<> ();
    List<ExperimentalMeasure> measures = new ArrayList<> ();

    double lifetime;
    double area;
    double eccentricity;
    double perimeter;
    double averageSpeed;



    public Cell(){
    }

    public Cell(int trackNumber) {
        this.trackNumber = trackNumber;
    }

    public Cell(int trackNumber, List<ExperimentalMeasure> measures) {
        this.trackNumber = trackNumber;
        this.measures = measures;
    }

    public Cell(double branchlength) {
        this.lifetime = branchlength;
    }

    public void setMeasures(CSVParser parser) {
        for(MeasureType measureType : MeasureType.values()) {
            if(measureType.isMeasureInCSV(parser)) {
                ExperimentalMeasure measure = new ExperimentalMeasure(measureType);
                measures.add(measure);
            }
        }
    }

    public void setMeasures(List<ExperimentalMeasure> measures) {
        this.measures = measures;
    }

    public void addTimePoint(CSVRecord record) {
        double[] timePoint = new double[measures.size()];
        for (int i = 0; i < measures.size(); i++) {
            timePoint[i] = Double.parseDouble(measures.get(i).getValueInRecord(record));
        }
        timePoints.add(timePoint);
    }

    public double getLifetime() {
        return this.lifetime;
    }

    public void setLifetime(double branchLength) {
        this.lifetime = branchLength;
    }

    public Fate getFate() {
        return this.fate;
    }

    public void setFate(Fate f) {
        this.fate = f;
    }




//    public double getArea() {
//        return area;
//    }
//
//    public void setArea(double area) {
//        this.area = area;
//    }
//
//    public double getEccentricity() {
//        return eccentricity;
//    }
//
//    public void setEccentricity(double eccentricity) {
//        this.eccentricity = eccentricity;
//    }
//
//    public double getAverageSpeed() {
//        return averageSpeed;
//    }
//
//    public void setAverageSpeed(double averageSpeed) {
//        this.averageSpeed = averageSpeed;
//    }
//
//    public double getPerimeter() {
//        return perimeter;
//    }
//
//    public void setPerimeter(double perimeter) {
//        this.perimeter = perimeter;
//    }
}
