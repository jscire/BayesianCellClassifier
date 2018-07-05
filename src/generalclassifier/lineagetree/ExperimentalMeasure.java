package generalclassifier.lineagetree;

import org.apache.commons.csv.CSVRecord;

import java.util.*;

public class ExperimentalMeasure {

    static Map<MeasureType, CalculationMethod> mapMeasureToMethod = createMeasureMethodMap();

    private static Map<MeasureType, CalculationMethod> createMeasureMethodMap() {
        Map<MeasureType, CalculationMethod> result = new HashMap<>();

        result.put(MeasureType.Eccentricity, CalculationMethod.averageValue);

        result.put(MeasureType.Area, CalculationMethod.averageRate);
        result.put(MeasureType.Perimeter, CalculationMethod.averageRate);
        result.put(MeasureType.Ig2afc, CalculationMethod.averageRate);
        result.put(MeasureType.Sca1, CalculationMethod.averageRate);
        result.put(MeasureType.CD41, CalculationMethod.averageRate);
        result.put(MeasureType.FcgRIII, CalculationMethod.averageRate);
        result.put(MeasureType.ROS, CalculationMethod.averageRate);

        result.put(MeasureType.TMRMmean, CalculationMethod.averageRate);
        result.put(MeasureType.TMRMmax, CalculationMethod.maxRate);

        result.put(MeasureType.CD71APC, CalculationMethod.averageRate);
        result.put(MeasureType.CD71PE, CalculationMethod.averageRate);

        result.put(MeasureType.cMycGFP, CalculationMethod.maxRate);

        result.put(MeasureType.ElapsedTime, CalculationMethod.differenceStartToEnd);
        result.put(MeasureType.Lifetime, CalculationMethod.differenceStartToEnd);

        result.put(MeasureType.XPosition, CalculationMethod.undefined);
        result.put(MeasureType.YPosition, CalculationMethod.undefined);

        result.put(MeasureType.InstantSpeed, CalculationMethod.averageInstantSpeed);

        return Collections.unmodifiableMap(result);
    }


    CalculationMethod calculationMethod = CalculationMethod.undefined;

    MeasureType measureType;

    List<Double> dataPoints = new ArrayList<>();

    double summaryValue = Double.NaN;

    ExperimentalMeasure() {
        // do nothing
    }

    ExperimentalMeasure(MeasureType measureType) {
        this.measureType = measureType;
        this.calculationMethod = mapMeasureToMethod.get(measureType);
    }

    ExperimentalMeasure(MeasureType measureType, CalculationMethod method) {
        this.measureType = measureType;
        this.calculationMethod = method;
    }

    ExperimentalMeasure(CalculationMethod method) {
        this.calculationMethod = method;
    }

    CalculationMethod getCalculationMethod(){
        return this.calculationMethod;
    }


    void summarize(List<Double> timePoints) {
        double res = 0;

        if(dataPoints.size() != timePoints.size())
            throw new IllegalArgumentException("The number of measure points does not match the number of time points.");

        int lastIndex = timePoints.size()-1;

        switch (this.calculationMethod) {
            case maxRate :
                double instantRate;
                for (int i = 1; i <= lastIndex; i++) {
                    instantRate = (dataPoints.get(i) - dataPoints.get(i-1))/(timePoints.get(i) - timePoints.get(i-1));
                    if(instantRate > res)
                        res = instantRate;
                }
                break;

            case averageRate :
                res = (dataPoints.get(lastIndex) - dataPoints.get(0))/(timePoints.get(lastIndex) - timePoints.get(0));
                break;

            case averageValue:
                // sum all the observed values
                for (int i = 0; i <= lastIndex ; i++) {
                    res += dataPoints.get(i);
                }
                res /= (lastIndex + 1); // take the mean
                break;

            case differenceStartToEnd:
                res = dataPoints.get(lastIndex) - dataPoints.get(0);
                break;

            case averageInstantSpeed:
            case undefined:
            default:
                throw new IllegalArgumentException("Cannot summarize this measure. Method is: " + this.calculationMethod);

        }
        summaryValue = res;
    }




    /**
     *
     * @param record
     * @return value in record.
     * return null if field is not in csv file
     */
    String getValueInRecord(CSVRecord record){
        for (String name : this.measureType.getNames()) {
            if(record.isMapped(name))
                return record.get(name);
        }
        return null;
    }

    void addDataPointFromRecord(CSVRecord record) {
        for (String name : this.measureType.getNames()) {
            if(record.isMapped(name)) {
                dataPoints.add(Double.parseDouble(record.get(name)));
                break;
            }
        }
    }

//    boolean isMeasureInRecord(CSVRecord record) {
//        for (String name : this.namesInInputFile) {
//            if(record.isMapped(name))
//                return true;
//        }
//        return false;
//    }

    enum CalculationMethod {
        undefined,
        maxRate,
        averageRate,
        differenceStartToEnd,
        averageValue,
        averageInstantSpeed
    }
}


