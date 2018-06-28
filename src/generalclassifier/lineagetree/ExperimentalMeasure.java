package generalclassifier.lineagetree;

import org.apache.commons.csv.CSVRecord;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class ExperimentalMeasure {

    static Map<MeasureType, CalculationMethod> mapMeasureToMethod = createMeasureMethodMap();

    private static Map<MeasureType, CalculationMethod> createMeasureMethodMap() {
        Map<MeasureType, CalculationMethod> result = new HashMap<MeasureType, CalculationMethod>();

        result.put(MeasureType.Eccentricity, CalculationMethod.averageValue);

        result.put(MeasureType.Area, CalculationMethod.averageRate);
        result.put(MeasureType.Perimeter, CalculationMethod.averageRate);
        result.put(MeasureType.Ig2afc, CalculationMethod.averageRate);
        result.put(MeasureType.Sca1, CalculationMethod.averageRate);
        result.put(MeasureType.CD41, CalculationMethod.averageRate);
        result.put(MeasureType.FcgRIII, CalculationMethod.averageRate);
        result.put(MeasureType.ROS, CalculationMethod.averageRate);
        result.put(MeasureType.TMRM, CalculationMethod.averageRate);
        result.put(MeasureType.CD71APC, CalculationMethod.averageRate);
        result.put(MeasureType.CD71PE, CalculationMethod.averageRate);

        result.put(MeasureType.cMycGFP, CalculationMethod.maxRate);

        result.put(MeasureType.ElapsedTime, CalculationMethod.differenceStartToEnd);
        result.put(MeasureType.Lifetime, CalculationMethod.differenceStartToEnd);

        result.put(MeasureType.XPosition, CalculationMethod.undefined);
        result.put(MeasureType.YPosition, CalculationMethod.undefined);

        result.put(MeasureType.AverageSpeed, CalculationMethod.averageInstantSpeed);

        return Collections.unmodifiableMap(result);
    }


    CalculationMethod calculationMethod = CalculationMethod.undefined;

    MeasureType measureType;

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


