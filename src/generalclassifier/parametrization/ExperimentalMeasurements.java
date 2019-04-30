package generalclassifier.parametrization;

import beast.core.BEASTObject;
import beast.core.Input;

import java.util.HashMap;

//TODO rename class
public class ExperimentalMeasurements extends BEASTObject {

    Input<String> tagInput = new Input<>("measurementTag", "Tag specifying the type of measurement.",
            Input.Validate.REQUIRED);

    Input<String> valuesInput = new Input<>("values", "Values format: " +
            "'..., N:XXXX.XX, N+1:YYYYY.YYY,...'" +
            " with N the cell number and XXXX.XXX the summary value for cell N. " +
            "If no value is given for a given cell, the cell is considered to have no observed " +
            "value for the particular phenotype trait.",
            Input.Validate.REQUIRED);

    HashMap<Integer, Double> measuredValues;

    String measurementTag;

    @Override
    public void initAndValidate() {
        measurementTag = tagInput.get();
        parseMeasuresInput();
    }

    public void parseMeasuresInput(){
        //remove all whitespaces and split on commas
        String[] parts = valuesInput.get().replaceAll("\\s","").split(",");
        for(String s : parts) {
            if(s.contains(":")) {
                try{
                    Integer cellNumber = Integer.parseInt(s.split(":")[0]);
                    Double measuredValue = Double.parseDouble(s.split(":")[1]);
                    measuredValues.put(cellNumber, measuredValue);
                }
                catch(Exception e) {
                    throw new IllegalArgumentException("Wrong input format, unable to parse 'values' input.");
                }
            }
            else throw new IllegalArgumentException("Wrong input format, unable to parse 'values' input, missing ':' character.");
        }
    }

    public String getMeasurementTag() {
        return measurementTag;
    }

    public HashMap<Integer, Double> getMeasuredValues(){
        return measuredValues;
    }
}
