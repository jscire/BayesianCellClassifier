package generalclassifier.utils;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Function;
import beast.core.Input.Validate;
import beast.core.Loggable;
import generalclassifier.core.LineageTree;
import generalclassifier.core.LineageTreeProb;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.BitSet;
import java.util.List;

@Description("Logger to report the probabilities associated with each possible configuration.")
public class ConfigurationsProbabilities extends CalculationNode implements Loggable,Function{

    public Input<LineageTreeProb> lineageTreeProbInput = new Input<>("lineageTreeProb",
            "Probability distr of the lineage tree of interest.", Validate.REQUIRED);


    static DecimalFormat df = new DecimalFormat("######0.00");
    static DecimalFormat df2 = new DecimalFormat("0.00E0");

    @Override
    public void initAndValidate(){
       // ltProb = lineageTreeProbInput.get();
    } //TODO, I think there's nothing to do in initValidate

    @Override
    public void init(PrintStream out) {
        out.print(logCorrespondanceCellLabelsToConfig()+ System.getProperty("line.separator"));
    }

    @Override
    public void log(int nSample, PrintStream out) {
        LineageTreeProb ltProb = lineageTreeProbInput.get();

        if(!ltProb.isHSCtree()) { // tree has MPP root
            out.print(nSample + "\t");
            out.print(ltProb.getTreeID() + "\t");

            double probaConfig = ltProb.getCurrentLogP();
            out.print(df.format(probaConfig) + "\t0\t"); // represent the MPP config by a single zero
        } else { // root is HSC and then there are multipleconfigurations to take into account
            List<BitSet> configsHSC = ltProb.getTypeConfigurationsHSC();
            double[] configsLogProbsHSC = ltProb.getCurrentConfigsLogProbsHSC();

            for (int i=0; i<configsHSC.size(); i++) {
                //if(i != 0) out.print(nSample + "\t");
//                out.print(nSample + "\t");
//                out.print(ltProb.getTreeID() + "\t");
//                out.print(df.format(ltProb.getProbabilityForConfiguration(configsHSC.get(i))) + "\t");
//                out.print(logConfig(configsHSC.get(i)) + System.getProperty("line.separator"));
                out.print(nSample + "\t");
                out.print(ltProb.getTreeID() + "\t");
                out.print(df.format(configsLogProbsHSC[i]) + "\t");
                out.print(logConfig(configsHSC.get(i)) + System.getProperty("line.separator"));

            }
        }
    }

    private String logConfig(BitSet configuration) {
        return MiscUtils.toString(configuration) + "\t";
    }

    public String logCorrespondanceCellLabelsToConfig() {
        LineageTreeProb ltProb = lineageTreeProbInput.get();
        return ltProb.getCorrespondanceLabelsToNodeArray();
    }


    @Override
    public void close(PrintStream out) { }


    @Override
    public int getDimension() {
        return 1;
    } // TODO implement method

    @Override
    public double getArrayValue() {
        return 0; //TODO implement method
    }

    @Override
    public double getArrayValue(int iDim) {
        return 0; //TODO implement method
    }

}
