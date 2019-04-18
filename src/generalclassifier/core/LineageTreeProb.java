package generalclassifier.core;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

import generalclassifier.lineagetree.Cell;
import generalclassifier.lineagetree.ExperimentalMeasure;
import generalclassifier.lineagetree.LineageTree;
import generalclassifier.lineagetree.MeasureType;
import generalclassifier.utils.InputGroup;
import generalclassifier.utils.Utils;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

import java.io.*;
import java.util.*;


public class  LineageTreeProb extends Distribution {

    public Input<LineageTree> lineageTreeInput = new Input<>("tree",
            "Lineage tree.",
            Input.Validate.REQUIRED);

    public Input<List<RealParameter>> fateProbabilitiesInput = new Input<>("fateProbabilities",
            "List of vectors containing probabilities for each cell type of the cell being a divider, apoptoser, non-divider or type-transitioner (if defined).",
            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);

    public Input<List<RealParameter>> scaleWeibullInput = new Input<>("scaleWeibull",
            "List of vectors containing Weibull scale parameters for each type for division, apoptosis and type-transition (if defined) times.",
            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);

    public Input<List<RealParameter>> shapeWeibullInput = new Input<>("shapeWeibull",
            "List of vectors containing Weibull shape parameters for each type for division, apoptosis and type-transition (if defined) times.",
            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);

    public Input<List<RealParameter>> transitionUponDivisionProbsInput = new Input<>("transitionUponDivisionProbs",
            "List of flattened 2-dimensional arrays containing the probabilites of transition between states upon division." +
                    "If mkj in array i represents the probability that a cell of type i divides into cells of type k and j (k<= j)," +
                    "the element mkj has index n*k + j - k(k+1)/2 in the input vector." +
                    "If j<k, invert k and j. We assume here that the daughter cells are not ordered and so that the arrays are symmetric before being truncated (the triangle below the diagonal is removed)." +
                    "This renumbering is equivalent to flattening this truncated array by rows." +
                    "Note that all vectors in this list must have elements which sum up to 1.",
            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);


    public Input<RealParameter> lossProbInput = new Input<>("lossProb",
            "Probability of losing any given cell",
            Input.Validate.REQUIRED);

    public Input<IntegerParameter> cellTypeInput = new Input<>("cellType", "Each number corresponds to the type of a cell. The index of a cell type is its tracknumber times the index of its tree." +
            "The cell type defined here refers to the type at the beginning of the tree edge that carries the cell. This is important if type transitions are to be allowed along tree edges." +
            "It is not the case yet.");

    //TODO replace in all xmls "rootTypeIdx" by "treeIdx" (do it when rerunning all the analysis of validation and RNAseq exp data)
    // and rootType by cellType and add "rootTypeOnly" if needed
    public Input<Integer> treeIdxInput = new Input<>("treeIdx",
            "If provided, this is an index into the rootIsHSC boolean " +
                    "parameter specifying which element corresponds to the root " +
                    "HSC status.");

    public Input<Boolean> rootTypeOnlyInput = new Input<>("rootTypeOnly", "If true, only the types of the root cells are inferred, the types of other cells are marginalized over. Default: true", true);

    public Input<BooleanParameter> matrixOfAllowedTransitionsInput = new Input<> ("allowedTransitions",
            "Flattened matrix of allowed transitions between cell types. " +
                    "If mij is the boolean that characterizes transitions from type i to j, the element of index i*(n-1)+j if i>j, or i*(n-1) + j-1 otherwise, gives the index of element mij. n is the number of types.");

    // TODO possibly rework on the way inputs are written in. This is a first draft
    // what about an additional class that takes care of it?
    public Input<RealParameter> meanLogNormalAreaGrowthRateInput = new Input<>("meanAreaGrowthRateLogNormal", "");
    public Input<RealParameter> sdLogNormalAreaGrowthRateInput = new Input<>("sdAreaGrowthRateLogNormal", "");

    public Input<RealParameter> meanLogNormalPerimeterGrowthRateInput = new Input<>("meanPerimeterGrowthRateLogNormal", "");
    public Input<RealParameter> sdLogNormalPerimeterGrowthRateInput = new Input<>("sdPerimeterGrowthRateLogNormal", "");

    public Input<RealParameter> meanNormalEccentricityInput = new Input<>("meanEccentricity", "");
    public Input<RealParameter> sdNormalEccentricityInput = new Input<>("sdEccentricity", "");

    public Input<RealParameter> meanLogNormalInstantSpeedInput = new Input<>("meanInstantSpeedLogNormal", "");
    public Input<RealParameter> sdLogNormalInstantSpeedInput = new Input<>("sdInstantSpeedLogNormal", "");

    public Input<RealParameter> meanLogNormalCD41ProductionRateInput = new Input<>("meanCD41ProductionRateLogNormal", "");
    public Input<RealParameter> sdLogNormalCD41ProductionRateInput = new Input<>("sdCD41ProductionRateLogNormal", "");

    public Input<RealParameter> meanLogNormalFcgRIIIProductionRateInput = new Input<>("meanFcgRIIIProductionRateLogNormal", "");
    public Input<RealParameter> sdLogNormalFcgRIIIProductionRateInput = new Input<>("sdFcgRIIIProductionRateLogNormal", "");

    public Input<RealParameter> meanLogNormalROSProductionRateInput = new Input<>("meanROSProductionRateLogNormal", "");
    public Input<RealParameter> sdLogNormalROSProductionRateInput = new Input<>("sdROSProductionRateLogNormal", "");

    public Input<RealParameter> meanLogNormalTMRMProductionRateInput = new Input<>("meanTMRMProductionRateLogNormal", "");
    public Input<RealParameter> sdLogNormalTMRMProductionRateInput = new Input<>("sdTMRMProductionRateLogNormal", "");

    public Input<RealParameter> meanLogNormalTMRMMaxRateInput = new Input<>("meanTMRMMaxRateLogNormal", "");
    public Input<RealParameter> sdLogNormalTMRMMaxRateInput = new Input<>("sdTMRMMaxRateLogNormal", "");

    public Input<RealParameter> meanLogNormalSca1MeanRateInput = new Input<>("meanSca1MeanRateLogNormal", "");
    public Input<RealParameter> sdLogNormalSca1MeanRateInput = new Input<>("sdSca1MeanRateLogNormal", "");

    public Input<RealParameter> meanLogNormalSca1MeanValueInput = new Input<>("meanSca1MeanValueLogNormal", "");
    public Input<RealParameter> sdLogNormalSca1MeanValueInput = new Input<>("sdSca1MeanValueLogNormal", "");

    public Input<RealParameter> meanLogNormalIg2afcProductionRateInput = new Input<>("meanIg2afcProductionRateLogNormal", "");
    public Input<RealParameter> sdLogNormalIg2afcProductionRateInput = new Input<>("sdIg2afcProductionRateLogNormal", "");

    public Input<RealParameter> meanLogNormalCD71APCProductionRateInput = new Input<>("meanCD71APCProductionRateLogNormal", "");
    public Input<RealParameter> sdLogNormalCD71APCProductionRateInput = new Input<>("sdCD71APCProductionRateLogNormal", "");

    public Input<RealParameter> meanLogNormalCD71PEProductionRateInput = new Input<>("meanCD71PEProductionRateLogNormal", "");
    public Input<RealParameter> sdLogNormalCD71PEProductionRateInput = new Input<>("sdCD71PEProductionRateLogNormal", "");

    public Input<RealParameter> meanLogNormalcMycGFPMaxRateInput = new Input<>("meancMycGFPMaxRateLogNormal", "");
    public Input<RealParameter> sdLogNormalcMycGFPMaxRateInput = new Input<>("sdcMycGFPMaxRateLogNormal", "");


    public Input<RealParameter> zeroFractionAreaGrowthRateInput = new Input<>("zeroFractionAreaGrowthRate", "");

    public Input<RealParameter> zeroFractionPerimeterGrowthRateInput = new Input<>("zeroFractionPerimeterGrowthRate", "");

    public Input<RealParameter> zeroFractionInstantSpeedInput = new Input<>("zeroFractionInstantSpeed", "");

    public Input<RealParameter> zeroFractionCD41ProductionRateInput = new Input<>("zeroFractionCD41ProductionRate", "");

    public Input<RealParameter> zeroFractionLysobriteProductionRateInput = new Input<>("zeroFractionLysobriteProductionRate", "");

    public Input<RealParameter> zeroFractionFcgRIIIProductionRateInput = new Input<>("zeroFractionFcgRIIIProductionRate", "");

    public Input<RealParameter> zeroFractionTMRMProductionRateInput = new Input<>("zeroFractionTMRMProductionRate", "");

    public Input<RealParameter> zeroFractionTMRMMaxRateInput = new Input<>("zeroFractionTMRMMaxRate", "");

    public Input<RealParameter> zeroFractionROSProductionRateInput = new Input<>("zeroFractionROSProductionRate", "");

    public Input<RealParameter> zeroFractionSca1MeanRateInput = new Input<>("zeroFractionSca1ProductionRate", "");

    public Input<RealParameter> zeroFractionSca1MeanValueInput = new Input<>("zeroFractionSca1MeanValue", "");

    public Input<RealParameter> zeroFractionIg2afcProductionRateInput = new Input<>("zeroFractionIg2afcProductionRate", "");

    public Input<RealParameter> zeroFractionCD71APCProductionRateInput = new Input<>("zeroFractionCD71APCProductionRate", "");

    public Input<RealParameter> zeroFractionCD71PEProductionRateInput = new Input<>("zeroFractionCD71PEProductionRate", "");

    public Input<RealParameter> zeroFractioncMycGFPMaxRateInput = new Input<>("zeroFractioncMycGFPMaxRate", "");




    public Input<RealParameter> shapeGammaAreaGrowthRateInput = new Input<>("shapeAreaGrowthRate",
            "");
    public Input<RealParameter> meanGammaAreaGrowthRateInput = new Input<>("meanAreaGrowthRate",
            "");
    public Input<RealParameter> shapeGammaPerimeterGrowthRateInput = new Input<>("shapePerimeterGrowthRate",
            "");
    public Input<RealParameter> meanGammaPerimeterGrowthRateInput = new Input<>("meanPerimeterGrowthRate",
            "");
    public Input<RealParameter> shapeGammaInstantSpeedInput = new Input<>("shapeInstantSpeed",
            "");
    public Input<RealParameter> meanGammaInstantSpeedInput = new Input<>("meanInstantSpeed",
            "");
    public Input<RealParameter> shapeGammaCD41ProductionRateInput = new Input<>("shapeCD41ProductionRate",
            "");
    public Input<RealParameter> meanGammaCD41ProductionRateInput = new Input<>("meanCD41ProductionRate",
            "");
    public Input<RealParameter> shapeGammaLysobriteProductionRateInput = new Input<>("shapeLysobriteProductionRate",
            "");
    public Input<RealParameter> meanGammaLysobriteProductionRateInput = new Input<>("meanLysobriteProductionRate",
            "");
    public Input<RealParameter> shapeGammaFcgRIIIProductionRateInput = new Input<>("shapeFcgRIIIProductionRate",
            "");
    public Input<RealParameter> meanGammaFcgRIIIProductionRateInput = new Input<>("meanFcgRIIIProductionRate",
            "");
    public Input<RealParameter> shapeGammaROSProductionRateInput = new Input<>("shapeROSProductionRate",
            "");
    public Input<RealParameter> meanGammaROSProductionRateInput = new Input<>("meanROSProductionRate",
            "");
    public Input<RealParameter> shapeGammaTMRMProductionRateInput = new Input<>("shapeTMRMProductionRate",
            "");
    public Input<RealParameter> meanGammaTMRMProductionRateInput = new Input<>("meanTMRMProductionRate",
            "");
    public Input<RealParameter> shapeGammaTMRMMaxRateInput = new Input<>("shapeTMRMMaxRate",
            "");
    public Input<RealParameter> meanGammaTMRMMaxRateInput = new Input<>("meanTMRMMaxRate",
            "");
    public Input<RealParameter> shapeGammaSca1ProductionRateInput = new Input<>("shapeSca1ProductionRate",
            "");
    public Input<RealParameter> meanGammaSca1ProductionRateInput = new Input<>("meanSca1ProductionRate",
            "");
    public Input<RealParameter> shapeGammaSca1MeanValueInput = new Input<>("shapeSca1MeanValue",
            "");
    public Input<RealParameter> meanGammaSca1MeanValueInput = new Input<>("meanSca1MeanValue",
            "");
    public Input<RealParameter> shapeGammaIg2afcProductionRateInput = new Input<>("shapeIg2afcProductionRate",
            "");
    public Input<RealParameter> meanGammaIg2afcProductionRateInput = new Input<>("meanIg2afcProductionRate",
            "");
    public Input<RealParameter> shapeGammaCD71APCProductionRateInput = new Input<>("shapeCD71APCProductionRate",
            "");
    public Input<RealParameter> meanGammaCD71APCProductionRateInput = new Input<>("meanCD71APCProductionRate",
            "");
    public Input<RealParameter> shapeGammaCD71PEProductionRateInput = new Input<>("shapeCD71PEProductionRate",
            "");
    public Input<RealParameter> meanGammaCD71PEProductionRateInput = new Input<>("meanCD71PEProductionRate",
            "");
    public Input<RealParameter> shapeGammacMycGFPMaxRateInput = new Input<>("shapecMycGFPMaxRate",
            "");
    public Input<RealParameter> meanGammacMycGFPMaxRateInput = new Input<>("meancMycGFPMaxRate",
            "");


    //mot50
    public Input<RealParameter> alphaBetaMot50Input = new Input<>("alphaMot50", "");
    public Input<RealParameter> betaBetaMot50Input = new Input<>("betaMot50", "");

    Map<InputGroup, MeasureType> mapMeasureTypeToInput = createMapOfMeasureInputs();

    //TODO change how the map is created (too much duplication between the distributions)
    private Map<InputGroup, MeasureType> createMapOfMeasureInputs(){
        Map<InputGroup, MeasureType> result = new HashMap<>();

        InputGroup areaGrowthRateLogNormInput = new InputGroup(meanLogNormalAreaGrowthRateInput, sdLogNormalAreaGrowthRateInput, zeroFractionAreaGrowthRateInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup perimeterGrowthRateLogNormInput = new InputGroup(meanLogNormalPerimeterGrowthRateInput, sdLogNormalPerimeterGrowthRateInput, zeroFractionPerimeterGrowthRateInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup eccentricityInput = new InputGroup(meanNormalEccentricityInput, sdNormalEccentricityInput, InputGroup.DistributionType.NORMAL);

        InputGroup instantSpeedLogNormInput = new InputGroup(meanLogNormalInstantSpeedInput, sdLogNormalInstantSpeedInput, zeroFractionInstantSpeedInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup CD41ProductionRateLogNormInput = new InputGroup(meanLogNormalCD41ProductionRateInput, sdLogNormalCD41ProductionRateInput, zeroFractionCD41ProductionRateInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup FcgRIIIProductionRateLogNormInput = new InputGroup(meanLogNormalFcgRIIIProductionRateInput, sdLogNormalFcgRIIIProductionRateInput, zeroFractionFcgRIIIProductionRateInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup ROSProductionRateLogNormInput = new InputGroup(meanLogNormalROSProductionRateInput, sdLogNormalROSProductionRateInput, zeroFractionROSProductionRateInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup TMRMProductionRateLogNormInput = new InputGroup(meanLogNormalTMRMProductionRateInput, sdLogNormalTMRMProductionRateInput, zeroFractionTMRMProductionRateInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup TMRMMaxRateLogNormInput = new InputGroup(meanLogNormalTMRMMaxRateInput, sdLogNormalTMRMMaxRateInput, zeroFractionTMRMMaxRateInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup Sca1MeanRateLogNormInput = new InputGroup(meanLogNormalSca1MeanRateInput, sdLogNormalSca1MeanRateInput, zeroFractionSca1MeanRateInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup Sca1MeanValueLogNormInput = new InputGroup(meanLogNormalSca1MeanValueInput, sdLogNormalSca1MeanValueInput, zeroFractionSca1MeanValueInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup Ig2afcProductionRateLogNormInput = new InputGroup(meanLogNormalIg2afcProductionRateInput, sdLogNormalIg2afcProductionRateInput, zeroFractionIg2afcProductionRateInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup CD71APCProductionRateLogNormInput = new InputGroup(meanLogNormalCD71APCProductionRateInput, sdLogNormalCD71APCProductionRateInput, zeroFractionCD71APCProductionRateInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup CD71PEProductionRateLogNormInput = new InputGroup(meanLogNormalCD71PEProductionRateInput, sdLogNormalCD71PEProductionRateInput, zeroFractionCD71PEProductionRateInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup cMycGFPMaxRateLogNormInput = new InputGroup(meanLogNormalcMycGFPMaxRateInput, sdLogNormalcMycGFPMaxRateInput, zeroFractioncMycGFPMaxRateInput,
                InputGroup.DistributionType.LOGNORMAL);

        InputGroup areaGrowthRateGammaInput = new InputGroup(shapeGammaAreaGrowthRateInput, meanGammaAreaGrowthRateInput, zeroFractionAreaGrowthRateInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup perimeterGrowthRateGammaInput = new InputGroup(shapeGammaPerimeterGrowthRateInput, meanGammaPerimeterGrowthRateInput, zeroFractionPerimeterGrowthRateInput,
                InputGroup.DistributionType.GAMMA);
        
        InputGroup instantSpeedGammaInput = new InputGroup(shapeGammaInstantSpeedInput, meanGammaInstantSpeedInput, zeroFractionInstantSpeedInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup CD41ProductionRateGammaInput = new InputGroup(shapeGammaCD41ProductionRateInput, meanGammaCD41ProductionRateInput, zeroFractionCD41ProductionRateInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup LysobriteProductionRateGammaInput = new InputGroup(shapeGammaLysobriteProductionRateInput, meanGammaLysobriteProductionRateInput, zeroFractionLysobriteProductionRateInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup FcgRIIIProductionRateGammaInput = new InputGroup(shapeGammaFcgRIIIProductionRateInput, meanGammaFcgRIIIProductionRateInput, zeroFractionFcgRIIIProductionRateInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup ROSProductionRateGammaInput = new InputGroup(shapeGammaROSProductionRateInput, meanGammaROSProductionRateInput, zeroFractionROSProductionRateInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup TMRMProductionRateGammaInput = new InputGroup(shapeGammaTMRMProductionRateInput, meanGammaTMRMProductionRateInput, zeroFractionTMRMProductionRateInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup TMRMMaxRateGammaInput = new InputGroup(shapeGammaTMRMMaxRateInput, meanGammaTMRMMaxRateInput, zeroFractionTMRMMaxRateInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup Sca1MeanRateGammaInput = new InputGroup(shapeGammaSca1ProductionRateInput, meanGammaSca1ProductionRateInput, zeroFractionSca1MeanRateInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup Sca1MeanValueGammaInput = new InputGroup(shapeGammaSca1MeanValueInput, meanGammaSca1MeanValueInput, zeroFractionSca1MeanValueInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup Ig2afcProductionRateGammaInput = new InputGroup(shapeGammaIg2afcProductionRateInput, meanGammaIg2afcProductionRateInput, zeroFractionIg2afcProductionRateInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup CD71APCProductionRateGammaInput = new InputGroup(shapeGammaCD71APCProductionRateInput, meanGammaCD71APCProductionRateInput, zeroFractionCD71APCProductionRateInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup CD71PEProductionRateGammaInput = new InputGroup(shapeGammaCD71PEProductionRateInput, meanGammaCD71PEProductionRateInput, zeroFractionCD71PEProductionRateInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup cMycGFPMaxRateGammaInput = new InputGroup(shapeGammacMycGFPMaxRateInput, meanGammacMycGFPMaxRateInput, zeroFractioncMycGFPMaxRateInput,
                InputGroup.DistributionType.GAMMA);

        InputGroup mot50BetaInput = new InputGroup(alphaBetaMot50Input, betaBetaMot50Input, InputGroup.DistributionType.BETA);
        

        result.put(areaGrowthRateLogNormInput, MeasureType.Area);
        result.put(perimeterGrowthRateLogNormInput, MeasureType.Perimeter);
        result.put(eccentricityInput, MeasureType.Eccentricity);
        result.put(instantSpeedLogNormInput, MeasureType.InstantSpeed);
        result.put(CD41ProductionRateLogNormInput, MeasureType.CD41);
        result.put(FcgRIIIProductionRateLogNormInput, MeasureType.FcgRIII);
        result.put(ROSProductionRateLogNormInput, MeasureType.ROS);
        result.put(TMRMProductionRateLogNormInput, MeasureType.TMRMmean);
        result.put(TMRMMaxRateLogNormInput, MeasureType.TMRMmax);
        result.put(Sca1MeanRateLogNormInput, MeasureType.Sca1meanRate);
        result.put(Sca1MeanValueLogNormInput, MeasureType.Sca1meanValue);
        result.put(Ig2afcProductionRateLogNormInput, MeasureType.Ig2afc);
        result.put(CD71APCProductionRateLogNormInput, MeasureType.CD71APC);
        result.put(CD71PEProductionRateLogNormInput, MeasureType.CD71PE);
        result.put(cMycGFPMaxRateLogNormInput, MeasureType.cMycGFP);

        result.put(areaGrowthRateGammaInput, MeasureType.Area);
        result.put(perimeterGrowthRateGammaInput, MeasureType.Perimeter);
        result.put(instantSpeedGammaInput, MeasureType.InstantSpeed);
        result.put(CD41ProductionRateGammaInput, MeasureType.CD41);
        result.put(LysobriteProductionRateGammaInput, MeasureType.Lysobrite);
        result.put(FcgRIIIProductionRateGammaInput, MeasureType.FcgRIII);
        result.put(ROSProductionRateGammaInput, MeasureType.ROS);
        result.put(TMRMProductionRateGammaInput, MeasureType.TMRMmean);
        result.put(TMRMMaxRateGammaInput, MeasureType.TMRMmax);
        result.put(Sca1MeanRateGammaInput, MeasureType.Sca1meanRate);
        result.put(Sca1MeanValueGammaInput, MeasureType.Sca1meanValue);
        result.put(Ig2afcProductionRateGammaInput, MeasureType.Ig2afc);
        result.put(CD71APCProductionRateGammaInput, MeasureType.CD71APC);
        result.put(CD71PEProductionRateGammaInput, MeasureType.CD71PE);
        result.put(cMycGFPMaxRateGammaInput, MeasureType.cMycGFP);

        result.put(mot50BetaInput, MeasureType.mot50);


        return Collections.unmodifiableMap(result);
    }

    LineageTree tree;
    
    boolean transitionUponDivisionIsAllowed;
    boolean transitionDuringLifetimeIsAllowed;
    boolean sumOverDaughterCellTypes;

    int numberOfTypes;

    // TrapezoidIntegrator numericalIntegrator = new TrapezoidIntegrator(1e-6, 1e-6, 1, 60);
    IterativeLegendreGaussIntegrator numericalIntegrator = new IterativeLegendreGaussIntegrator(2, 1e-6, 1e-6);

    int maxEval = 10000;


    //TODO change code so that there is a matrix of transition probabilities for transitions on branches

//    public LineageTreeProb() {
//    }

    @Override
    public void initAndValidate() {

        tree = lineageTreeInput.get();

        numberOfTypes = fateProbabilitiesInput.get().size();

        //TODO remove when tested with more than 2 types.
//        if(numberOfTypes != 2)
//            throw new IllegalStateException("Number of types is not two. Such a configuration is not tested yet.");

        if(shapeWeibullInput.get().size() != numberOfTypes || scaleWeibullInput.get().size() != numberOfTypes)
            throw new IllegalStateException("Inputs fateProbabilities, shapeWeibull and scaleWeibull should have the same number of elements: the number of cell types.");

        if(matrixOfAllowedTransitionsInput.get() != null) {
            transitionDuringLifetimeIsAllowed = true;
            throw new IllegalStateException("Transitions on branches are not fully implemented yet");

            //TODO uncomment check below when exception above is not thrown anymore (when transitions on branches are implemented)
//            if(matrixOfAllowedTransitionsInput.get().getDimension() != numberOfTypes * (numberOfTypes -1))
//                throw new IllegalStateException("Incorrect size of matrix of allowed transitions.");
        }
        else {
            transitionDuringLifetimeIsAllowed = false;
        }

        if(transitionUponDivisionProbsInput.get().get(0) != null) {
            transitionUponDivisionIsAllowed=true;

            if(transitionUponDivisionProbsInput.get().size() != numberOfTypes)
                throw new IllegalStateException("Incorrect size of list of probabilities of transition upon division.");
        }
        else {
            transitionUponDivisionIsAllowed = false;
        }

        sumOverDaughterCellTypes = rootTypeOnlyInput.get().booleanValue();

        if (lossProbInput.get() == null)
            throw new IllegalStateException("lossProb input is missing. Set to zero if cells are never lost.");

        //TODO add checks for parameter input format. In particular, check that fateProbabilities, transitionProbsUponDivision, shapeWeibull and scaleWeibull have correct dimension for each cell type.

    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        //for now prior is that frequency is equal for all types TODO allow for more flexibility
        try {
            double p = 0;
            double[] resultPruning = calculatePruningProb((Cell) tree.getRoot());
            for (int i = 0; i < numberOfTypes; i++) {
                p += (resultPruning[i] * 1.0/numberOfTypes);
            }
            logP = Math.log(p);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return logP;
    }

    /**
     * Calculate the probability for current node conditional on a type nodeType using pruning algorithm
     * Assumes the node is not a leaf
     * Type transitions on a branch are possible. So method makes a difference between type at start and end of branch.
     * @param node
     * @return
     */
    public double[] calculatePruningProb(Cell node) { // here nodeType refers to the type of the node at the beginning of the branch

        int nodeType = getFixedCellType(node);

        if(node.isLeaf()) {
            // if a node is a leaf here then it is the root.
            // Otherwise leaves are caught below in a previous recursive call of calculatePruningProb
            return getProbaAtLeaves(node,-1);
        }

        int childIndex = 0;
        double[] pruningProbaFirstChild;

        if (node.getChild(childIndex).isLeaf())
            pruningProbaFirstChild = getProbaAtLeaves((Cell) node.getChild(childIndex),  -1);
        else
            pruningProbaFirstChild = calculatePruningProb((Cell) node.getChild(childIndex));

        childIndex = 1;
        double[] pruningProbaSecondChild;

        if (node.getChild(childIndex).isLeaf())
            pruningProbaSecondChild = getProbaAtLeaves((Cell) node.getChild(childIndex), -1);
        else
            pruningProbaSecondChild = calculatePruningProb((Cell) node.getChild(childIndex));


        // for each type, compute the conditional probability at the end of the edge, using the pruning algorithm
        double[] probsEndBranch = new double[numberOfTypes];
        for (int i = 0; i < numberOfTypes; i++) {
            for (int j = 0; j < numberOfTypes; j++) {
                probsEndBranch[i] += 2 * getProbabilityAtDivisionNode(i, j, j) * pruningProbaFirstChild[j] * pruningProbaSecondChild[j];
                for (int k = (j + 1); k < numberOfTypes; k++) {
                    probsEndBranch[i] += getProbabilityAtDivisionNode(i, j, k) *
                            (pruningProbaFirstChild[j] * pruningProbaSecondChild[k] + pruningProbaFirstChild[k] * pruningProbaSecondChild[j]);
                }
            }
        }

        // now compute conditional probability at the start of the edge,
        // given the probabilities at the end of edges that were calculated immediately above
        double[] probsStartBranch = new double[numberOfTypes];

        if(nodeType == -1) { // type of node is not fixed
            for (int i = 0; i < numberOfTypes; i++) {
                for (int j = 0; j < numberOfTypes; j++) {
                    probsStartBranch[i] += getProbabilityCellBranch(node, i, j) * probsEndBranch[j];
                }
            }
        }
        else  { // type of this node is fixed
            for (int j = 0; j < numberOfTypes; j++) {
                probsStartBranch[nodeType] += getProbabilityCellBranch(node, nodeType, j) * probsEndBranch[j];
            }
        }

        return probsStartBranch;
    }

    /**
     * For each potential cell type, compute conditional probability cell is in that type at beginning of branch above a leaf.
     * 2nd argument: type of the cell at the end of the tip (-1 if not known).
     * For now there is no way to specify in the inputs of an analysis that it is not -1.
     * I have not yet encountered a case where soecifying the type at the end was necessary (but could be useful if other type of data (result of RNAseq for instance is provided)
     * @param leaf
     * @param endBranchType
     * @return
     */
    public double[] getProbaAtLeaves(Cell leaf, int endBranchType) {

        double[] startingProbas  = new double[numberOfTypes];
        double[] resultProbs = new double[numberOfTypes];

        int startBranchType = getFixedCellType(leaf);

        if(startBranchType >= numberOfTypes || startBranchType < -1 || endBranchType >= numberOfTypes || endBranchType < -1) {
            throw new IllegalStateException("Undefined start or end branch type." +
                    " Value must be -1 if unspecified, and -1 < i < numberOfTypes otherwise");
        }

        if(endBranchType == -1) { // type of leaf at end of branch is not fixed
            // initialize all the starting probabilities (at end of branch)
            for (int i = 0; i < numberOfTypes; i++) {
                startingProbas[i] = 1.0;
            }

            if(startBranchType == -1) { // type of leaf at start of branch is not fixed
                for (int i = 0; i < numberOfTypes; i++) {
                    resultProbs[i] = 0;
                    for (int j = 0; j < numberOfTypes; j++) {
                        resultProbs[i] += getProbabilityCellBranch(leaf, i, j) * startingProbas[j];
                    }
                }
            }
            else { // type of leaf at start of branch is fixed
                for (int j = 0; j < numberOfTypes; j++) {
                    resultProbs[startBranchType] += getProbabilityCellBranch(leaf, startBranchType, j) * startingProbas[j];
                }
            }

        }
        else { // type of this leaf at end of branch is fixed
            startingProbas[endBranchType] = 1.0;

            if(startBranchType == -1) { // type of leaf at start of branch is not fixed
                for (int i = 0; i < numberOfTypes; i++) {
                    resultProbs[i] = 0;
                    for (int j = 0; j < numberOfTypes; j++) {
                        resultProbs[i] += getProbabilityCellBranch(leaf, i, j) * startingProbas[j];
                    }
                }
            }
            else {
                for (int j = 0; j < numberOfTypes; j++) {
                    resultProbs[startBranchType] += getProbabilityCellBranch(leaf, startBranchType, j) * startingProbas[j];
                }
            }
        }
        return resultProbs;
    }

    public double getProbabilityCellBranch(Cell node, int typeStartBranch, int typeEndBranch) {

        double branchProb=0;

        if(node.getFate() == Cell.Fate.L) return lossProbInput.get().getValue();

        if(typeStartBranch == typeEndBranch) { // no type transition on branch

            if (node.getFate() == Cell.Fate.D || node.getFate() == Cell.Fate.A) {

                int cellFate = node.getFate() == Cell.Fate.D ? 0 : 1; // get the index of the cellFate of the cell (0 for dividers, 1 for cells which undergo apotosis)

                if(node.isRoot()) {
                    // if cell is root cell, take its observed lifetime as the minimal bound of its real lifetime
                    branchProb = fateProbabilitiesInput.get().get(typeEndBranch).getArrayValue(cellFate)
                            * Math.exp(-Math.pow(node.getEdgeLength() / scaleWeibullInput.get().get(typeEndBranch).getArrayValue(cellFate), shapeWeibullInput.get().get(typeEndBranch).getArrayValue(cellFate)));
                }
                else {
                    // observed lifetime is the real lifetime of the cell
                    branchProb = fateProbabilitiesInput.get().get(typeEndBranch).getArrayValue(cellFate)
                            * Utils.getWeibullDensity(node.getEdgeLength(),
                            scaleWeibullInput.get().get(typeEndBranch).getArrayValue(cellFate),
                            shapeWeibullInput.get().get(typeEndBranch).getArrayValue(cellFate));
                }
            }
            else if (node.getFate() == Cell.Fate.U) { // fate of the cell is not observed in the observation window, therefore lifetime is at least as long as observation window
                // two possibilities: cell either divides or dies after end of observation window.

                branchProb = fateProbabilitiesInput.get().get(typeEndBranch).getValue(0) // cell divides after end of branch
                        * Math.exp(-Math.pow(node.getEdgeLength() / scaleWeibullInput.get().get(typeEndBranch).getValue(0), shapeWeibullInput.get().get(typeEndBranch).getValue(0)));

                if(fateProbabilitiesInput.get().get(typeEndBranch).getValue(1) > 0) // possibility cell dies after end of branch
                    branchProb += fateProbabilitiesInput.get().get(typeEndBranch).getValue(1)
                        * Math.exp(-Math.pow(node.getEdgeLength() / scaleWeibullInput.get().get(typeEndBranch).getValue(1), shapeWeibullInput.get().get(typeEndBranch).getValue(1)));

                if (transitionDuringLifetimeIsAllowed && fateProbabilitiesInput.get().get(typeEndBranch).getDimension() > 3) // check if this state can transition
                    //TODO if more than 1 fate that the cell can transition to, take it into account, either by summing the different probas or by having one big proba for all
                    branchProb += fateProbabilitiesInput.get().get(typeEndBranch).getValue(3) // cell transitions after end of branch
                            * Math.exp(-Math.pow(node.getEdgeLength() / scaleWeibullInput.get().get(typeEndBranch).getValue(2), shapeWeibullInput.get().get(typeEndBranch).getValue(2)));

            }
            else {
                throw new IllegalStateException("Invalid cell fate.");
            }

            // account for the experimental measures performed on the cells
            branchProb *= getProbabilityCellMeasures(node, typeStartBranch, typeEndBranch);

        }
        else { // typeStartBranch != typeEndBranch

            //TODO transition on branches are not properly implemented yet.

            //get index of transition in the matrixOfAllowedTransitions
            int indexTransition = typeStartBranch > typeEndBranch ? typeStartBranch * (numberOfTypes - 1) + typeEndBranch : typeStartBranch * (numberOfTypes - 1) + typeEndBranch - 1;

            if (!transitionDuringLifetimeIsAllowed || !matrixOfAllowedTransitionsInput.get().getValue(indexTransition))
                return 0; // configuration has probability zero (impossible type transition)
            else {
                if (node.getFate() == Cell.Fate.D || node.getFate() == Cell.Fate.A) {
                    int cellEndFate = node.getFate() == Cell.Fate.D ? 0 : 1; // get the index of the cellFate of the cell (0 for dividers, 1 for apoptosers)

                    // integration over transition point
                    // the "2" below refers to the fate transitioner. TODO refactor this "2", make cleaner which fate is which, when generalising
                    DensityProbOfStateTransition f = new DensityProbOfStateTransition(
                            node.getEdgeLength(),
                            shapeWeibullInput.get().get(typeStartBranch).getValue(2),
                            scaleWeibullInput.get().get(typeStartBranch).getValue(2),
                            shapeWeibullInput.get().get(typeEndBranch).getValue(cellEndFate),
                            scaleWeibullInput.get().get(typeEndBranch).getValue(cellEndFate));

                    try {
                        double integrationRes = numericalIntegrator.integrate(maxEval, f, 0, 1);
                        // the "3" below refers to the fate transitioner. TODO refactor this "3", make cleaner which fate is which, when generalising
                        branchProb = integrationRes * fateProbabilitiesInput.get().get(typeStartBranch).getValue(3) * fateProbabilitiesInput.get().get(typeEndBranch).getValue(cellEndFate);
                    }
                    catch (Exception e) { //TODO solve the pb with some integrations that don't happen: it's a pb with properly exploring the state space
                        branchProb = 0;
                    }

                }
                else if (node.getFate() == Cell.Fate.U) {
                    //TODO implement.
                    // will be similar to fate N but removing the possibility that fate N appears.
                }
                else {
                    throw new IllegalStateException("Invalid cell fate.");
                }
            }


        }
        branchProb *= (1 - lossProbInput.get().getValue());
        return branchProb;
    }

    public double getProbabilityCellMeasures(Cell node, int typeStartBranch, int typeEndBranch) {
        double branchProb = 1;

        //TODO implement transitions on branches
        if(typeEndBranch != typeStartBranch) {
            return 0;
        }

        for (InputGroup input : mapMeasureTypeToInput.keySet()) {
            MeasureType measureType = mapMeasureTypeToInput.get(input);

            // if input is provided by user and the cell actually contains info about the measure
            if(input.getParm1().get() != null && node.hasMeasure(measureType)) {
                double x = node.getSummaryValue(measureType);

                if(Double.isNaN(x)) continue; // skip this measure if there is no summary value

                boolean summaryIsAverage = ExperimentalMeasure.isSummaryAverage(ExperimentalMeasure.getCalculationMethodOfMeasureType(measureType));

                //TODO do this better because this assumes we are interested in the maximum if not interested in average. In particular, nothing is done for mot50 at the moment
                if(node.getFate() == Cell.Fate.U && !summaryIsAverage)
                    branchProb *= input.getComplementaryCumulativeDistribution(typeEndBranch, x);
                else
                    branchProb *= input.getProbabilityDensity(typeEndBranch, x);
            }
        }
        return branchProb;
    }

    public double getProbabilityAtDivisionNode(int typeEndParentBranch, int typeStartChildBranch1, int typeStartChildBranch2) {

        if (!transitionUponDivisionIsAllowed) {
            if (typeEndParentBranch == typeStartChildBranch1 && typeStartChildBranch1 == typeStartChildBranch2)
                return 1;
            else return 0; // if one of the three states is different, local configuration is impossible without transition upon division
        }
        else {

            if(typeStartChildBranch1 > typeStartChildBranch2) {
                int temp = typeStartChildBranch2;
                typeStartChildBranch2 = typeStartChildBranch1;
                typeStartChildBranch1 = temp;
            }

            int indexInFlattenedArray = numberOfTypes * typeStartChildBranch1 + typeStartChildBranch2 - typeStartChildBranch1 * (typeStartChildBranch1+1)/2;

            return transitionUponDivisionProbsInput.get().get(typeEndParentBranch).getValue(indexInFlattenedArray);
        }
    }

    //TODO implement check that treeIdxInputis long enough to get value

    /**
     * Return the fixed cell type at beginning of tree edge
     * returns -1 if type of cell is not fixed
     * !!! The type of the root must always be fixed !!!! (with the current code, but not necessary to keep it that way)
     * @param cell
     * @return
     */
    private int getFixedCellType(Cell cell){
        if(!sumOverDaughterCellTypes)
            return cellTypeInput.get().getValue(treeIdxInput.get() * cell.getTrackNumber());
        else if(cell.getTrackNumber() == 1) // cell is root cell
            return cellTypeInput.get().getValue(treeIdxInput.get());
        else return -1;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        throw new RuntimeException("Not implemented.");
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    public String getTreeID() {
        if (treeIdxInput.get() != null)
            return treeIdxInput.get() + "";
        else return "";
    }


    public String getCorrespondanceLabelsToNodeArray() {

        String result = "";
        for (Node node : tree.getNodesAsArray()) {
            int newNodeNumber = Integer.parseInt(node.getID().substring(0, 1));
            result += newNodeNumber;
        }
        return result;
    }

    private class DensityProbOfStateTransition implements UnivariateFunction {

        double cellLifetime;

        double shapeLifetimeOriginState;
        double shapeLifetimeEndState;

        double scaleLifetimeOriginState;
        double scaleLifetimeEndState;

        DensityProbOfStateTransition(double cellLifetime, double shapeWeibullHSC, double scaleWeibullHSC, double shapeWeibullMPP, double scaleWeibullMPP) {
            this.cellLifetime = cellLifetime;

            this.shapeLifetimeOriginState = shapeWeibullHSC;
            this.scaleLifetimeOriginState = scaleWeibullHSC;

            this.shapeLifetimeEndState = shapeWeibullMPP;
            this.scaleLifetimeEndState = scaleWeibullMPP;
        }

        /**
         * @param x the fraction of the branch length for which the cell-state transition happens
         * @return the probability density value of the cell transitioning at proportion 'x' of the branch length
         */
        @Override
        public double value(double x) {

            double timeUntilTransition = x * cellLifetime;
            double beforeTransition = Utils.getWeibullDensity(timeUntilTransition, scaleLifetimeOriginState, shapeLifetimeOriginState);
            double afterTransition = Utils.getConditionalWeibullDensity(cellLifetime, timeUntilTransition, scaleLifetimeEndState, shapeLifetimeEndState);

            return beforeTransition * afterTransition;
        }
    }

    private class DensityProbTransitionBeforeEndAndDivisionDeathAfter implements UnivariateFunction {

        double cellLifetime;

        double shapeLifetimeOriginState;
        double shapeLifetimeEndState1;
        double shapeLifetimeEndState2;

        double scaleLifetimeOriginState;
        double scaleLifetimeEndState1;
        double scaleLifetimeEndState2;

        double probaEndState1;
        double probaEndState2;
        double probaEndState3;


        DensityProbTransitionBeforeEndAndDivisionDeathAfter(double cellLifetime,
                                                            double shapeWeibullTrans, double scaleWeibullTrans,
                                                            double shapeWeibullDiv, double scaleWeibullDiv,
                                                            double shapeWeibullDeath, double scaleWeibullDeath,
                                                            double probaDiv, double probaDeath, double probaNothing) {
            this.cellLifetime = cellLifetime;

            this.shapeLifetimeOriginState = shapeWeibullTrans;
            this.scaleLifetimeOriginState = scaleWeibullTrans;

            this.shapeLifetimeEndState1 = shapeWeibullDiv;
            this.scaleLifetimeEndState1 = scaleWeibullDiv;

            this.shapeLifetimeEndState2 = shapeWeibullDeath;
            this.scaleLifetimeEndState2 = scaleWeibullDeath;

            this.probaEndState1 = probaDiv;
            this.probaEndState2 = probaDeath;
            this.probaEndState3 = probaNothing;
        }

        /**
         * @param x the fraction of the branch length for which the cell-state transition happens
         * @return the probability density value of the cell transitioning at proportion 'x' of the branch length
         */
        @Override
        public double value(double x) {

            double timeUntilTransition = x * cellLifetime;
            double beforeTransition = Utils.getWeibullDensity(timeUntilTransition, scaleLifetimeOriginState, shapeLifetimeOriginState);
            double afterTransition1 = probaEndState1 * Utils.getConditionalWeibullProbaAboveT(cellLifetime, timeUntilTransition, scaleLifetimeEndState1, shapeLifetimeEndState1);
            double afterTransition2 = probaEndState2 * Utils.getConditionalWeibullProbaAboveT(cellLifetime, timeUntilTransition, scaleLifetimeEndState2, shapeLifetimeEndState2);
            double afterTransition3 = probaEndState3; // HSC transitions and MPP does nothing

            return beforeTransition * (afterTransition1 + afterTransition2 + afterTransition3);
        }
    }

    public static void main(String[] args) throws IOException {
        String fileName = "../Data/Examples/MinusInfinityLik4.csv";
        LineageTree tree =  new LineageTree();
        tree.setInputValue("measuresCSVFile", fileName);

        tree.initAndValidate();

        System.out.println(tree.toString());

        LineageTreeProb probTree = new LineageTreeProb();
        probTree.setInputValue("tree", tree);

        List<RealParameter> fateProbabilities = new ArrayList<>();
        fateProbabilities.add(new RealParameter("1.0 0. 0."));
        fateProbabilities.add(new RealParameter("1.0 0. 0."));
        fateProbabilities.add(new RealParameter("1.0 0. 0."));

        List<RealParameter> scaleWeibull = new ArrayList<>();
        scaleWeibull.add(new RealParameter("10 10"));
        scaleWeibull.add(new RealParameter("10 10"));
        scaleWeibull.add(new RealParameter("10 10"));

        List<RealParameter> shapeWeibull = new ArrayList<>();
        shapeWeibull.add(new RealParameter("1 1"));
        shapeWeibull.add(new RealParameter("1 1"));
        shapeWeibull.add(new RealParameter("1 1"));

        List<RealParameter> transitionUponDivisionProbs = new ArrayList<>();
        transitionUponDivisionProbs.add(new RealParameter("0.3 0.1 0.1 0.15 0.15 0.2"));
        transitionUponDivisionProbs.add(new RealParameter("0 0 0 1.0 0 0 "));
        transitionUponDivisionProbs.add(new RealParameter("0 0 0 0 0 1.0 "));


        RealParameter meanLogNormalAreaGrowthRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter sdLogNormalAreaGrowthRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter zeroFractionAreaGrowthRateInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanNormalEccentricityInput = new RealParameter("0.5 .5 .5");
        RealParameter sdNormalEccentricityInput = new RealParameter("0.3 0.3 .3");

        RealParameter meanLogNormalInstantSpeedInput = new RealParameter("1. 1. 1.");
        RealParameter sdLogNormalInstantSpeedInput = new RealParameter("1.0 1.0 1.");
        RealParameter zeroFractionInstantSpeedInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanLogNormalTMRMProductionRateInput = new RealParameter("3.0 3.0 3.0");
        RealParameter sdLogNormalTMRMProductionRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter zeroFractionTMRMProductionRateInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanLogNormalTMRMMaxRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter sdLogNormalTMRMMaxRateInput = new RealParameter("1.0 1.0 1.0");
        RealParameter zeroFractionTMRMMaxRateInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanLogNormalROSProductionRateInput = new RealParameter("1.0 1.0 1.0");
        RealParameter sdLogNormalROSProductionRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter zeroFractionROSProductionRateInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanLogNormalCD71APCProductionRateInput = new RealParameter("5.0 5.0 5.0");
        RealParameter sdLogNormalCD71APCProductionRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter zeroFractionCD71APCProductionRateInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanLogNormalCD71PEProductionRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter sdLogNormalCD71PEProductionRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter zeroFractionCD71PEProductionRateInput = new RealParameter("0.1 0.1 0.1");

        RealParameter meanLogNormalcMycGFPMaxRateInput = new RealParameter("1.0 1.0 1.0");
        RealParameter sdLogNormalcMycGFPMaxRateInput = new RealParameter("2.0 2.0 2.0");
        RealParameter zeroFractioncMycGFPMaxRateInput = new RealParameter("0.1 0.1 0.1");


//        RealParameter meanNormalPerimeterGrowthRateInput = new RealParameter("0. 1.");
//        RealParameter sdNormalPerimeterGrowthRateInput = new RealParameter("5.0 10.0");

        probTree.setInputValue("transitionUponDivisionProbs", transitionUponDivisionProbs);
        probTree.setInputValue("fateProbabilities", fateProbabilities);
        probTree.setInputValue("scaleWeibull",scaleWeibull);
        probTree.setInputValue("shapeWeibull",shapeWeibull);
        probTree.setInputValue("lossProb", new RealParameter("0."));

//        probTree.setInputValue("rootIsHSC", new BooleanParameter("false"));
        probTree.setInputValue("rootType", new IntegerParameter("2"));

        probTree.setInputValue("meanAreaGrowthRate", meanLogNormalAreaGrowthRateInput);
        probTree.setInputValue("sdAreaGrowthRate", sdLogNormalAreaGrowthRateInput);
        probTree.setInputValue("zeroFractionAreaGrowthRate", zeroFractionAreaGrowthRateInput);

        probTree.setInputValue("meanEccentricity", meanNormalEccentricityInput);
        probTree.setInputValue("sdEccentricity", sdNormalEccentricityInput);

        probTree.setInputValue("meanInstantSpeed", meanLogNormalInstantSpeedInput);
        probTree.setInputValue("sdInstantSpeed", sdLogNormalInstantSpeedInput);
        probTree.setInputValue("zeroFractionInstantSpeed", zeroFractionInstantSpeedInput);

        probTree.setInputValue("meanTMRMProductionRate", meanLogNormalTMRMProductionRateInput);
        probTree.setInputValue("sdTMRMProductionRate", sdLogNormalTMRMProductionRateInput);
        probTree.setInputValue("zeroFractionTMRMProductionRate", zeroFractionROSProductionRateInput);

        probTree.setInputValue("meanTMRMMaxRate", meanLogNormalTMRMMaxRateInput);
        probTree.setInputValue("sdTMRMMaxRate", sdLogNormalTMRMMaxRateInput);
        probTree.setInputValue("zeroFractionTMRMMaxRate", zeroFractionTMRMMaxRateInput);

        probTree.setInputValue("meanROSProductionRate", meanLogNormalROSProductionRateInput);
        probTree.setInputValue("sdROSProductionRate", sdLogNormalROSProductionRateInput);
        probTree.setInputValue("zeroFractionROSProductionRate", zeroFractionROSProductionRateInput);

        probTree.setInputValue("meanCD71APCProductionRate", meanLogNormalCD71APCProductionRateInput);
        probTree.setInputValue("sdCD71APCProductionRate", sdLogNormalCD71APCProductionRateInput);
        probTree.setInputValue("zeroFractionCD71APCProductionRate", zeroFractionCD71APCProductionRateInput);

        probTree.setInputValue("meanCD71PEProductionRate", meanLogNormalCD71PEProductionRateInput);
        probTree.setInputValue("sdCD71PEProductionRate", sdLogNormalCD71PEProductionRateInput);
        probTree.setInputValue("zeroFractionCD71PEProductionRate", zeroFractionCD71PEProductionRateInput);

        probTree.setInputValue("meancMycGFPMaxRate", meanLogNormalcMycGFPMaxRateInput);
        probTree.setInputValue("sdcMycGFPMaxRate", sdLogNormalcMycGFPMaxRateInput);
        probTree.setInputValue("zeroFractioncMycGFPMaxRate", zeroFractioncMycGFPMaxRateInput);

//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);
//
//        probTree.setInputValue("",);
//        probTree.setInputValue("",);


        //probTree.setInputValue("meanPerimeterGrowthRate", meanNormalPerimeterGrowthRateInput);
//        probTree.setInputValue("sdPerimeterGrowthRate", sdNormalPerimeterGrowthRateInput);
//        probTree.setInputValue("meanCD41ProductionRate", meanNormalCD41ProductionRateInput);
//        probTree.setInputValue("sdCD41ProductionRate", sdNormalCD41ProductionRateInput);
//
//        probTree.setInputValue("meanFcgRIIIProductionRate", meanNormalFcgRIIIProductionRateInput);
//        probTree.setInputValue("sdFcgRIIIProductionRate", sdNormalFcgRIIIProductionRateInput);
        
        probTree.initAndValidate();
        double logP;
        logP = probTree.calculateLogP();

        System.out.println("logP = " + logP);

        System.out.println(Math.log(Double.MAX_VALUE));
    }


}



