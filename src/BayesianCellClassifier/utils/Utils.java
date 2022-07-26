package BayesianCellClassifier.utils;

import BayesianCellClassifier.parametrization.DistributionForMeasurement;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Erf;

public class Utils {

    /**
     * Compute log density at x under a Weibull distribution with scale lambda
     * and shape parameter k.
     *
     * @param x parameter with distribution
     * @param lambda scale parameter
     * @param k shape parameter
     * @return log density at x
     */
    public static double getWeibullLogDensity(double x, double lambda, double k) {
        double loglambda = Math.log(lambda);
        double logk = Math.log(k);
        double logx = Math.log(x);

        return x<0.0
                ? Double.NEGATIVE_INFINITY
                : logk - loglambda + (k-1)*(logx - loglambda) - Math.pow(x/lambda, k);
    }

    /**
     * Compute density at x under a Weibull distribution with scale lambda
     * and shape parameter k.
     *
     * @param x parameter with distribution
     * @param lambda scale parameter
     * @param k shape parameter
     * @return density at x
     */
    public static double getWeibullDensity(double x, double lambda, double k) {

        return x<0.0
                ? 0.0
                : k/lambda * Math.pow(x/lambda, k-1) * Math.exp(-Math.pow(x/lambda, k));
    }

    /**
     * Compute density at x under a Weibull distribution with median m
     * and shape parameter k.
     *
     * @param x parameter with distribution
     * @param m median parameter
     * @param k shape parameter
     * @return density at x
     */
    public static double getWeibullDensityMedianShapeParam(double x, double m, double k) {

        if(x < 0.0 || m <= 0 || k <= 0)
            return 0.0;
        else{
            double lambda = m/Math.pow(Math.log(2), 1/k);
            return k/lambda * Math.pow(x/lambda, k-1) * Math.exp(-Math.pow(x/lambda, k));
        }
    }

    public static double getWeibullCumulativeDistribution(double x, double lambda, double k) {

        return x<0.0
                ? 0.0
                : 1 - Math.exp(-Math.pow(x/lambda, k));
    }

    public static double getWeibullCumulativeDistributionMedianShapeParam(double x, double m, double k) {

        if(x < 0.0 || m <= 0 || k <= 0)
            return 0.0;
        else{
            double lambda = m/Math.pow(Math.log(2), 1/k);
            return 1 - Math.exp(-Math.pow(x/lambda, k));
        }
    }



    /**
     * Compute log density at x under a Weibull distribution with scale lambda
     * and shape parameter k, conditioned on the lifetime x being higher than y.
     *
     * @param x parameter with distribution
     * @param lambda scale parameter
     * @param k shape parameter
     * @return density at x
     */
    public static double getConditionalWeibullLogDensity(double x, double y, double lambda, double k) {
        double logOneMinusCDF = -Math.pow(y/lambda, k);

        return (x<y || y<0)
                ? Double.NEGATIVE_INFINITY
                : getWeibullLogDensity(x, lambda, k) - logOneMinusCDF;
    }

    /**
     * Compute density at x under a Weibull distribution with scale lambda
     * and shape parameter k, conditioned on the lifetime x being higher than y.
     *
     * @param x parameter with distribution
     * @param lambda scale parameter
     * @param k shape parameter
     * @return density at x
     */
    public static double getConditionalWeibullDensity(double x, double y, double lambda, double k) {
        return (x<y || y<0)
                ? 0.0
                : getWeibullDensity(x, lambda, k)/Math.exp(-Math.pow(y/lambda, k));
    }

    /**
     * Compute probability that event happens after T under a Weibull distribution with scale lambda
     * and shape parameter k, conditioned on the lifetime being higher than y.
     *
     * @param T minimum lifetime
     * @param y conditional lifetime
     * @param lambda scale parameter
     * @param k shape parameter
     * @return conditional probability that event happens after time T
     */
    public static double getConditionalWeibullProbaAboveT(double T, double y, double lambda, double k) {
        return (T<y || y<0)
                ? 0.0
                : Math.exp(-Math.pow(T/lambda, k))/Math.exp(-Math.pow(y/lambda, k));
    }

    /**
     * Compute log of cumulative density function at x under a Weibull distribution with scale lambda
     * and shape parameter k.
     *
     * @param x parameter with distribution
     * @param lambda scale parameter
     * @param k shape parameter
     * @return density at x
     */
    public static double getWeibullLogCDF(double x, double lambda, double k) {

        return (x<0)
                ? Double.NEGATIVE_INFINITY
                : Math.log(1 - Math.exp(-Math.exp(k * Math.log(x/lambda))));
    }

    public static double getLogNormalDensity(double x, double mu, double sigma) {
        if(x <= 0)
            return 0;
        else
            return 1.0/ (x * sigma * Math.sqrt(2*Math.PI)) * Math.exp(-(Math.log(x) - mu)*(Math.log(x)-mu)/(2*sigma*sigma));
    }

    public static double getLogNormalCumulativeDistribution(double x, double mu, double sigma) {
        if(x <= 0)
            return 0;
        else
            return 0.5 * Erf.erfc((- Math.log(x) - mu)/(sigma * Math.sqrt(2))); // erfc is the complementary error function
    }

    public static double getGammaDensityShapeRateParam(double x, double alpha, double beta) {
        if (x < 0 || alpha <= 0 || beta <= 0)
            return 0;
        else
            return Math.pow(x * beta, alpha - 1) * beta * Math.exp(-x * beta) / Gamma.gamma(alpha);
    }

    public static double getGammaDensityShapeMeanParam(double x, double k, double mu) {
        if (x < 0 || k <= 0 || mu <= 0)
            return 0;
        else
            return Math.pow(x * k/mu, k - 1) * k/mu * Math.exp(-x * k/mu) / Gamma.gamma(k);
    }

    public static double getGammaDensityMeanShapeParam(double x, double mu, double k) {
        if (x < 0 || k <= 0 || mu <= 0)
            return 0;
        else
            return Math.pow(x * k/mu, k - 1) * k/mu * Math.exp(-x * k/mu) / Gamma.gamma(k);
    }

    public static double getGammaCumulativeDistributionShapeMeanParam(double x, double k, double mu) {
        if (x <= 0 || k <= 0 || mu <= 0)
            return 0;
        else
            return Gamma.regularizedGammaP(k,k * x/mu);
    }

    public static double getGammaCumulativeDistributionMeanShapeParam(double x, double mu, double k) {
        if (x <= 0 || k <= 0 || mu <= 0)
            return 0;
        else
            return Gamma.regularizedGammaP(k,k * x/mu);
    }


    public static double getNormalDensity(double x, double mu, double sigma) {
        return 1.0/Math.sqrt(2*Math.PI*sigma*sigma)*Math.exp(-(x - mu)*(x-mu)/(2*sigma*sigma));
    }

    public static double getNormalCumulativeDistribution(double x, double mu, double sigma) {
        return 0.5*(1 + Erf.erf((x - mu)/(sigma * Math.sqrt(2)))) ; // erf is the error function
    }

    public static double getNormalLogDensity(double x, double mu, double sigma) {
        return - 1.0/2 * Math.log(2*Math.PI*sigma*sigma) - (x - mu)*(x-mu)/(2*sigma*sigma);
    }

    public static double getBetaDensity(double x, double alpha, double beta) {

        if(x <0 || x>1 || alpha <= 0 || beta <= 0)
            return 0;
        else {
            double logBetaDensity = (alpha - 1) * Math.log(x) + (beta - 1) * Math.log(1 - x) - Beta.logBeta(alpha, beta);
            return Math.exp(logBetaDensity);
        }

    }

    public static double getBetaCumulativeDistribution(double x, double alpha, double beta) {
        if(x <0 || x>1 || alpha <= 0 || beta <= 0)
            return 0;
        else
            return Beta.regularizedBeta(x, alpha, beta);
    }

    public static double getRandomValueFromDistribution(DistributionForMeasurement.DistributionType distrType, double parm1, double parm2) {
        switch (distrType) {
            case NORMAL:
                return new NormalDistribution(parm1, parm2).sample();
            case LOGNORMAL:
                return new LogNormalDistribution(parm1, parm2).sample();
            case BETA:
                return new BetaDistribution(parm1, parm2).sample();
            case GAMMA_MEAN_SHAPE:
                if(parm2 == 0)
                    throw new IllegalArgumentException("Shape parameter of gamma distribution cannot be 0.");
                double shape = parm2;
                double scale = parm1/parm2;
                return new GammaDistribution(shape, scale).sample();
            default:
                throw new IllegalArgumentException("Distribution type not implemented yet in simulation.");
        }
    }


    public static void main(String[] args){
        double alpha = 80;
        double beta = 80;
        double gammaAlphaPlusBeta = Gamma.gamma(alpha + beta);

        double betaValue = Beta.logBeta(alpha, beta);

        System.out.println(betaValue);

        System.out.println(gammaAlphaPlusBeta);
    }

}
