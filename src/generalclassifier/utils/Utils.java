package generalclassifier.utils;

import generalclassifier.lineagetree.MeasureType;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

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


    public static double getNormalDensity(double x, double mu, double sigma) {
        return 1.0/Math.sqrt(2*Math.PI*sigma*sigma)*Math.exp(-(x - mu)*(x-mu)/(2*sigma*sigma));
    }

    public static double getTruncatedNormalDensity(double x, int nodeType, InputPair input, MeasureType measureType) {

        NormalDistribution dist = new NormalDistributionImpl(input.getMean().get().getArrayValue(nodeType),
                                                            input.getStandardDev().get().getArrayValue(nodeType));

        double upperBound = measureType.getHardUpperBound();
        double lowerBound = measureType.getHardLowerBound();

        if(lowerBound == Double.NEGATIVE_INFINITY && upperBound == Double.POSITIVE_INFINITY) return dist.density(x);

        double scaleFactor = 1.0;

        try{
            if(upperBound == Double.POSITIVE_INFINITY) {
                scaleFactor = 1.0 / (1.0 - dist.cumulativeProbability(lowerBound));
            }
            else if(lowerBound == Double.NEGATIVE_INFINITY) {
                scaleFactor = 1.0 / dist.cumulativeProbability(upperBound);
            }
            else {
                scaleFactor = 1.0/(dist.cumulativeProbability(upperBound) - dist.cumulativeProbability(lowerBound));
            }
        }
        catch(Exception e) {
            System.out.println("Cumulative probability of normal distribution crashed." +
                    "\nLower bound" + lowerBound +
                    "\nUpper bound" + upperBound +
                    "\nx: " + x +
                    "\nmean: " + dist.getMean() +
                    "\nsd: " + dist.getStandardDeviation());
        }

        return scaleFactor * dist.density(x);
    }







}
