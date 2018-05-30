package generalclassifier.utils;

import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.UnivariateRealFunction;

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






}
