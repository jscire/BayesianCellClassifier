package generalclassifier.lineagetree;

public class Cell {

    public enum Fate {
        // D: divides, A: apoptoses, L: lost cell, N: does nothing (not sure we'll keep that), U: unobserved
        D, A, L, N, U
    }

    Fate fate;

    double lifetime;
    double area;
    double eccentricity;
    double perimeter;
    double averageSpeed;



    public Cell(){
    }



    public Cell(double branchlength) {
        this.lifetime = branchlength;
    }

    public double getArea() {
        return area;
    }

    public void setArea(double area) {
        this.area = area;
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

    public double getEccentricity() {
        return eccentricity;
    }

    public void setEccentricity(double eccentricity) {
        this.eccentricity = eccentricity;
    }

    public double getAverageSpeed() {
        return averageSpeed;
    }

    public void setAverageSpeed(double averageSpeed) {
        this.averageSpeed = averageSpeed;
    }

    public double getPerimeter() {
        return perimeter;
    }

    public void setPerimeter(double perimeter) {
        this.perimeter = perimeter;
    }
}
