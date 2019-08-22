package generalclassifier.utils;

public class Pair {

    int firstInt;
    int secondInt;

    public Pair(int i1, int i2){
        firstInt = i1;
        secondInt = i2;
    }

    public int[] getPairAsArray(){
        return new int[]{firstInt, secondInt};
    }

    public int getFirstInt(){
        return firstInt;
    }

    public int getSecondInt(){
        return secondInt;
    }

    public void setFirstInt(int firstInt) {
        this.firstInt = firstInt;
    }

    public void setSecondInt(int secondInt) {
        this.secondInt = secondInt;
    }
}
