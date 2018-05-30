package generalclassifier.utils;

import java.io.ByteArrayOutputStream;
import java.util.BitSet;

public class MiscUtils {
    public static boolean contains(int[] array, int element){
        for (int i = 0; i < array.length; i++) {
            if(array[i] == element) return true;
        }
        return false;
    }

    public static int match(int[] array, int element){
        for (int i = 0; i < array.length; i++) {
            if(array[i] == element) return i;
        }
        return -1;
    }


    public static String toString(BitSet bitSet) {
        String interResult = Long.toString(bitSet.toLongArray()[0], 2);
        String result = new StringBuilder(interResult).reverse().toString(); // reverse the string to get bit number 0 on the left and the highest index bit on the right end.
        return result;
    }

}
