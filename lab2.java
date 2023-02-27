package com.lab;

import Jama.LUDecomposition;
import Jama.Matrix;

public class Main {

    static final int N = 8;

    static double[][] signal = {
            {2.0, 3.0, 4.0, 0.0, 4.0, 5.0, 3.0, 3.0},
            {1.0, 2.0, 0.0, 3.0, 3.0, 0.0, 6.0, 5.0},
            {2.0, 3.0, 4.0, 0.0, 5.0, 4.0, 5.0, 6.0},
            {3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0},
            {4.0, 4.0, 0.0, 2.0, 5.0, 5.0, 3.0, 3.0},
            {2.0, 2.0, 0.0, 3.0, 4.0, 3.0, 5.0, 7.0},
            {0.0, 4.0, 4.0, 4.0, 0.0, 4.0, 0.0, 5.0},
            {3.0, 2.0, 0.0, 5.0, 5.0, 7.0, 0.0, 3.0}
    };

    static double[][] noisySignal1 = {
            {2.0, 3.0, 4.0, 0.0, 4.0, 5.0, 3.0, 3.0},
            {1.0, 2.0, 0.0, 3.0, 3.0, 0.0, 6.0, 5.0},
            {2.0, 3.0, 4.0, 0.0, 5.0, 4.0, 5.0, 6.0},
            {3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0},
            {4.0, 4.0, 0.0, 2.0, 5.0, 5.0, 3.0, 3.0},
            {2.0, 2.0, 0.0, 3.0, 4.0, 3.0, 9.0, 7.0},
            {0.0, 4.0, 4.0, 4.0, 0.0, 4.0, 0.0, 5.0},
            {3.0, 2.0, 0.0, 5.0, 5.0, 7.0, 0.0, 3.0}
    };

    static double[][] noisySignal2 = {
            {2.0, 3.0, 4.0, 0.0, 4.0, 5.0, 3.0, 3.0},
            {1.0, 2.0, 0.0, 3.0, 3.0, 0.0, 6.0, 5.0},
            {2.0, 3.0, 9.0, 0.0, 5.0, 4.0, 5.0, 6.0},
            {3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0},
            {4.0, 4.0, 0.0, 2.0, 5.0, 5.0, 3.0, 3.0},
            {2.0, 2.0, 0.0, 3.0, 4.0, 3.0, 5.0, 7.0},
            {0.0, 4.0, 4.0, 4.0, 0.0, 4.0, 0.0, 5.0},
            {3.0, 2.0, 0.0, 5.0, 5.0, 7.0, 0.0, 3.0}
    };

    public static void main(String[] args) {
        System.out.println("signal:");
        printMatrix(signal, 0);

        System.out.println("noisySignal1:");
        printMatrix(noisySignal1, 0);

        System.out.println("noisySignal2:");
        printMatrix(noisySignal2, 0);

        double[][] spectorOfSignal = computeBPA2N(signal);
        spectorOfSignal = computeNormalizeMatrix(spectorOfSignal);
        System.out.println("spectorOfSignal:");
        printMatrix(spectorOfSignal, 3);

        double[][] spectorOfNoisySignal1 = computeBPA2N(noisySignal1);
        spectorOfNoisySignal1 = computeNormalizeMatrix(spectorOfNoisySignal1);
        System.out.println("spectorOfNoisySignal1:");
        printMatrix(spectorOfNoisySignal1, 3);

        double[][] spectorOfNoisySignal2 = computeBPA2N(noisySignal2);
        spectorOfNoisySignal2 = computeNormalizeMatrix(spectorOfNoisySignal2);
        System.out.println("spectorOfNoisySignal2:");
        printMatrix(spectorOfNoisySignal2, 3);

        double[][] filterMatrix = computeFilter(spectorOfSignal,
                spectorOfNoisySignal1);
        System.out.println("filterMatrix:");
        printMatrix(filterMatrix, 3);

        double[][] mult1 =computeMultiplication(
                spectorOfNoisySignal1,filterMatrix);
        System.out.println("Multiplication spectorOfNoisySignal2 matrix on filter");
        printMatrix(mult1,3);
        double[][] mult2 =computeMultiplication(spectorOfNoisySignal2,filterMatrix);
        System.out.println("Multiplication spectorOfNoisySignal3 matrix on filter");
        printMatrix(mult2,3);

        double[][] filteredNoisySignal1 = computeReversBPA2N(mult1);
        System.out.println("filteredNoisySignal1:");
        printMatrix(filteredNoisySignal1, 1);


        double[][] filteredNoisySignal2 = computeReversBPA2N(mult2);
        System.out.println("filteredNoisySignal2:");
        printMatrix(filteredNoisySignal2, 1);
    }

    static double[][] computeBPA2N(double[][] matrix) {
        matrix = matrix.clone();
        for (int i = 0; i < N; i++) {
            matrix[i] = computeBPA(matrix[i]);
        }
        for (int i = 0; i < N; i++) {
            double[] result = computeBPA(new double[]{matrix[0][i], matrix[1][i],
                    matrix[2][i], matrix[3][i],
                    matrix[4][i], matrix[5][i], matrix[6][i], matrix[7][i]});
            for (int j = 0; j < N; j++) {
                matrix[j][i] = result[j];
            }
        }
        return matrix;
    }

    static double[] computeBPA(double[] line) {
        line = new double[]{line[0], line[4], line[2], line[6], line[1], line[5],
                line[3], line[7]};
        double[] result = new double[N];
        for (int i = 0; i < 4; i++) {
            result[i] = line[i] + line[i + 4];
        }
        for (int i = 4; i < 8; i++) {
            result[i] = -line[i] + line[i - 4];
        }
        line = result.clone();
        for (int i = 0; i < 2; i++) {
            result[i] = line[i] + line[i + 2];
        }
        for (int i = 2; i < 4; i++) {
            result[i] = -line[i] + line[i - 2];
        }
        for (int i = 4; i < 6; i++) {
            result[i] = line[i] - line[i + 2];
        }
        for (int i = 6; i < 8; i++) {
            result[i] = line[i] + line[i - 2];
        }
        line = result.clone();
        result[0] = line[0] + line[1];
        result[1] = line[0] - line[1];
        result[2] = line[2] - line[3];
        result[3] = line[2] + line[3];
        result[4] = line[4] + line[5];
        result[5] = line[4] - line[5];
        result[6] = line[6] - line[7];
        result[7] = line[6] + line[7];
        return result;
    }

    static double[][] computeReversBPA2N(double[][] matrix) {
        matrix = matrix.clone();
        for (int i = 0; i < N; i++) {
            matrix[i] = computeReverseBPA(matrix[i]);
        }
        for (int i = 0; i < N; i++) {
            double[] result = computeReverseBPA(new double[]{matrix[0][i], matrix[1][i],
                    matrix[2][i], matrix[3][i], matrix[4][i], matrix[5][i],
                    matrix[6][i], matrix[7][i]});
            for (int j = 0; j < N; j++) {
                matrix[j][i] = result[j];
            }
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                matrix[i][j] = matrix[i][j];
            }
        }
        return matrix;
    }

    static double[] computeReverseBPA(double[] line) {
        line = line.clone();
        double[] result = new double[N];
        for (int i = 0; i < 4; i++) {
            result[i] = line[i] + line[i + 4];
        }
        for (int i = 4; i < 8; i++) {
            result[i] = -line[i] + line[i - 4];
        }
        line = result.clone();
        for (int i = 0; i < 2; i++) {
            result[i] = line[i] + line[i + 2];
        }
        for (int i = 2; i < 4; i++) {
            result[i] = -line[i] + line[i - 2];
        }
        for (int i = 4; i < 6; i++) {
            result[i] = line[i] - line[i + 2];
        }
        for (int i = 6; i < 8; i++) {
            result[i] = line[i] + line[i - 2];
        }
        line = result.clone();
        result[0] = line[0] + line[1];
        result[1] = line[0] - line[1];
        result[2] = line[2] - line[3];
        result[3] = line[2] + line[3];
        result[4] = line[4] + line[5];
        result[5] = line[4] - line[5];
        result[6] = line[6] - line[7];
        result[7] = line[6] + line[7];
        return new double[]{result[0], result[7], result[4], result[3],
                result[2], result[5], result[6], result[1]};
    }

    static double [] [] computeFilter(double[][] signal, double[] [] noisySignal) {
        LUDecomposition luDecomposition = new LUDecomposition(
                new Matrix(noisySignal));
        Matrix result = luDecomposition.solve(new Matrix(signal));
        return result.getArray();
    }

    static public double[][] computeNormalizeMatrix(double[][] matrix) {
        matrix = matrix.clone();
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                matrix[i][j] /= N * N;
            }
        }
        return matrix;
    }

    static double[][] computeMultiplication(double[][] lval, double[][] rval) {
        double[][] result = new double[lval.length][rval[0].length];
        for (int i = 0; i < lval.length; i++) {
            for (int j = 0; j < rval[0].length; j++) {
                for (int k = 0; k < rval.length; k++) {
                    result[i][j] += lval[i][k] * rval[k][j];
                }
            }
        }
        return result;
    }

    static void printMatrix(double[][] matrix, int places) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                System.out.printf("%." + places + "f", matrix[i][j]);
                if (j != N - 1) {
                    System.out.print("\t");
                }
            }
            System.out.print("\n");
        }
        System.out.print("\n");
    }
}
