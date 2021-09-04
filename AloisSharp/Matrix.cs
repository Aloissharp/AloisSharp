using System;

namespace AloisSharp
{
    public static class Matrix
    {
        //Inverse of a matrix
        public static double[][] inverse(double[][] matrix)
        {
            double[][] inverse = new double[matrix.Length][];

            for (int j = 0; j < matrix.Length; j++)
                inverse[j] = new double[matrix[j].Length];

            // minors and cofactors
            for (int i = 0; i < matrix.Length; i++)
                for (int j = 0; j < matrix[i].Length; j++)
                    inverse[i][j] = Math.Pow(-1.0, i + j) * determinant(minor(matrix, i, j));

            // adjugate and determinant
            double det = 1.0 / determinant(matrix);

            for (int i = 0; i < inverse.Length; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    double temp = inverse[i][j];
                    inverse[i][j] = inverse[j][i] * det;
                    inverse[j][i] = temp * det;
                }
            }

            return inverse;
        }

        //Determinant of a matrix
        public static double determinant(double[][] matrix)
        {

            if (matrix.Length == 2)
                return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

            double det = 0;
            for (int i = 0; i < matrix[0].Length; i++)
                det += Math.Pow(-1.0, i) * matrix[0][i] * determinant(minor(matrix, 0, i));


            return det;
        }

        //Multiply matrix m1 and m2
        public static double[][] Multiply(double[][] m1, double[][] m2)
        {

            int i, j, k, r1, c1, r2, c2;

            double sum = 0;

            r1 = m1.Length;
            c1 = m1[0].Length;

            r2 = m2.Length;
            c2 = m2[0].Length;


            double[][] m3 = new double[r1][];

            for (i = 0; i < r1; i++)
                m3[i] = new double[c2];

            for (i = 0; i < r1; i++)
                for (j = 0; j < c2; j++)
                    m3[i][j] = 0;

            for (i = 0; i < r1; i++)    //row of first matrix
            {
                for (j = 0; j < c2; j++)    //column of second matrix
                {
                    sum = 0;
                    for (k = 0; k < c1; k++)
                        sum = sum + m1[i][k] * m2[k][j];
                    m3[i][j] = sum;
                }
            }

            return m3;

        }

        //Multiply matrix and array
        public static double[] Multiply(double[][] mat, double[] arr)
        {

            int i, j, r, c;

            double sum = 0.0;

            r = mat.Length;
            c = mat[0].Length;

            double[] arm = new double[r];

            for (i = 0; i < r; i++)
                arm[i] = 0.0;

            for (i = 0; i < c; i++)
            {
                sum = 0.0;

                for (j = 0; j < r; j++)
                    sum = sum + mat[i][j] * arr[j];

                arm[i] = sum;

            }

            return arm;

        }

        //Traspose of a matrix
        public static double[][] Traspose(double[][] m)
        {
            int i, j;


            int r = m.Length;
            int c = m[0].Length;


            double[][] mt = new double[c][];

            for (i = 0; i < c; i++)
                mt[i] = new double[r];


            for (i = 0; i < c; i++)
                for (j = 0; j < r; j++)
                    mt[i][j] = 0;


            for (i = 0; i < r; i++)
            {
                for (j = 0; j < c; j++)
                {

                    mt[j][i] = m[i][j];
                }
            }

            return mt;

        }

        //Submatrix
        public static double[][] minor(double[][] matrix, int row, int column)
        {
            double[][] minor = new double[matrix.Length - 1][];

            for (int j = 0; j < matrix.Length - 1; j++)
                minor[j] = new double[matrix[j].Length - 1];


            for (int i = 0; i < matrix.Length; i++)
                for (int j = 0; i != row && j < matrix[i].Length; j++)
                    if (j != column)
                        minor[i < row ? i : i - 1][j < column ? j : j - 1] = matrix[i][j];
            return minor;
        }

        //Duplicate of a matrix
        public static double[][] MatrixDuplicate(double[][] matrix)
        {
            // assumes matrix is not null.
            double[][] result = new double[matrix.Length][];

            for (int i = 0; i < matrix.Length; ++i)
            {
                result[i] = new double[matrix.Length];

                for (int j = 0; j < matrix[i].Length; ++j)
                    result[i][j] = matrix[i][j];
            }

            return result;
        }

        //Duplicate of an array
        public static double[] ArrayDuplicate(double[] array)
        {

            double[] result = new double[array.Length];

            for (int i = 0; i < array.Length; ++i)
            {
                result[i] = array[i];

            }

            return result;
        }

        //Switch rows
        public static void SwitchRows(int n, double[] y, double[][] m)
        {
            double tempD;
            int i, j;

            int maxOrder = m.Length;

            for (i = n; i <= maxOrder - 2; i++)
            {
                for (j = 0; j <= maxOrder - 1; j++)
                {
                    tempD = m[i][j];
                    m[i][j] = m[i + 1][j];
                    m[i + 1][j] = tempD;
                }
                tempD = y[i];
                y[i] = y[i + 1];
                y[i + 1] = tempD;
            }
        }

    }
}
