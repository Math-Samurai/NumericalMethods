using System;

static class QRFactorization
{
    public static double[,] GramSchmidt(ref double[,] Q, double[,] matrix)
    {
        double[,] R = new double[matrix.GetLength(1), matrix.GetLength(1)];
        for (int i = 0; i < matrix.GetLength(1); ++i)
        {
            double[] vec = new double[matrix.GetLength(0)];
            for (int k = 0; k < matrix.GetLength(0); ++k)
            {
                vec[k] = matrix[k, i];
            }
            double[] t = new double[matrix.GetLength(0)];
            vec.CopyTo(t, 0);
            for (int j = 0; j < i; ++j)
            {
                double temp = 0;
                for (int k = 0; k < matrix.GetLength(0); ++k)
                {
                    temp += vec[k] * Q[k, j];
                }
                R[j, i] = temp;
                for (int k = 0; k < matrix.GetLength(0); ++k)
                {
                    t[k] -= temp * Q[k, j];
                }
            }
            t.CopyTo(vec, 0);
            double norm = EuclideanNorm(vec);
            R[i, i] = norm;
            for (int k = 0; k < matrix.GetLength(0); ++k)
            {
                vec[k] /= norm;
            }
            for (int k = 0; k < matrix.GetLength(0); ++k)
            {
                Q[k, i] = vec[k];
            }
        }
        return R;
    }
    public static double[,] ModifiedGramSchmidt(ref double[,] Q, double[,] matrix)
    {
        double[,] R = new double[matrix.GetLength(1), matrix.GetLength(1)];
        double[,] V = new double[matrix.GetLength(0), matrix.GetLength(1)];
        for (int i = 0; i < matrix.GetLength(0); ++i)
        {
            for (int j = 0; j < matrix.GetLength(1); ++j)
            {
                V[i, j] = matrix[i, j];
            }
        }
        for (int i = 0; i < matrix.GetLength(1); ++i)
        {
            double temp = 0;
            for (int k = 0; k < matrix.GetLength(0); ++k)
            {
                temp += V[k, i] * V[k, i];
            }
            R[i, i] = Math.Sqrt(temp);
            for (int k = 0; k < matrix.GetLength(0); ++k)
            {
                Q[k, i] = V[k, i] / R[i, i];
            }
            for (int j = i + 1; j < matrix.GetLength(1); ++j)
            {
                R[i, j] = 0;
                for (int k = 0; k < matrix.GetLength(0); ++k)
                {
                    R[i, j] += Q[k, i] * V[k, j];
                }
                for(int k = 0; k < matrix.GetLength(0); ++k)
                {
                    V[k, j] -= R[i, j] * Q[k, i];
                }
            }
        }
        return R;
    }
    public static double EuclideanNorm(double[] vec)
    {
        double temp = 0;
        for (int i = 0; i < vec.Length; ++i)
        {
            temp += vec[i] * vec[i];
        }
        return Math.Sqrt(temp); ;
    }
}

static class SLAE
{
    public static double[] GaussPivot(double[,] matrix, double[] right)
    {
        if (matrix.GetLength(0) != matrix.GetLength(1) || matrix.GetLength(0) != right.Length)
        {
            throw new DivideByZeroException();
        }
        int length = matrix.GetLength(0);
        double[] res = new double[length];
        int[] permutation = new int[length - 1];
        for (int i = 0; i < length - 1; ++i)
        {
            permutation[i] = FindMax(matrix, i);
            SwapColumns(matrix, i, permutation[i]);
            for (int j = i + 1; j < length; ++j)
            {
                double divider = matrix[j, i] / matrix[i, i];
                for (int k = i; k < length; ++k)
                {
                    matrix[j, k] -= matrix[i, k] * divider;
                }
                right[j] -= right[i] * divider;
            }
        }
        for (int i = length - 1; i >= 0; --i)
        {
            double temp = 0;
            for (int j = length - 1; j > i; --j)
            {
                temp += matrix[i, j] * res[j];
            }
            res[i] = (right[i] - temp) / matrix[i, i];
        }
        for (int i = length - 2; i >= 0; --i)
        {
            double temp = res[i];
            res[i] = res[permutation[i]];
            res[permutation[i]] = temp;
        }
        return res;
    }

    private static int FindMax(double[,] matrix, int StrNum)
    {
        double max = Math.Abs(matrix[StrNum, StrNum]);
        int maxNum = StrNum;
        for (int i = StrNum + 1; i < matrix.GetLength(0); ++i)
        {
            if (Math.Abs(matrix[StrNum, i]).CompareTo(Math.Abs(max)) > 0)
            {
                max = matrix[StrNum, i];
                maxNum = i;
            }
        }
        return maxNum;
    }
    private static void SwapColumns(double[,] matrix, int FirstCol, int SecondCol)
    {
        if (FirstCol == SecondCol)
        {
            return;
        }
        for (int i = 0; i < matrix.GetLength(0); ++i)
        {
            double temp = matrix[i, FirstCol];
            matrix[i, FirstCol] = matrix[i, SecondCol];
            matrix[i, SecondCol] = temp;
        }
    }
}

static class Operations
{
    public static double[] multiplication(double[,] mx, double[] vec)
    {
        double[] res = new double[mx.GetLength(0)];
        for (int i = 0; i < mx.GetLength(0); ++i)
        {
            double temp = 0;
            for (int j = 0; j < mx.GetLength(1); ++j)
            {
                temp += mx[i, j] * vec[j];
            }
            res[i] = temp;
        }

        return res;
    }
    public static double[,] transpose(double[,] mx)
    {
        double[,] res = new double[mx.GetLength(1), mx.GetLength(0)];
        for (int i = 0; i < mx.GetLength(0); ++i)
        {
            for (int j = 0; j < mx.GetLength(1); ++j)
            {
                res[j, i] = mx[i, j];
            }
        }
        return res;
    }
} 

class Demo
{
    static int Main()
    {
        int PointCount = 0;
        int PolynomDegree = 0;
        double step;
        Console.WriteLine("Введите количество точек: ");
        PointCount = Convert.ToInt32(Console.ReadLine());
        Console.WriteLine("Введите степень многочлена");
        PolynomDegree = Convert.ToInt32(Console.ReadLine());
        step = 1.0f / PointCount;
        double[,] matrix = new double[PointCount, PolynomDegree];
        double[] function = new double[PointCount];
        for (int i = 0; i < PointCount; ++i)
        {
            for (int j = 0; j < PolynomDegree; ++j)
            {
                matrix[i, j] = Math.Pow(i * step, j);
            }
            function[i] = Math.Cos(4 * i * step);
        }
        double[,] GSQ = new double[PointCount, PolynomDegree];
        double[,] GSMQ = new double[PointCount, PolynomDegree];
        double[,] GSR = QRFactorization.GramSchmidt(ref GSQ, matrix);
        double[,] GSMR = QRFactorization.ModifiedGramSchmidt(ref GSMQ, matrix);
        double[] GSNewFunction = new double[PointCount];
        double[] GSMNewFunction = new double[PointCount];
        function.CopyTo(GSNewFunction, 0);
        function.CopyTo(GSMNewFunction, 0);
        GSNewFunction = Operations.multiplication(Operations.transpose(GSQ), GSNewFunction);
        GSMNewFunction = Operations.multiplication(Operations.transpose(GSMQ), GSMNewFunction);
        double[] GSRes = SLAE.GaussPivot(GSR, GSNewFunction);
        double[] GSMRes = SLAE.GaussPivot(GSMR, GSMNewFunction);
        double[] t1 = Operations.multiplication(matrix, GSRes);
        double[] t2 = Operations.multiplication(matrix, GSMRes);
        for(int i = 0; i < t1.Length; ++i)
        {
            t1[i] -= function[i];
            t2[i] -= function[i];
        }
        Console.WriteLine(QRFactorization.EuclideanNorm(t1));
        Console.WriteLine(QRFactorization.EuclideanNorm(t2));
        Console.ReadLine();
        return 0;
    }
}