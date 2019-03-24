using System;
using System.Linq;

namespace SystemOfEquations
{
    internal class Program
    {
        public static double[][] _matrix;
        public static double[] _row;

        private static void Main()
        {
            Console.WriteLine("Epsilon is");
            var eps = Convert.ToDouble("0,0001");

            var random = new Random();
            const int n = 1000;

           // _matrix = MatrixHelper.Matrix.GetRandom(n, n, () => random.Next(100));
            SolveMethod parallel = new LUDecomposition(_matrix, _row, true);
            SolveMethod serial = new LUDecomposition(_matrix, _row, false);

            var time = MatrixHelper.Matrix.GetExecutionTime(parallel.Solve);
            var timeSerial = 0; // GetExecutionTime(serial.Solve);

            Console.WriteLine(
                "{0} error\ntime {1} ms",
                MatrixHelper.Matrix.GetError(_matrix, _row, parallel.Answer).Aggregate(0.0, (d, _d) => d + _d), time);

            Console.WriteLine(
                "{0} error\ntime {1} ms",
                MatrixHelper.Matrix.GetError(_matrix, _row, serial.Answer).Aggregate(0.0, (d, _d) => d + _d), timeSerial);
            Console.ReadLine();
        }
    }
}