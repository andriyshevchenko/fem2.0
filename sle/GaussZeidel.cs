using System;
using System.Collections.Generic;

namespace SystemOfEquations
{
    public class GaussSeidel
    {
        public double Epsilon { get; set; }

        public GaussSeidel(double[][] leftPart, double[] rightPart, double epsilon, bool isParallel)
        {
            if (Epsilon > 0.1)
                throw new ArgumentException("epsilon > 0.1", nameof(epsilon));

            Epsilon = epsilon;
            Answer = new double[rightPart.GetLength(0)];
            if (isParallel)
            {
                TrySolveParallel(leftPart, rightPart);
                return;
            }

            TrySolve(leftPart, rightPart);
        }
        public long Iterations { get; private set; }
        public double[] Answer { get; private set; }

        private static bool IsConverge(
            IReadOnlyList<double> curr,
            IReadOnlyList<double> prev, int n, double eps)
        {
            double norm = 0;
            for (var i = 0; i < n; i++)
            {
                norm += Math.Pow(curr[i] - prev[i], 2);
            }
            return !(Math.Sqrt(norm) >= eps);
        }

        private void TrySolve(double[][] leftPart, IReadOnlyList<double> rightPart)
        {
            var curr = new double[rightPart.Count];
            for (var i = 0; i < curr.Length; i++)
            {
                curr[i] = 0;
            }
            var prev = new double[rightPart.Count];

            var n = leftPart.GetLength(0);
            long counter = 0;
            do
            {
                if (counter > 1000000)
                    return;

                ++counter;

                for (var i = 0; i < n; i++)
                    prev[i] = curr[i];

                for (var i = 0; i < n; i++)
                {
                    double var = 0;
                    for (var j = 0; j < i; j++)
                        var += leftPart[i][j] * curr[j];
                    for (var j = i + 1; j < n; j++)
                        var += leftPart[i][j] * prev[j];
                    curr[i] = (rightPart[i] - var) / leftPart[i][i];
                }
            } while (!IsConverge(curr, prev, n, Epsilon));
            Answer = curr;
            Iterations = counter;
        }

        private void TrySolveParallel(double[][] leftPart, double[] rightPart)
        {
        }
    }
}