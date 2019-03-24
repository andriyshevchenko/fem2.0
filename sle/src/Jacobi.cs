using System;

namespace SystemOfEquations
{
    public class Jacobi : IterationMethod
    {
        public Jacobi(double[][] leftPart, double[] rightPart, double epsilon, bool isParallel)
            : base(leftPart, rightPart, epsilon, isParallel)
        {
        }

        protected override void TrySolve(double[][] leftPart, double[] rightPart)
        {
            var prev = new double[rightPart.Length];
            for (var i = 0; i < N; i++)
            {
                prev[i] = 0;
            }

            long counter = 0;
            while (true)
            {
                if (counter > 10000)
                    //    return;
                    ++counter;

                var curr = new double[rightPart.Length];

                for (var i = 0; i < N; i++)
                {
                    curr[i] = rightPart[i];
                    for (var j = 0; j < N; j++)
                    {
                        if (i == j) continue;
                        curr[i] -= leftPart[i][j]*prev[j];
                    }
                    curr[i] /= leftPart[i][i];
                }

                double error = 0;

                for (var i = 0; i < N; i++)
                {
                    error += Math.Abs(curr[i] - prev[i]);
                }

                if (error < Epsilon)
                    break;

                prev = curr;
            }

            Answer = prev;
            Iterations = counter;
        }

        protected override void TrySolveParallel(double[][] leftPart, double[] rightPart)
        {
        }
    }
}