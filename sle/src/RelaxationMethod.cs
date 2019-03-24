using System;

namespace SystemOfEquations
{
    public class RelaxationMethod : IterationMethod
    {
        public RelaxationMethod(double[][] leftPart, double[] rightPart, double omega, double epsilon, bool isParallel)
            : base(leftPart, rightPart, epsilon, isParallel)
        {
            Omega = omega;
        }

        public double Omega { get; set; }

        protected override void TrySolve(double[][] leftPart, double[] rightPart)
        {
            var x = new double[N];
            var xn = new double[N];
            double norm = 0;
            long counter = 0;
            do
            {
                ++counter;
                if (counter > 5000
                    || double.IsNaN(x[0])
                    || double.IsInfinity(x[0]))
                {
                    Answer[0] = double.NaN;
                    return;
                }

                for (var i = 0; i < N; i++)
                {
                    x[i] = rightPart[i];
                    for (var j = 0; j < N; j++)
                    {
                        if (i != j)
                            x[i] -= leftPart[i][j]*x[j];
                    }
                    x[i] /= leftPart[i][i];
                    x[i] = Omega*x[i] + (1 - Omega)*xn[i];

                    var error = Math.Abs(x[i] - xn[i]);
                    if (error > norm)
                        norm = error;

                    xn[i] = x[i];
                }
            } while (norm > Epsilon);
            Answer = x;
        }

        protected override void TrySolveParallel(double[][] leftPart, double[] rightPart)
        {
        }
    }
}