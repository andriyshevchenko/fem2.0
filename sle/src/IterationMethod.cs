using System;

namespace SystemOfEquations
{
    public abstract class IterationMethod : SolveMethod
    {
        protected IterationMethod(double[][] leftPart, double[] rightPart, double epsilon, bool isParallel)
            : base(leftPart, rightPart, isParallel)
        {
            if (Epsilon > 0.1)
                throw new ArgumentOutOfRangeException(nameof(epsilon), "epsilon > 0.1");
            Epsilon = epsilon;
        }

        public double Epsilon { get; set; }
        public long Iterations { get; protected set; }

        public bool Converge => !Answer[0].Equals(double.NaN)
                              && !double.IsInfinity(Answer[0]);
    }
}