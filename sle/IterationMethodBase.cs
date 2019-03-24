using System;
using System.Collections.Generic;

namespace SystemOfEquations
{
    public abstract class IterationMethodBase
    {
        protected IterationMethodBase(double[][] leftPart, double[] rightPart, bool isParallel, double epsilon)
        {
            LeftPart = leftPart;
            RightPart = rightPart;
            IsParallel = isParallel;

            if (Epsilon > 0.1)
                throw new ArgumentException("epsilon > 0.1", nameof(epsilon));

            Answer = new double[rightPart.GetLength(0)];
        }
        public double Epsilon { get; set; }
        public long Iterations { get; protected set; }
        public double[] Answer { get; protected set; }

        public bool Converge => !Answer[0].Equals(double.NaN)
                                && !double.IsInfinity(Answer[0]);

        protected double[][] LeftPart { get; }
        protected double[] RightPart { get; }
        protected bool IsParallel { get; }

        protected abstract void TrySolve(double[][] leftPart, IReadOnlyList<double> rightPart);
        protected abstract void TrySolveParallel(double[][] leftPart, double[] rightPart);

        public void Solve()
        {
            if (IsParallel)
            {
                TrySolveParallel(LeftPart, RightPart);
                return;
            }
            TrySolve(LeftPart, RightPart);
        }
    }
}