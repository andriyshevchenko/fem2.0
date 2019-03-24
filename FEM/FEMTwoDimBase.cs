using static System.Collections.Generic.Create;
using static System.Collections.Generic.SeriesCreate;
using System.Collections.Generic;

namespace FEM
{
    public abstract class FEMTwoDimBase : FEMBase
    {
        public int Nt;
        public double TMax;
        public double[] T;

        protected double[] U0;
        protected double[][] UiCalculated;

        public FEMTwoDimBase(double tMax, int nT, double[] elements) 
            : base(elements)
        {
            Nt = nT;
            TMax = tMax;
            T = Series(0.0, TMax, Nt + 1);
            UiCalculated = array<double[]>(Nt);
        }

        public double GetApproximateTime(int tIndex)
        {
            return T[tIndex] + 0.5 * Step(tIndex);
        }

        public IReadOnlyList<double> Time => T;

        protected double Step(int i)
        {
            return T[i + 1] - T[i];
        }
    }
}
