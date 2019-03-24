using System;

namespace FEM
{
    public partial class SmallTransverseStringVibrations : FEMTwoDimBase
    {
        public struct BoundaryCondition
        {
            public Func<double, double> U0 { get; }
            public Func<double, double> V0 { get; }
            public BoundaryCondition(Func<double, double> u0, Func<double, double> v0)
            {
                U0 = u0;
                V0 = v0;
            }
        }
    }
}
