using System;

namespace FEM
{
    public partial class DiffusionConvectionReactionTwoDim :FEMTwoDimBase
    {
        public struct BoundaryCondition
        {
            public Func<double, double> U0 { get; }
            public Func<double, double> U1 { get; }
            public Func<double, double> U_0 { get; }
            public BoundaryCondition(
                Func<double, double> u0,
                Func<double, double> u1,
                Func<double, double> u_0)
            {
                U0 = u0;
                U1 = u1;
                U_0 = u_0;
            }
        }
    }
}
