namespace FEM
{
    public partial class DiffusionConvectionReaction : FEMBase
    {
        public struct BoundaryCondition
        {
            public BoundaryCondition(double u0, double u1)
            {
                U0 = u0;
                U1 = u1;
            }
            public double U1 { get; }
            public double U0 { get; }
        }
    }
}