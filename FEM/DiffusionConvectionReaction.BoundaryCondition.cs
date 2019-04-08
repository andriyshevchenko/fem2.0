namespace FEM
{
    public partial class DiffusionConvectionReaction : FEMBase
    {
        public struct BoundaryCondition
        {
            public BoundaryCondition(double u_, double u1)
            {
                U_ = u_;
                U1 = u1;
            }
            public double U1 { get; }
            public double U_ { get; }
        }
    }
}