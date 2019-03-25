using System; 
using static System.Collections.Generic.SeriesCreate;
using static System.Collections.Generic.Create;
using static System.Functional.Func;
using SystemOfEquations;
using System.Linq;
using MathNet.Numerics.Integration;

namespace FEM
{
    public partial class DiffusionConvectionReactionTwoDim : FEMTwoDimBase
    {

        public DiffusionConvectionReactionTwoDim(double tMax, int nT, double[] elements, BoundaryCondition condition)
            : base(tMax, nT, elements)
        {
            Condition = condition;
            U0 = array(partEnd(Elements.ToArray(), 1).Select(Condition.U_0));
        }

        public Func<double, double, double> Mu;
        public Func<double, double, double> Beta;
        public Func<double, double, double> Omega;
        public Func<double, double, double> F;
        public BoundaryCondition Condition;
        const double Theta = 0.5;

        public Func<double, int, double> Calculate()
        {
            double[] uPrev = U0;
            for (int i = 0; i < Nt; i++)
            {
                UiCalculated[i] = Sum(uPrev, mul(Step(i), Solve(i, uPrev)));
                uPrev = UiCalculated[i];
            }

            return (x, tIndex) =>
            {
                return CalculateFromBasis(x, UiCalculated[tIndex]);
            };
        }

        double[] Solve(int tIndex, double[] u)
        {
            var method = new TridiagonalSolve(FillMatrix(tIndex), FillRightPartVector(tIndex, u));
            method.Solve();
            return method.Answer;
        }

        double[][] FillMatrix(int tIndex)
        {
            var ti = GetApproximateTime(tIndex);
            var step = Step(tIndex);
            return array(Series(2, N - 2, i => Product(i, i - 1) + Theta * step * BillinearForm(i, i - 1, ti)),
                         Series(1, N - 1, i => Product(i, i) + Theta * step * BillinearForm(i, i, ti)),
                         Series(1, N - 2, i => Product(i, i + 1) + Theta * step * BillinearForm(i, i + 1, ti)));
        }

        public double[] FillRightPartVector(int tIndex, double[] u)
        {
            return Series(1, N - 1, i => L(i, tIndex, u));
        }

        public double L(int s, int tIndex, double[] ui)
        {
            var t = GetApproximateTime(tIndex);
            var (a, b) = GetIntegrationBoundsForLinearFunctional(s);

            double l = GaussLegendreRule.Integrate(x => (F(x, t) - Condition.U0(t) * Omega(x, t)) * Fi(x, s), a, b, 5)
                                    + Mu(1, t) * Fi(1, s) * Condition.U1(t);

            double sum = 0.0;
            for (int j = 1; j < N; j++)
            {
                sum += ui[j - 1] * BillinearForm(j, s, t);
            }
            return l - sum;
        }

        public double BillinearForm(int k, int s, double t)
        {
            var integrand = fun((double x) => Mu(x, t) * FiDx(k, x) * FiDx(s, x)
                                      + (Beta(x, t) * FiDx(k, x) + Omega(x, t) * Fi(x, k))
                                      * Fi(x, s));

            (double a, double b) = GetIntegrationBounds(k, s);
            return GaussLegendreRule.Integrate(integrand, a, b, 5);
        }
    }
}
