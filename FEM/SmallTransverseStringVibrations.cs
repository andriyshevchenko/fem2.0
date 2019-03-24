using MathNet.Numerics.Integration;
using System;
using System.Linq;
using static System.Collections.Generic.Create;
using static System.Collections.Generic.SeriesCreate;

namespace FEM
{
    public partial class SmallTransverseStringVibrations : FEMTwoDimBase
    {
        public double W1;
        public double W2;

        double Ro;
        double T0;
        public Func<double, double, double> P;
        public BoundaryCondition Condition;

        double[] V0;

        public SmallTransverseStringVibrations(double tMax, int nT, double[] elements, BoundaryCondition condition, double ro, double t0)
            : base(tMax, nT, elements)
        {
            Ro = ro;
            T0 = t0;
            Condition = condition;

            U0 = partEnd(Elements, 1).Select(Condition.U0).ToArray();
            V0 = partEnd(Elements, 1).Select(Condition.V0).ToArray();
        }

        public Func<double, int, double> Calculate()
        {
            double[] uPrev = U0;
            double[] vPrev = V0;

            for (int i = 0; i < Nt; i++)
            {
                double[] uii = Solve(i, uPrev, vPrev);
                double step = Step(i);
                UiCalculated[i] = Sum(uPrev,
                                      mul((step * step) / 2, uii),
                                      mul(step, vPrev));
                vPrev = Sum(mul(step, uii), vPrev);
                uPrev = UiCalculated[i];
            }

            return (x, tIndex) =>
            {
                return CalculateFromBasis(x, UiCalculated[tIndex]);
            };
        }

        double[][] values;
        double[] Solve(int tIndex, double[] u, double[] v)
        {
            var method = new SystemOfEquations.TridiagonalSolve(values ?? (values = FillMatrix()), FillRightPartVector(tIndex, u, v));
            method.Solve();
            return method.Answer;
        }

        double BillinearForm(int k, int s)
        {
            (double a, double b) = GetIntegrationBounds(k, s);
            return -T0 * GaussLegendreRule.Integrate(x => FiDx(k, x) * FiDx(s, x), a, b, 5);
        }

        public double[] FillRightPartVector(int tIndex, double[] u, double[] v)
        {
            return Series(1, N - 2, i => L(i, tIndex, u, v));
        }

        double[][] FillMatrix()
        {
            var step = Step(0);
            double step2 = W2 * ((step * step) / 2);
            return array(Series(2, N - 3, i => M(i, i - 1) - step2 * BillinearForm(i, i - 1)),
                         Series(1, N - 2, i => M(i, i) - step2 * BillinearForm(i, i)),
                         Series(1, N - 3, i => M(i, i + 1) - step2 * BillinearForm(i, i + 1)));
        }

        double M(int k, int s)
        {
            return Ro * Product(k, s);
        }

        double L(int s, int tIndex, double[] ui, double[] vi)
        {
            var t = GetApproximateTime(tIndex);
            var step = Step(tIndex);
            var (a, b) = GetIntegrationBoundsForLinearFunctional(s);

            double l = GaussLegendreRule.Integrate(x => P(x, t) * Fi(x, s), a, b, 5);

            double sum = 0;
            for (int j = 1; j < N - 1; j++)
            {
                double bf = BillinearForm(j, s);
                sum += bf * ui[j - 1];
                sum += step * W1 * bf * vi[j - 1];
            }

            return l + sum;
        }
    }
}
