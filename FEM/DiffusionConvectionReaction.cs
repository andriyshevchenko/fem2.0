using System;
using System.Linq;
using static System.Collections.Generic.SeriesCreate;
using static System.Collections.Generic.Create;
using static System.Functional.Func;
using static System.Math;
using SystemOfEquations;
using MathNet.Numerics.Integration;

namespace FEM
{
    public partial class DiffusionConvectionReaction : FEMBase
    {
        public double Mu;
        public double Beta;
        public double Omega;
        public double Alpha;
        public Func<double, double> F;
        public BoundaryCondition Condition;

        private double[] u;
        private double[] e_coefficients;

        public DiffusionConvectionReaction(
            double[] elements,
            BoundaryCondition condition,
            Func<double, double> f,
            double mu,
            double beta,
            double omega,
            double alpha
          ) : base(elements)
        {
            Condition = condition;
            F = f;
            Mu = mu;
            Beta = beta;
            Omega = omega;
            Alpha = alpha;
        }

        public double ErrorEstimate(double x)
        {
            double sum = 0.0;
            for (int i = 0; i < e_coefficients.Length; i++)
            {
                sum += e_coefficients[i] * B_BasisFunc.B_i(x, i, Elements);
            }
            return sum;
        }

        public void Calc_Eh()
        {
            u = u ?? Solve();
            e_coefficients = new double[N-1];
            for (int i = 0; i < e_coefficients.Length; i++)
            {
                e_coefficients[i] = Calc_RO_Uh(i) / Calc_a_Bi_Bj(i);
            } 
        }

        public double U(double x)
        {
            u = u ?? Solve();
             
            return CalculateFromBasis(x, u) + Condition.U0; //враховано заміну
        }

        public double[] Solve()
        {
            var method = new TridiagonalSolve(FillMatrix(), FillRightPartVector(), false);
            method.Solve();
            return method.Answer;
        }

        double[][] FillMatrix()
        {
            return array(Series(2, N - 2, i => BillinearForm(i, i - 1)),
                         Series(1, N - 1, i => BillinearForm(i, i)),
                         Series(1, N - 2, i => BillinearForm(i, i + 1)));
        }

        public double[] FillRightPartVector()
        {
            return Series(1, N - 1, L);
        }


        public double L(int i)
        {
            var (a, b) = GetIntegrationBoundsForLinearFunctional(i);
            return GaussLegendreRule.Integrate(x => (F(x) - Condition.U0 * Omega) * Fi(x, i), a, b, 5) + Alpha * Fi(1, i) * Condition.U1;
        }

        public double BillinearForm(int i, int j)
        {
            var func = fun((double x) => Mu * FiDx(i, x) * FiDx(j, x)
                                      + (Beta * FiDx(i, x) + Omega * Fi(x, i))
                                      * Fi(x, j));

            var (a, b) = GetIntegrationBounds(i, j);

            return GaussLegendreRule.Integrate(func, a, b, 5) + Alpha * Fi(b, i) * Fi(b, j);
        }

        public double Calc_a_Bi_Bj(int i)
        {
            var func = fun((double x) => Mu * B_BasisFunc.d_dx(x, i, Elements) * B_BasisFunc.d_dx(x, i, Elements)
                                      + (Beta * B_BasisFunc.d_dx(x, i, Elements) + Omega * B_BasisFunc.B_i(x, i, Elements))
                                      * B_BasisFunc.B_i(x, i, Elements));

            var (a, b) = (Elements[i], Elements[i+1]);

            return GaussLegendreRule.Integrate(func, a, b, 5) + Alpha * B_BasisFunc.B_i(b, i, Elements) * B_BasisFunc.B_i(b, i, Elements);
        }

        /// <summary>
        /// Функціонал <p(Uh),v>
        /// </summary>
        public double Calc_RO_Uh(int j)
        {
            var l_bj = GaussLegendreRule.Integrate(
                          x => (F(x) - Condition.U0 * Omega) * B_BasisFunc.B_i(x, j, Elements), Elements[j], Elements[j+1], 5
                       ) +
                       Alpha * B_BasisFunc.B_i(1, j, Elements) * Condition.U1;

            double a_uh_v = 0.0;
            for (int ind = 0; ind < u.Length; ind++)
            {
                double item = u[ind];
                double a1 = Elements[j];
                double b1 = Elements[j+1];
                a_uh_v += item * GaussLegendreRule.Integrate(
                   x => Mu * FiDx(ind, x) * B_BasisFunc.d_dx(x, j, Elements) + (Beta * FiDx(ind, x) + Omega * Fi(x, ind))
                                 * B_BasisFunc.B_i(x, j, Elements), a1, b1, 5
                ) + 
                Alpha * B_BasisFunc.B_i(b1, ind, Elements) * B_BasisFunc.B_i(b1, j, Elements);
            }
            return l_bj - a_uh_v;
        }
    }
}