using System;
using System.Linq;
using static System.Collections.Generic.SeriesCreate;
using static System.Collections.Generic.Create;
using static System.Functional.Func;
using static System.Math;
using SystemOfEquations;
using MathNet.Numerics.Integration;
using System.Collections.Generic;

namespace FEM
{
    public partial class DiffusionConvectionReaction : FEMBase
    {
        public double Mu;
        public double Beta;
        public double Sigma;
        public double Alpha;
        public Func<double, double> F;
        public BoundaryCondition Condition;

        private double[] u;
        private double[][] K;
        private List<double> e_coefficients;

        public DiffusionConvectionReaction(
            double[] elements,
            BoundaryCondition condition,
            Func<double, double> f,
            double mu,
            double beta,
            double sigma,
            double alpha
          ) : base(elements)
        {
            Condition = condition;
            F = f;
            Mu = mu;
            Beta = beta;
            Sigma = sigma;
            Alpha = alpha;
        }

        double ErrorVNorm()
        {
            double sum = 0.0;

            for (int k = 0; k < N - 1; k++)
            {
                sum += ErrorVNorm_at_i_element(k);
            }

            return sum;
        }

        double ErrorVNorm_at_i_element(int i)
        {
            return Pow(e_coefficients[i], 2) * Calc_a_Bi_Bj(i);
        }

        double Uh_VNorm()
        {
            double sum = 0.0;
            sum += u[0] * (u[0] * K[1][0] + u[1] * K[2][0]);
            for (int k = 1; k < N - 2; k++)
            {
                sum += u[k] * (u[k - 1] * K[0][k - 1] + u[k] * K[1][k] + u[k + 1] * K[2][k]);
            }
            sum += u[N - 2] * (u[N - 3] * K[0][N - 3] + u[N - 2] * K[1][N - 2]);
            return sum; 
        }

        public List<double> StartAdaptationAlgorithm(double allowedErrorInPercents)
        {
            List<double> eta;
            while (true)
            { 
                eta = new List<double>(this.N - 1);
                double errorVNorm = ErrorVNorm();
                double denom = Sqrt(Uh_VNorm() + errorVNorm);
                for (int i = 0; i < this.N - 1; i++)
                {
                    double localError = Sqrt(this.N) * ErrorVNorm_at_i_element(i) * 100.0 / denom;
                    eta.Add(localError);
                }

                int k = 0;
                int i1 = 0;
                int h_adaptation_count = 0;
                while (k < eta.Count)
                {
                    if (eta[k] > allowedErrorInPercents)
                    {
                        h_adaptation_count++;
                        double x = (Elements[i1 + 1] + Elements[i1]) / 2.0;
                        InsertFiniteElement(x, i1);
                        i1++;
                    }
                    k++;
                    i1++;
                }

                if (h_adaptation_count == 0)
                {
                    break;
                }

                Solve();
                Calc_Eh();
            }
            return eta;
        }

        public double Error(double x)
        {
            double sum = 0.0;
            for (int i = 0; i < e_coefficients.Count; i++)
            {
                sum += e_coefficients[i] * B_BasisFunc.B_i(x, i, Elements);
            }
            return sum;
        }

        void InsertFiniteElement(double x, int order)
        {
            Elements.Insert(order + 1, x);
        }

        public void Calc_Eh()
        {
            e_coefficients = new List<double>(N - 1);
            for (int i = 0; i < N - 1; i++)
            {
                e_coefficients.Add(Calc_RO_Uh(i) / Calc_a_Bi_Bj(i));
            }
        }

        public double U(double x)
        {
            return CalculateFromBasis(x, u);  
        }

        public void Solve()
        {
            K = FillMatrix();
            var method = new TridiagonalSolve(K, FillRightPartVector(), false);
            method.Solve();
            u = method.Answer;
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
            return SimpsonRule.IntegrateComposite(x => F(x) * Fi(x, i), a, b, 100) + Alpha * Fi(1, i) * Condition.U_;
        }

        public double BillinearForm(int i, int j)
        { 
            var (a, b) = GetIntegrationBounds(i, j);

            return SimpsonRule.IntegrateComposite(
                x => Mu * FiDx(i, x) * FiDx(j, x) +
                    (Beta * FiDx(i, x) + Sigma * Fi(x, i)) * Fi(x, j),
                a, b, 100
            ) + Alpha * Fi(1,i) * Fi(1, j);
        }

        public double Calc_a_Bi_Bj(int i)
        { 
            return Mu * B_BasisFunc.integrate_d_dx_2(i, Elements)
                + Beta * B_BasisFunc.integrate_d_dx_fx(i, Elements)
                + Sigma * B_BasisFunc.integrate_fx_fx(i, Elements);
        }

        /// <summary>
        /// Функціонал <p(Uh),v>
        /// </summary>
        public double Calc_RO_Uh(int j)
        {
            var l_bj = GaussLegendreRule.Integrate(
                          x => F(x) * B_BasisFunc.B_i(x, j, Elements), Elements[j], Elements[j + 1], 2
                       ) +
                       Alpha * B_BasisFunc.B_i(1, j, Elements) * Condition.U_;

            double a_uh_v = 0.0;
            //for (int ind = 0; ind < u.Length; ind++)
            //{
            //    double item = u[ind];
            //    double a1 = Elements[j];
            //    double b1 = Elements[j + 1];
            //    a_uh_v += item * GaussLegendreRule.Integrate(
            //       x => Mu * FiDx(ind, x) * B_BasisFunc.d_dx(x, j, Elements) + (Beta * FiDx(ind, x) + Omega * Fi(x, ind))
            //                     * B_BasisFunc.B_i(x, j, Elements), a1, b1, 2
            //    ) +
            //    Alpha * B_BasisFunc.B_i(b1, ind, Elements) * B_BasisFunc.B_i(b1, j, Elements);
            //}
            a_uh_v = GaussLegendreRule.Integrate(
                x => Mu*CalculateFromBasisDx(x,u)*B_BasisFunc.d_dx(x,j,Elements) +
                     (Beta*CalculateFromBasisDx(x,u) + Sigma*U(x)) * B_BasisFunc.B_i(x, j, Elements),
                Elements[j], Elements[j+1], 2
            ) + Alpha * U(1) * B_BasisFunc.B_i(1, j, Elements);
            return l_bj - a_uh_v;
        }
    }
}