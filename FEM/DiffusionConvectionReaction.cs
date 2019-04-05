using System;
using static System.Collections.Generic.SeriesCreate;
using static System.Collections.Generic.Create;
using static System.Math;
using SystemOfEquations;
using MathNet.Numerics.Integration;
using System.Collections.Generic;
using System.Diagnostics;
using System.Threading.Tasks;

namespace FEM
{
    public partial class DiffusionConvectionReaction : FEMBase
    {
        public double Mu;
        public double Beta;
        public double Sigma;
        public double Alpha;
        public IList<double> Eta;
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
                sum += Pow(ErrorVNorm_at_i_element(k), 2);
            }

            return sum;
        }

        double ErrorVNorm_at_i_element(int i)
        {
            return Abs(e_coefficients[i]) * Sqrt(Calc_a_Bi_Bj(i));
        }

        double Uh_VNorm_Sqr()
        {
            return SimpsonRule.IntegrateComposite(x =>
                Mu * Pow(Udx(x), 2) + Beta * Udx(x) * U(x) + Sigma * Pow(U(x), 2), 0, 1, 10000
            ) + Alpha * Pow(U(1), 2);
            double sum = 0.0;
            sum += u[0] * (u[0] * K[1][0] + u[1] * K[2][0]);
            for (int k = 1; k < N - 3; k++)
            {
                sum += u[k] * (u[k - 1] * K[0][k - 1] + u[k] * K[1][k] + u[k + 1] * K[2][k]);
            }
            sum += u[N - 2] * (u[N - 3] * K[0][N - 3] + u[N - 2] * K[1][N - 2]);
            return sum;
        }

        public async Task StartAdaptationAlgorithm(double allowedErrorInPercents, Func<Task> progress)
        {
            while (true)
            {
                int k = 0;
                int i1 = 0;
                int h_adaptation_count = 0;
                while (k < Eta.Count)
                {
                    if (Eta[k] > allowedErrorInPercents)
                    {
                        h_adaptation_count++;
                        double x = (Elements[i1 + 1] + Elements[i1]) / 2.0;
                        InsertFiniteElement(x, i1);
                        i1++;
                        //h_adaptation_count++;

                        //int numInsert = (int) Floor(Sqrt(eta[k] / allowedErrorInPercents) / 2.0);
                        //if (numInsert == 0)
                        //{
                        //    numInsert = 1;
                        //}

                        //double step = (Elements[i1 + 1] - Elements[i1]) / (numInsert + 1);

                        //for (int p = 0; p < numInsert; p++)
                        //{
                        //    InsertFiniteElement(Elements[i1] + step, i1);
                        //    i1++;
                        //}
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
                CalcEta();

                await progress();
            }
        }

        public void CalcEta()
        {
            Eta = new List<double>(this.N - 1);
            double errorVNorm = ErrorVNorm();
            double denom = Sqrt(Uh_VNorm_Sqr() + errorVNorm);
            for (int i = 0; i < this.N - 1; i++)
            {
                double localError = Sqrt(this.N - 1) * ErrorVNorm_at_i_element(i) * 100.0 / denom;
                Eta.Add(localError);
            }
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

        public double Udx(double x)
        {
            return CalculateFromBasisDx(x, u);
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
            return array(Series(2, N - 3, i => BillinearForm(i, i - 1)),
                         Series(1, N - 2, i => BillinearForm(i, i)),
                         Series(1, N - 3, i => BillinearForm(i, i + 1)));
        }

        public double[] FillRightPartVector()
        {
            return Series(1, N - 2, L);
        }

        public double L(int i)
        {
            double a = Elements[i - 1];
            double b = Elements[i];
            double c = Elements[i + 1];
            double fp12(double l, double r) => F((l + r) / 2.0);
            return fp12(a, b) * (b - a) / 2.0 + fp12(b, c) * (c - b) / 2.0;
        }

        public double BillinearForm(int i, int j)
        {
            return Mu * CourantFunction.integrate_dx_dx(i, j, Elements) +
                       Beta * CourantFunction.integrate_dx_fx(i, j, Elements) +
                       Sigma * CourantFunction.integate_fx_fx(i, j, Elements) +
                       Alpha * Fi(1, i) * Fi(1, j);
        }

        public double Calc_a_Bi_Bj(int i)
        {
            return Mu * B_BasisFunc.integrate_d_dx_2(i, Elements)
                + Beta * B_BasisFunc.integrate_d_dx_fx(i, Elements)
                + Sigma * B_BasisFunc.integrate_fx_fx(i, Elements)
                + Alpha * Pow(B_BasisFunc.B_i(1, i, Elements), 2);
        }

        /// <summary>
        /// Функціонал <p(Uh),v>
        /// </summary>
        public double Calc_RO_Uh(int j)
        {
            var l_bj = (2.0 * (Elements[j + 1] - Elements[j]) * F((Elements[j + 1] + Elements[j]) / 2.0)) / 3.0 +
                       Alpha * Condition.U_;

            double a_uh_v = 0.0;

            /////Білінійна форма c(Fi_i,B_i+12)
            double c(int i, int k)
            {
                double ret = 0.0;
                if (i == k - 1)
                {
                    ret = (Beta * 2.0) / 3.0 + Sigma * (Elements[i + 1] - Elements[i]) / 3.0;
                }
                if (i == k)
                {
                    ret = (-Beta * 2.0) / 3.0 + Sigma * (Elements[i + 1] - Elements[i]) / 3.0;
                }
                return ret;
            }

            if (j == 0)
            {
                a_uh_v = u[0] * c(0, 1) + Alpha * Condition.U1 * B_BasisFunc.B_i(1, j, Elements);
            }
            else if (j == N - 2)
            {
                a_uh_v = u[j - 1] * c(j, j) + Alpha * Condition.U1 * B_BasisFunc.B_i(1, j, Elements);
            }
            else
            {
                a_uh_v = u[j - 1] * c(j, j) + u[j] * c(j, j + 1) + Alpha * Condition.U1 * B_BasisFunc.B_i(1, j, Elements);
            }

            return l_bj - a_uh_v;
        }
    }
}