using System;
using static System.Collections.Generic.SeriesCreate;
using static System.Collections.Generic.Create;
using static System.Math;
using SystemOfEquations;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace FEM
{
    public partial class DiffusionConvectionReaction : FEMBase
    {
        public Func<double, double> Mu;
        public Func<double, double> Beta;
        public Func<double, double> Sigma;
        public double Alpha;
        public IList<double> Eta;
        public Func<double, double> F;
        public BoundaryCondition Condition;

        /// <summary>
        /// К-сть інтегралів.
        /// </summary>
        private double nIntegrals;
        private double[] u;
        private double[][] K;
        private List<double> e_coefficients;
        private double lastAddedElements;
        private int N_prev;
        private double E_prev;

        public DiffusionConvectionReaction(
            double[] elements,
            BoundaryCondition condition,
            Func<double, double> f,
            Func<double, double> mu,
            Func<double, double> beta,
            Func<double, double> sigma,
            double alpha
          ) : base(elements)
        {
            N_prev = elements.Length - 1;
            Condition = condition;
            F = f;
            Mu = mu;
            Beta = beta;
            Sigma = sigma;
            Alpha = alpha;
            nIntegrals = N - 1;
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

        double Uh_VNorm_Sqr()
        {
            double sum = 0.0;
            if (N == 3)
            {
                sum = u[0] * (u[0] * K[1][0]);
            }
            else
            {
                sum += u[0] * ((u[0] * K[1][0]) + (u[1] * K[2][0]));
                for (int k = 1; k < N - 3; k++)
                {
                    sum += u[k] * ((u[k - 1] * K[0][k - 1]) + (u[k] * K[1][k]) + (u[k + 1] * K[2][k]));
                }
                sum += u[N - 3] * ((u[N - 4] * K[0][N - 4]) + (u[N - 3] * K[1][N - 3]));
            }
            return sum;
        }

        public async Task StartAdaptationAlgorithm(double allowedErrorInPercents, Func<List<string>, Task> progress)
        {
            int iter = 0;
            while (true)
            {
                iter++;
                int k = 0;
                int i1 = 0;
                int h_adaptation_count = 0;
                while (k < Eta.Count)
                {
                    if (Eta[k] > allowedErrorInPercents * 1.61)
                    {
                        //h_adaptation_count++;
                        //double x = (Elements[i1 + 1] + Elements[i1]) / 2.0;
                        //InsertFiniteElement(x, i1);
                        //i1++;
                        int numInsert = (int)Floor((Eta[k] / allowedErrorInPercents) / 5.0);
                        if (numInsert == 0)
                        {
                            numInsert = 1;
                        }

                        double step = (Elements[i1 + 1] - Elements[i1]) / (numInsert + 1);

                        for (int p = 0; p < numInsert; p++)
                        {
                            h_adaptation_count++;
                            InsertFiniteElement(Elements[i1] + step, i1);
                            i1++;
                        }
                    }
                    k++;
                    i1++;
                }
                Console.WriteLine("Iteration: " + iter);
                lastAddedElements = h_adaptation_count;
                if (h_adaptation_count == 0)
                {
                    break;
                }

                Solve();
                Calc_Eh();
                List<string> outputData = CalcEta();
                nIntegrals += N - 1;

                await progress(outputData);
            }
        }

        public List<string> CalcEta()
        {
            Eta = new List<double>(this.N - 1);

            double errorVNorm = ErrorVNorm();
            double uhNorm = Uh_VNorm_Sqr();
            double denom = Sqrt(uhNorm + errorVNorm);
            double maxerror = 0.0;

            for (int i = 0; i < this.N - 1; i++)
            {
                double localError = Sqrt((N - 1) * ErrorVNorm_at_i_element(i)) * 100.0 / denom;
                if (maxerror < localError)
                {
                    maxerror = localError;
                }
                Eta.Add(localError);
            }

            string convergence;
            if (N_prev == 0)
            {
                convergence = "-";
                N_prev = N - 1;
                E_prev = Sqrt(errorVNorm);
            }
            else
            {
                double e_norm = Sqrt(errorVNorm);
                double convergence1 = (Log(E_prev) - Log(e_norm)) / (Log(N - 1) - Log(N_prev));
                E_prev = e_norm;
                N_prev = N - 1;
                convergence = string.Format("{0:0.0000}", convergence1);
            }
            return new List<string>()
            {
                (N-1).ToString(),
                lastAddedElements.ToString(),
                Sqrt(uhNorm).ToString(),
                Sqrt(errorVNorm).ToString(),
                string.Format("{0:0.0000}%", maxerror),
                convergence,
                (9* nIntegrals).ToString()
            };
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
            return (fp12(a, b) * (b - a) / 2.0) + (fp12(b, c) * (c - b) / 2.0);
        }

        private double F_x12(Func<double, double> f, int i)
        {
            return f((Elements[i + 1] - Elements[i]) / 2.0);
        }

        private double Mu12(int i) => F_x12(Mu, i);
        private double Beta12(int i) => F_x12(Mu, i);
        private double Sigma12(int i) => F_x12(Mu, i);

        public double BillinearForm(int i, int j)
        {
            return Mu12(i) * CourantFunction.integrate_dx_dx(i, j, Elements) +
                       Beta(i) * CourantFunction.integrate_dx_fx(i, j, Elements) +
                       Sigma(i) * CourantFunction.integate_fx_fx(i, j, Elements) +
                       Alpha * Fi(1, i) * Fi(1, j);
        }

        public double Calc_a_Bi_Bj(int i)
        {
            return (Mu(i) * B_BasisFunc.integrate_d_dx_2(i, Elements))
                + (Beta(i) * B_BasisFunc.integrate_d_dx_fx(i, Elements))
                + (Sigma(i) * B_BasisFunc.integrate_fx_pow2(i, Elements))
                + (Alpha * Pow(B_BasisFunc.B_i(1, i, Elements), 2));
        }

        /// <summary>
        /// Функціонал <p(Uh),v>
        /// </summary>
        public double Calc_RO_Uh(int j)
        {
            double h = Elements[j + 1] - Elements[j];
            double xj12 = (Elements[j + 1] + Elements[j]) / 2.0;
            var l_bj = ((2.0 * h * F(xj12)) / 3.0) +
                       (Alpha * Condition.U_);


            /////Білінійна форма c(Fi_i,B_i+12)
            double c(int i, int k)
            {
                double ret = 0.0;
                if (i == k - 1)
                {
                    ret = ((Beta(i) * 2.0) / 3.0) + (Sigma(i) * (Elements[i + 1] - Elements[i]) / 3.0);
                }
                if (i == k)
                {
                    ret = ((-Beta(i) * 2.0) / 3.0) + (Sigma(i) * (Elements[i + 1] - Elements[i]) / 3.0);
                }
                return ret;
            }

            double a_uh_v = 0.0;
            if (j == 0)
            {
                a_uh_v = (u[0] * c(0, 1)) + (Alpha * Condition.U1 * B_BasisFunc.B_i(1, j, Elements));
            }
            else if (j == N - 2)
            {
                a_uh_v = (u[j - 1] * c(j, j)) + (Alpha * Condition.U1 * B_BasisFunc.B_i(1, j, Elements));
            }
            else
            {
                a_uh_v = (u[j - 1] * c(j, j)) + (u[j] * c(j, j + 1)) + (Alpha * Condition.U1 * B_BasisFunc.B_i(1, j, Elements));
            }

            return l_bj - a_uh_v;
        }
    }
}