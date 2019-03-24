using MathNet.Numerics.Integration;
using System;
using static System.Math;
namespace FEM
{
    public class FEMBase
    {
        public double[] Elements;
        protected int N => Elements.Length;

        public FEMBase(double[] elements)
        {
            Elements = elements;
        }

        protected double CalculateFromBasis(double x, double[] basis)
        {
            double sum = 0.0;
            for (int i = 0; i < basis.Length; i++)
            {
                sum += basis[i] * Fi(x, i+1);
            }
            return sum;
        }

        protected double Fi(double x, int i)
        {
            return CurantFunction.Fi(x, i, Elements);
        }

        protected double FiDx(int i, double x)
        {
            return CurantFunction.Derivative(x, i, Elements);
        }
        
        protected double Product(int k, int s)
        {
            (double a, double b) = GetIntegrationBounds(k, s);
            return GaussLegendreRule.Integrate(x => Fi(x, k) * Fi(x, s), a, b, 5);
        }

        protected double L2Product(Func<double, double> f1, Func<double, double> f2)
        {
            return GaussLegendreRule.Integrate(x => f1(x) * f2(x), 0, 1, 5);
        }

        static double L2Norm(Func<double, double> f)
        {
            return Sqrt(GaussLegendreRule.Integrate(x => Pow(Abs(f(x)), 2), 0, 1, 5));
        }

        public static double ErrorEstimate(Func<double, double> f1, Func<double, double> f2)
        {
            return L2Norm(x => f1(x) - f2(x));
        }

        protected (double, double) GetIntegrationBoundsForLinearFunctional(int i)
        {
            double a = Elements[i - 1];
            double b;
            if (i == N - 1)
            {
                b = Elements[i];
            }
            else
            {
                b = Elements[i + 1];
            }

            return (a, b);
        }

        protected (double, double) GetIntegrationBounds(int i, int j)
        {
            double a;
            double b;

            if (i > j)
            {
                a = Elements[j];
                b = Elements[i];
            }
            else if (i < j)
            {
                a = Elements[i];
                b = Elements[j];
            }
            else if (i == j)
            {
                a = Elements[i - 1];
                b = Elements[i == N - 1 ? i : i + 1];
            }
            else throw new ArgumentOutOfRangeException(nameof(i), nameof(j));
            return (a, b);
        }
    }
}