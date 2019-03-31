using MathNet.Numerics.Integration;
using System;
using System.Collections.Generic;
using static System.Math;
namespace FEM
{
    public class FEMBase
    {
        public List<double> Elements;
        protected int N => Elements.Count;

        public FEMBase(double[] elements)
        {
            Elements = new List<double>(elements);
        }

        protected double CalculateFromBasisDx(double x, double[] basis)
        {
            double sum = 0.0;
            for (int i = 0; i < basis.Length; i++)
            {
                double fi = FiDx(i + 1, x);
                sum += basis[i] * fi;
            }
            return sum;
        }

        protected double CalculateFromBasis(double x, double[] basis)
        {
            double sum = 0.0;
            for (int i = 0; i < basis.Length; i++)
            {
                double fi = Fi(x, i + 1);
                sum += basis[i] * fi;
            }
            return sum;
        }

        protected double Fi(double x, int i)
        {
            return CourantFunction.Fi(x, i, Elements);
        }

        protected double FiDx(int i, double x)
        {
            return CourantFunction.Derivative(x, i, Elements);
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