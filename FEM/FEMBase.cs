using System;
using System.Collections.Generic;
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
    }
}