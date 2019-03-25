using System;
using System.Collections;
using System.Collections.Generic;

namespace FEM
{
    /// <summary>
    /// Базисні функції Bi+1/2 для апроксимації похибки
    /// </summary>
    public class B_BasisFunc
    {
        public static double B_i(double x, int i, IList<double> xi)
        {
            double ret = 0.0;
            if (x > xi[i] && x < xi[i + 1])
            {
                double h = xi[i + 1] - xi[i];
                ret = -4 * (x - xi[i]) * (x - xi[i + 1]) / (h * h);
            }
            return ret;
        }

        public static double d_dx(double x, int i, IList<double> xi)
        {
            double ret = 0.0;
            if (x > xi[i] && x < xi[i + 1])
            {
                double h = xi[i + 1] - xi[i];
                ret = -4 * (2 * x - xi[i] - xi[i + 1]) / (h * h);
            }
            return ret;
        }
    }

    public class CurantFunction
    {
        public static double Fi(double arg, int i, IList<double> x)
        {
            double Step(int j) => x[j] - x[j - 1];

            if (i == x.Count - 1)
            {
                if ((int)arg == 1)
                {
                    return 1;
                }

                if (arg > x[i - 1] && arg <= x[i])
                {
                    return (arg - x[i - 1]) / Step(i);
                }
                else
                {
                    return 0;
                }
            }

            if (i == 0)
            {
                if (arg >= x[0] && arg <= x[1])
                {
                    return ((x[1] - arg) / Step(1));
                }
                else
                {
                    return 0;
                }
            }

            if (x[i - 1] < arg && arg <= x[i])
            {
                return (arg - x[i - 1]) / Step(i);
            }

            if (x[i] < arg && arg <= x[i + 1])
            {
                return (x[i + 1] - arg) / Step(i + 1);
            }

            else
            {
                return 0;
            }

            throw new ArgumentOutOfRangeException(nameof(arg));
        }

        public static double Derivative(double arg, int i, IList<double> x)
        {
            double Step(int j) => x[j] - x[j - 1];

            if (i == 0)
            {
                if (arg >= x[0] && arg <= x[1])
                {
                    return -1 / Step(1);
                }
                else
                {
                    return 0;
                }
            }

            if (i == x.Count - 1)
            {
                if (arg > x[i - 1] && arg <= x[i])
                {
                    return 1 / Step(i);
                }
                else
                {
                    return 0;
                }
            }

            if (x[i - 1] < arg && arg <= x[i])
            {
                return 1 / Step(i);
            }

            if (x[i] < arg && arg <= x[i + 1])
            {
                return -1 / Step(i + 1);
            }
            else
            {
                return 0;
            }

            throw new ArgumentOutOfRangeException(nameof(arg));
        }
    }
}
