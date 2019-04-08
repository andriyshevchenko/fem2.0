using System;
using System.Collections.Generic;
using static System.Math;

namespace FEM
{
    public class CourantFunction
    {
        public static double integrate_dx_dx(int i, int j, IList<double> x)
        {
            double ret = 0.0;
            if (i == j)
            {
                ret = 1.0 / (x[i] - x[i - 1]);
                if (i + 1 < x.Count)
                {
                    ret += 1.0 / (x[i + 1] - x[i]);
                }
            }
            else if (i == j - 1)
            {
                ret = -1.0 / (x[i + 1] - x[i]);
            }
            else if (i == j + 1)
            {
                ret = -1.0 / (x[i] - x[i - 1]);
            }
            return ret;
        }

        public static double integrate_dx_fx(int i, int j, IList<double> x)
        {
            double ret = 0.0;

            if (i == j)
            {
                ret = 0.0;
            }
            else if (i == j - 1)
            {
                ret = 0.5;
            }
            else if (i == j + 1)
            {
                ret = -0.5;
            }
            return ret;
        }

        public static double integate_fx_fx(int i, int j, IList<double> x)
        {
            double ret = 0.0;
            double antiderivate(double value, double a) =>
                (a * a * value) - (a * value * value) + (value * value * value / 3.0);
            double antiderivate2(double value, double a, double a1) =>
                (-value * value * value / 3.0) + (value * value * (a + a1) / 2.0) - (a * a1 * value);

            if (i == j)
            {
                ret = (antiderivate(x[i], x[i - 1]) - antiderivate(x[i - 1], x[i - 1])) / Pow(x[i] - x[i - 1], 2);
                if (i + 1 < x.Count)
                {
                    ret += (antiderivate(x[i + 1], x[i + 1]) - antiderivate(x[i], x[i + 1])) / Pow(x[i + 1] - x[i], 2);
                }
            }
            else if (i == j - 1)
            {
                ret = (antiderivate2(x[i + 1], x[i], x[i + 1]) - antiderivate2(x[i], x[i], x[i + 1])) / Pow(x[i + 1] - x[i], 2);
            }
            else if (i == j + 1)
            {
                ret = (antiderivate2(x[i], x[i - 1], x[i]) - antiderivate2(x[i - 1], x[i - 1], x[i])) / Pow(x[i] - x[i - 1], 2);
            }
            return ret;
        }

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
                if (arg > x[0] && arg <= x[1])
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
                if (arg > x[0] && arg <= x[1])
                {
                    return -1.0 / Step(1);
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
                    return 1.0 / Step(i);
                }
                else
                {
                    return 0;
                }
            }

            if (x[i - 1] < arg && arg <= x[i])
            {
                return 1.0 / Step(i);
            }

            if (x[i] < arg && arg <= x[i + 1])
            {
                return -1.0 / Step(i + 1);
            }
            else
            {
                return 0;
            }

            throw new ArgumentOutOfRangeException(nameof(arg));
        }
    }
}
