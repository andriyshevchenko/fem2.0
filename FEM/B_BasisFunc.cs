using System.Collections.Generic;
using static System.Math;

namespace FEM
{
    /// <summary>
    /// Базисні функції Bi+1/2 для апроксимації похибки
    /// </summary>
    public class B_BasisFunc
    {
        public static double integrate_fx_fx(int i, IList<double> x)
        {
            double h = x[i + 1] - x[i];
            double a = x[i] + x[i + 1];
            double c = x[i] * x[i + 1];
            double antiderivative(double val)
            {
                var x2 = val * val;
                return c * c * val + x2 * (-a * c + (a * a + 2 * c) * val / 3.0 - a * x2 / 2.0 + x2 * val / 5.0);
            }
            return 16.0 * (antiderivative(x[i + 1]) - antiderivative(x[i])) / Pow(h, 4);
        }

        public static double integrate_d_dx_fx(int i, IList<double> x)
        {
            double h = x[i + 1] - x[i];
            double a = x[i] + x[i + 1];
            double c = x[i] * x[i + 1];
            double antiderivative(double val)
            {
                var multiplier1 = val * (val - a);
                return multiplier1 * (multiplier1 + 2 * c);
            }
            return 8.0 * (antiderivative(x[i + 1]) - antiderivative(x[i])) / Pow(h, 4);
        }

        public static double integrate_d_dx_2(int i, IList<double> x)
        {
            double a = x[i];
            double b = x[i + 1];
            double h = b - a;
            double a1 = x[i] + x[i + 1];
            double antiderivative(double val) =>
                4.0 * Pow(val, 3) / 3.0 - val * val * 2 * a1 + a1 * a1 * val;
            return 16.0 * (antiderivative(b) - antiderivative(a)) / Pow(h, 4);
        }

        public static double B_i(double x, int i, IList<double> xi)
        {
            double ret = 0.0;
            if (x >= xi[i] && x <= xi[i + 1])
            {
                double h = xi[i + 1] - xi[i];
                ret = -4 * (x - xi[i]) * (x - xi[i + 1]) / (h * h);
            }
            return ret;
        }

        public static double d_dx(double x, int i, IList<double> xi)
        {
            double ret = 0.0;
            if (x >= xi[i] && x <= xi[i + 1])
            {
                double h = xi[i + 1] - xi[i];
                ret = -4 * (2 * x - xi[i] - xi[i + 1]) / (h * h);
            }
            return ret;
        }
    }
}
