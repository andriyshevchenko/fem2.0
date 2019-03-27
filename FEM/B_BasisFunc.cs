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
