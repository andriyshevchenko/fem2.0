using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace SystemOfEquations
{
    public class LUDecomposition : SolveMethod
    {
        public LUDecomposition(double[][] leftPart, double[] rightPart, bool isParallel)
            : base(leftPart, rightPart, isParallel)
        {
        }

        protected override void TrySolve(double[][] leftPart, double[] rightPart)
        {
            Initialize(out double[][] u, out double[][] l);

            for (var i = 0; i < N; i++)
            {
                l[i][i] = 1;
            }
            for (var j = 0; j < N; j++)
            {
                u[0][j] = leftPart[0][j];
            }
            for (var j = 1; j < N; j++)
            {
                l[j][0] = leftPart[j][0] / u[0][0];
            }

            for (var i = 1; i < N; i++)
            {
                for (var j = i - 1; j < N; j++)
                {
                    double sum = 0;
                    for (var k = 0; k < i; k++)
                    {
                        sum += l[i][k] * u[k][j];
                    }
                    u[i][j] = leftPart[i][j] - sum;
                }

                for (var j = i; j < N; j++)
                {
                    double sum = 0;
                    for (var k = 0; k < i; k++)
                    {
                        sum += l[j][k] * u[k][i];
                    }
                    l[j][i] = 1 / u[i][i] * (leftPart[j][i] - sum);
                }
            }
            //Ly = B
            var y = ForwardSubstitution(l, rightPart);

            //Ux = y;
            Answer = ReverseSubstitution(u, y);
        }

        private void Initialize(out double[][] u, out double[][] l)
        {
            u = new double[N][];
            l = new double[N][];
            for (var i = 0; i < N; i++)
            {
                u[i] = new double[N];
                l[i] = new double[N];
            }
        }

        private static double[] ReverseSubstitution(IReadOnlyList<double[]> leftPart, double[] rightPart)
        {
            var N = rightPart.GetLength(0);
            var answer = new double[N];
            answer[N - 1] = rightPart[N - 1] / leftPart[N - 1][N - 1];

            for (var i = N - 2; i >= 0; i--)
            {
                double sum = 0;
                for (var j = N - 1; j > i; j--)
                {
                    sum += answer[j] * leftPart[i][j];
                }
                answer[i] = (rightPart[i] - sum) / leftPart[i][i];
            }
            var t = answer;
            return t;
        }

        private static double[] ForwardSubstitution(IReadOnlyList<double[]> leftPart, double[] rightPart)
        {
            var N = rightPart.Length;
            var y = new double[N];
            y[0] = rightPart[0] / leftPart[0][0];

            for (var i = 1; i < N; i++)
            {
                double sum = 0;
                for (var j = 0; j < i; j++)
                {
                    sum += y[j] * leftPart[i][j];
                }
                y[i] = (rightPart[i] - sum) / leftPart[i][i];
            }

            return y;
        }

        protected override void TrySolveParallel(double[][] leftPart, double[] rightPart)
        {
            Initialize(out double[][] u, out double[][] l);

            Parallel.For(0, N, i =>
            {
                l[i][i] = 1;
                u[0][i] = leftPart[0][i];
            });

            Parallel.For(1, N, j =>
            {
                l[j][0] = leftPart[j][0] / u[0][0];
            });

            Algorithm(leftPart, u, l, N);

            //Ly = B
            var y = ForwardSubstitution(l, rightPart);

            //Ux = y;
            Answer = ReverseSubstitution(u, y);
        }

        private static void Algorithm(
                    IReadOnlyList<double[]> leftPart,
                    IReadOnlyList<double[]> u,
                    IReadOnlyList<double[]> l, int N)
        {

            for (var i = 1; i < N; i++)
            {
                var _i = i;

                for (var j = i - 1; j < N; j++)
                {
                    var _j = j;
                    var temp = new double[_i];

                    Parallel.For(0, _i, k =>
                    {
                        temp[k] = l[_i][k]*u[k][_j];
                    });

                    var sum = temp.AsParallel().Sum(); 

                    u[i][j] = leftPart[i][j] - sum;
                }

                for (var j = i; j < N; j++)
                {
                    double sum = 0;

                    var _j = j;

                    var temp = new double[_i];
                    Parallel.For(0, i, k =>
                    {
                        temp[k] = l[_j][k] * u[k][_i];
                    });

                    sum = temp.AsParallel().Sum();

                    l[j][i] = 1 / u[i][i] * (leftPart[j][i] - sum);
                }
            }
        }
    }
}