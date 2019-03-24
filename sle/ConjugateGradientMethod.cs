using System;
using System.Collections.Generic;

namespace SystemOfEquations
{
    public class ConjugateGradientMethod : IterationMethod
    {
        public ConjugateGradientMethod(
            double[][] leftPart,
            double[] rightPart,
            double[] initial,
            double epsilon,
            bool isParallel)
            : base(leftPart, rightPart, epsilon, isParallel)
        {
            Initial = initial;
        }

        public double[] Initial { get; set; }

        protected override void TrySolve(double[][] leftPart, IReadOnlyList<double> rightPart)
        {
            double alpha0 = 0, pq = 0;
            double beta0 = 0;
            double rho0 = 0, rho1 = 0;
            var p0 = new double[N];
            var p1 = new double[N];
            var q0 = new double[N];
            var q1 = new double[N];
            var r0 = new double[N];
            var r1 = new double[N];
            var x1 = new double[N];
            var z0 = new double[N];
            var z1 = new double[N];
            var k = 1;

            double residual = 0;

            for (var i = 0; i < N; i++)
            {
                double sum = 0;

                for (var j = 0; j < N; j++)
                    sum += LeftPart[i][j]*Initial[j];

                r0[i] = sum - RightPart[i];
                residual += r0[i]*r0[i];
            }

            residual = Math.Sqrt(residual/N);

            if (residual < Epsilon)
                Answer = Initial;

            while (true)
            {
                // choose the point-preconditioner 

                for (var i = 0; i < N; i++)
                    z0[i] = r0[i]/LeftPart[i][i];

                rho1 = 0;

                for (var i = 0; i < N; i++)
                    rho1 += r0[i]*z0[i];

                if (k == 1)
                {
                    for (var i = 0; i < N; i++)
                        p0[i] = -z0[i];

                    rho0 = rho1;
                }

                else
                {
                    beta0 = rho1/rho0;

                    for (var i = 0; i < N; i++)
                    {
                        p1[i] = -z0[i] - beta0*p0[i];
                        p0[i] = p1[i];
                    }

                    rho0 = rho1;
                }

                for (var i = 0; i < N; i++)
                {
                    double sum = 0;

                    for (var j = 0; j < N; j++)
                        sum += LeftPart[i][j]*p0[j];

                    q1[i] = sum;
                }

                pq = 0;

                for (var i = 0; i < N; i++)
                    pq += p0[i]*q1[i];

                alpha0 = rho0/pq;

                if (k == 1)
                {
                    for (var i = 0; i < N; i++)
                        x1[i] = Initial[i] + alpha0*p0[i];
                }

                else
                {
                    for (var i = 0; i < N; i++)
                        x1[i] = Initial[i] + alpha0*p1[i];
                }

                for (var i = 0; i < N; i++)
                    r1[i] = r0[i] + alpha0*q1[i];

                residual = 0;

                for (var i = 0; i < N; i++)
                {
                    double sum = 0;

                    for (var j = 0; j < N; j++)
                        sum += LeftPart[i][j]*x1[j];

                    r0[i] = sum - RightPart[i];
                    Initial[i] = x1[i];
                    residual += Math.Pow(r0[i], 2);
                }

                residual /= N;

                if (residual < Epsilon
                    || double.IsNaN(Initial[0]))
                    break;

                k++;
            }

            Answer = Initial;
        }

        protected override void TrySolveParallel(double[][] leftPart, double[] rightPart)
        {
            throw new NotImplementedException();
        }
    }
}