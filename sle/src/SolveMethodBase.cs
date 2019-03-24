namespace SystemOfEquations
{
    public abstract class SolveMethod
    {
        protected SolveMethod(double[][] leftPart, double[] rightPart, bool isParallel)
        {
            LeftPart = leftPart;
            RightPart = rightPart;
            IsParallel = isParallel;

            N = RightPart.GetLength(0);
            Answer = new double[N];
        }


        public double[] Answer { get; protected set; }
        protected double[][] LeftPart { get; }
        protected double[] RightPart { get; }
        protected bool IsParallel { get; }
        public int N { get; }
        protected abstract void TrySolve(double[][] leftPart, double[] rightPart);
        protected abstract void TrySolveParallel(double[][] leftPart, double[] rightPart);

        public void Solve()
        {
            if (IsParallel)
            {
                TrySolveParallel(LeftPart, RightPart);
                return;
            }
            TrySolve(LeftPart, RightPart);
        }
    }
}